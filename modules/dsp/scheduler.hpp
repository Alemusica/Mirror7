#pragma once

#include "aureonoise3beta/core/state.hpp"
#include "aureonoise3beta/dsp/spatial_utils.hpp"
#include "aureo_core/aureo_math.hpp"

#include <cmath>
#include <algorithm>

#include "../harmony/harmony.hpp"

namespace aureonoise {
namespace scheduler {

inline double manual_rate(const t_aureonoise* x)
{
  return aureonoise::control::tempo_manual_rate(x);
}

inline int find_free_grain(t_aureonoise* x)
{
  for (int i = 0; i < aureo::kMaxGrains; ++i) {
    if (!x->grains[i].on) return i;
  }
  return -1;
}

inline int schedule_gap_samples(t_aureonoise* x)
{
  if (x && x->p_dialogue_on) {
    auto& dq = x->dialogue;
    if (dq.fib_gap_queue_pos < dq.fib_gap_queue_size) {
      const double override_gap = std::max(1.0, dq.fib_gap_queue[dq.fib_gap_queue_pos++]);
      if (dq.fib_gap_queue_pos >= dq.fib_gap_queue_size) {
        dq.fib_gap_queue_pos = 0;
        dq.fib_gap_queue_size = 0;
        dq.fib_gap_queue.fill(0.0);
      }
      return static_cast<int>(std::max(1.0, std::round(override_gap)));
    }
  }

  const bool sync_active = aureonoise::control::tempo_sync_active(x);
  double base_rate = sync_active ? aureonoise::control::tempo_current_rate(x) : manual_rate(x);
#if AUREO_THERMO_LATTICE
  if (x->p_thermo && !sync_active) {
    const double ur = 0.5 + 0.5 * std::tanh(x->ou_rate.y);
    const double mapped = aureo::map_phi_range(std::max(0.001, base_rate / aureo::kPhi), base_rate * aureo::kPhi, ur);
    base_rate = aureo::clamp(mapped, 0.0, aureo::kMaxEventRateHz);
  }
#endif
  if (base_rate <= 1e-6) base_rate = manual_rate(x);

  if (sync_active) {
    const double gap_samples = std::max(1.0, x->sr / base_rate);
    return static_cast<int>(std::round(gap_samples));
  }

  const double t = static_cast<double>(x->sample_counter) / x->sr;
  double lambda = base_rate * (1.0 + 0.2 * std::sin(2.0 * aureo::kPi * (t * (1.0 / aureo::kPhi))));
  lambda = std::max(lambda, 1e-3);
#if AUREO_THERMO_LATTICE && AUREO_BURST_HAWKES
  if (x->p_burst) lambda += 0.3 * x->hawkes.lambda;
#endif

  const double U = std::max(1.0e-12, x->rng.uni01());
  double gap_sec = -std::log(U) / lambda;
  const double cap_sec = std::min(30.0, 4.0 / std::max(base_rate, 1.0e-6));
  gap_sec = std::min(gap_sec, cap_sec);
  double gap_samples = std::max(1.0, gap_sec * x->sr);
  gap_samples = quantize_gap_samples(x->harmony, gap_samples);
  return static_cast<int>(std::max(1.0, std::round(gap_samples)));
}

inline void reset_sequences(t_aureonoise* x)
{
  x->rng.seed(static_cast<uint64_t>(x->p_seed));
  x->w_phi.x = x->rng.uni01();
  x->w_phi.set_step(aureo::kInvPhi);
#if AUREO_WD_PHI_POWERS
  x->w_s2.x = x->rng.uni01();
  x->w_s2.set_step(1.0 / (aureo::kPhi * aureo::kPhi));
  x->w_pl.x = x->rng.uni01();
  x->w_pl.set_step(1.0 / (aureo::kPhi * aureo::kPhi * aureo::kPhi));
#else
  x->w_s2.x = x->rng.uni01();
  x->w_s2.set_step(aureo::kInvSqrt2);
  x->w_pl.x = x->rng.uni01();
  x->w_pl.set_step(aureo::kInvPlastic);
#endif
  x->noise.mode = static_cast<aureo::NoiseColor>(static_cast<int>(x->p_noise_mode));
  x->noise.classic.color = static_cast<aureo::NoiseColor>(x->p_noise_color);
  x->noise.classic.amount = aureo::clamp01(x->p_color_amt);
  x->noise.aureo.phi_pi_mix = aureo::clamp01(x->p_aureo_mix);
  x->noise.aureo.planck_decay = aureo::clamp01(x->p_aureo_decay);
  x->noise.aureo.prime_stride = std::max(0.1, x->p_aureo_stride);
  x->noise.aureo.active_partials = static_cast<int>(std::clamp(x->p_aureo_harmonics, 1.0, static_cast<double>(aureo::NoiseAureoState::kMaxPartials)));
  x->noise.aureo.velvet_on = (x->p_aureo_velvet != 0);
  x->noise.reset();
  x->tone.clear();
}

inline uint32_t map_len_samples(t_aureonoise* x, double u)
{
  const double base = aureo::clamp(x->p_baselen_ms, aureo::kMinBaseLengthMs, 2000.0) * 0.001 * x->sr;
  const double kexp = (2.0 * u - 1.0) * aureo::clamp01(x->p_len_phi);
  double L = base * std::pow(aureo::kPhi, kexp);
  L = aureo::clamp(L, aureo::kMinGrainSamples, x->sr * 4.0);
  L = quantize_length_samples(x->harmony, L);
  L = aureo::clamp(L, aureo::kMinGrainSamples, x->sr * 4.0);
  return static_cast<uint32_t>(L);
}

} // namespace scheduler
} // namespace aureonoise
