#include "aureonoise3beta/core/state.hpp"
#include "harmony.hpp"
#include "aureonoise3beta/control/tempo.hpp"

#include <algorithm>

namespace {

double compute_reference_gap(const t_aureonoise* x)
{
  if (!x || x->sr <= 0.0) return 0.0;
  double rate = aureonoise::control::tempo_current_rate(x);
  if (rate <= 1.0e-6) rate = aureonoise::control::tempo_manual_rate(x);
  rate = std::max(0.01, rate);
  return x->sr / rate;
}

double compute_reference_length(const t_aureonoise* x)
{
  if (!x || x->sr <= 0.0) return 0.0;
  const double base_ms = aureo::clamp(x->p_baselen_ms, aureo::kMinBaseLengthMs, 2000.0);
  return std::max(1.0, base_ms * 0.001 * x->sr);
}

} // namespace

void aureonoise_refresh_harmony(t_aureonoise* x)
{
  if (!x) return;
  aureonoise::HarmonyConfig cfg;
  cfg.prime_comb = (x->p_prime_comb != 0);
  cfg.time_quant = (x->p_time_quant != 0);
  cfg.time_strict = (x->p_time_strict != 0);
  cfg.reference_gap = compute_reference_gap(x);
  cfg.reference_length = compute_reference_length(x);
  cfg.active_partials = static_cast<int>(std::clamp(x->p_aureo_harmonics, 1.0, static_cast<double>(aureo::NoiseAureoState::kMaxPartials)));
  x->harmony = aureonoise::make_harmony_state(cfg);
  x->noise.aureo.configure_mask(x->harmony.partial_mask);
  x->noise.aureo.active_partials = std::max(1, x->harmony.active_count);
}
