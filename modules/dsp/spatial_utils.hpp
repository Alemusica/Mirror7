#pragma once

#include "aureonoise3beta/core/state.hpp"
#include "aureo_core/aureo_math.hpp"

#include <algorithm>
#include <cmath>

struct AureonoiseExternalState {
  bool active = false;
  double amount = 0.0;
  double feedback = 0.0;
  double cross = 0.0;
  double hf = 0.0;
  int delay_samples = 1;
};

inline void aureonoise_reset_external_buffers(t_aureonoise* x)
{
  if (!x) return;
  x->ext_delayL.fill(0.0);
  x->ext_delayR.fill(0.0);
  x->ext_write_idx = 0;
  x->ext_lpL = 0.0;
  x->ext_lpR = 0.0;
}

inline AureonoiseExternalState aureonoise_prepare_external_state(const t_aureonoise* x)
{
  AureonoiseExternalState state;
  if (!x) return state;
  state.amount = aureo::clamp01(x->p_spat_external);
  state.active = (state.amount > 1.0e-6);
  if (!state.active) return state;

  const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
  state.delay_samples = std::clamp(static_cast<int>(std::round((0.0007 + 0.0045 * state.amount) * sr)),
                                   1, t_aureonoise::kExternalDelaySize - 1);
  state.feedback = 0.07 * state.amount;
  state.cross = 0.32 * state.amount;
  state.hf = 0.55 * state.amount;
  return state;
}

inline void aureonoise_apply_external_mix(t_aureonoise* x,
                                          const AureonoiseExternalState& st,
                                          double& yL,
                                          double& yR,
                                          uint32_t& write_idx)
{
  if (!st.active || !x) return;
  const uint32_t mask = static_cast<uint32_t>(t_aureonoise::kExternalDelayMask);
  const uint32_t read_idx = (write_idx + static_cast<uint32_t>(t_aureonoise::kExternalDelaySize) -
                             static_cast<uint32_t>(st.delay_samples)) & mask;
  const double prevL = x->ext_delayL[read_idx];
  const double prevR = x->ext_delayR[read_idx];

  x->ext_delayL[write_idx] = yL + st.feedback * prevL;
  x->ext_delayR[write_idx] = yR + st.feedback * prevR;

  x->ext_lpL = prevL + (x->ext_lpL - prevL) * (1.0 - 0.32);
  x->ext_lpR = prevR + (x->ext_lpR - prevR) * (1.0 - 0.32);

  const double hfL = prevL - x->ext_lpL;
  const double hfR = prevR - x->ext_lpR;
  const double extL = prevL + st.hf * hfL + st.cross * prevR;
  const double extR = prevR + st.hf * hfR + st.cross * prevL;

  yL = (1.0 - st.amount) * yL + st.amount * extL;
  yR = (1.0 - st.amount) * yR + st.amount * extR;
  write_idx = (write_idx + 1u) & mask;
}

inline double aureonoise_apply_phi_pan(t_aureonoise* x, double pan)
{
  if (!x || !x->p_phi_pan) return aureo::clamp(pan, -1.0, 1.0);
  if (x->phi_last_mag > 1.0e-6) {
    const int next_sign = (x->phi_last_sign >= 0) ? -1 : 1;
    const double prev_mag = x->phi_last_mag;
    double target_mag = prev_mag * aureo::kPhi;
    if (target_mag > 0.999) {
      const double fallback = prev_mag / aureo::kPhi;
      if (fallback >= 0.05) {
        target_mag = fallback;
      } else {
        target_mag = std::max(0.05, std::abs(pan));
      }
    }
    target_mag = aureo::clamp(target_mag, 0.05, 1.0);
    return aureo::clamp(static_cast<double>(next_sign) * target_mag, -1.0, 1.0);
  }

  const double base_mag = std::max(0.05, std::abs(pan));
  const int first_sign = (pan >= 0.0) ? 1 : -1;
  return aureo::clamp(static_cast<double>(first_sign) * base_mag, -1.0, 1.0);
}

inline void aureonoise_commit_phi_state(t_aureonoise* x, double pan, bool ok)
{
  if (!x || !x->p_phi_pan) return;
  const double mag = std::abs(pan);
  if (!ok && mag < 1.0e-6) {
    x->phi_last_mag = 0.0;
    x->phi_last_sign = 1;
    return;
  }
  x->phi_last_mag = mag;
  x->phi_last_sign = (pan >= 0.0) ? 1 : -1;
}

inline double aureonoise_compute_burst_weight(const t_aureonoise* x)
{
#if AUREO_BURST_HAWKES
  if (x && x->p_burst) {
    const double base = x->hawkes.base;
    return aureo::clamp01((x->hawkes.lambda - base) / (base * 4.0));
  }
#endif
  return 0.0;
}

inline double aureonoise_apply_burst_position(t_aureonoise* x,
                                              double pan,
                                              double dur_norm,
                                              double gap_norm,
                                              double& amp_scale)
{
  amp_scale = 1.0;
  const double burst_w = aureonoise_compute_burst_weight(x);
  if (burst_w <= 1.0e-6) return pan;

  const double clamped_dur = aureo::clamp01(dur_norm);
  const double clamped_gap = aureo::clamp01(gap_norm);

  const double sign = (pan >= 0.0) ? 1.0 : -1.0;
  double mag = std::abs(pan);

  const double center_pull = burst_w * (0.35 + 0.45 * clamped_dur);
  mag = mag * (1.0 - center_pull) + center_pull * 0.12;

  const double edge_push = burst_w * (0.35 + 0.30 * (1.0 - clamped_dur));
  mag = aureo::clamp(mag + edge_push * (1.0 - clamped_gap), 0.0, 1.0);

  amp_scale = 1.0 + burst_w * (0.6 + 0.4 * clamped_dur);
  return aureo::clamp(sign * mag, -1.0, 1.0);
}
