#pragma once

#include "aureonoise3beta/core/state.hpp"
#include "aureo_core/aureo_math.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace aureonoise::dialogue {

namespace {
constexpr std::array<double, 8> kFibonacciWalk{{1.0, 1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0}};
constexpr int kFibonacciMaxIndex = static_cast<int>(kFibonacciWalk.size()) - 1;
}

struct Result {
  double pan = 0.0;
  double amp_scale = 1.0;
  double dur_scale = 1.0;
  double next_coherence = 1.0;
  bool   handshake = false;
  int    sign = 1;
  double base_amp = 1.0;
  double base_dur = 0.0;
  double gap_samples = 0.0;
  double timestamp = 0.0;
  double handshake_score = 0.0;
  double handshake_ratio_target = 1.0;
  int    fib_index = 1;
  int    fib_direction = 1;
};

inline void reset(t_aureonoise* x)
{
  if (!x) return;
  x->dialogue = {};
  x->dialogue.coherence = 1.0;
  x->dialogue.fib_index = 1;
  x->dialogue.fib_direction = 1;
  x->dialogue.fib_ratio_target = aureo::kPhi;
  x->dialogue.fib_score = 0.0;
  x->dialogue.fib_gap_queue.fill(0.0);
  x->dialogue.fib_gap_queue_size = 0;
  x->dialogue.fib_gap_queue_pos = 0;
}

inline Result evaluate(t_aureonoise* x,
                       double timestamp_sec,
                       double pan,
                       double amp,
                       double dur_samples,
                       double gap_samples)
{
  Result res;
  res.pan = aureo::clamp(pan, -1.0, 1.0);
  res.amp_scale = 1.0;
  res.dur_scale = 1.0;
  res.next_coherence = 1.0;
  res.handshake = false;
  res.base_amp = amp;
  res.base_dur = dur_samples;
  res.gap_samples = gap_samples;
  res.timestamp = timestamp_sec;
  res.sign = (res.pan >= 0.0) ? 1 : -1;
  res.handshake_score = 0.0;
  res.handshake_ratio_target = aureo::kPhi;
  res.fib_index = 1;
  res.fib_direction = 1;

  if (!x || !x->p_dialogue_on) {
    return res;
  }

  auto& st = x->dialogue;
  const int fib_count = kFibonacciMaxIndex + 1;
  const auto clamp_index = [&](int idx) {
    return std::clamp(idx, 0, fib_count - 1);
  };

  int current_index = clamp_index(st.fib_index);
  int current_direction = (st.fib_direction >= 0) ? 1 : -1;
  if (current_direction == 0) current_direction = 1;
  int candidate_direction = current_direction;
  int raw_candidate_index = current_index + candidate_direction;
  if (raw_candidate_index < 0 || raw_candidate_index >= fib_count) {
    candidate_direction = -candidate_direction;
    raw_candidate_index = current_index + candidate_direction;
  }
  int candidate_index = clamp_index(raw_candidate_index);

  const double fib_prev = std::max(1.0, kFibonacciWalk[current_index]);
  const double fib_next = std::max(1.0, kFibonacciWalk[candidate_index]);
  double fib_ratio = fib_next / fib_prev;
  if (!std::isfinite(fib_ratio) || fib_ratio <= 0.0) {
    fib_ratio = aureo::kPhi;
  }

  res.handshake_ratio_target = fib_ratio;
  res.fib_index = current_index;
  res.fib_direction = current_direction;

  const double strength = aureo::clamp01(x->p_dialogue_strength);
  const double memory = aureo::clamp01(x->p_dialogue_memory);
  const double phi_mix = aureo::clamp01(x->p_dialogue_phi_mix);
  const double mag = std::max(1.0e-6, std::abs(res.pan));
  const bool has_prev = (std::abs(st.last_sign) > 0.5);

  const auto ratio_score = [](double ratio, double target) -> double {
    if (!std::isfinite(ratio) || ratio <= 0.0 || !std::isfinite(target) || target <= 0.0) return 0.0;
    const double log_diff = std::abs(std::log(ratio / target));
    constexpr double tolerance = 0.4;
    return std::clamp(1.0 - log_diff / tolerance, 0.0, 1.0);
  };

  double coherence_target = aureo::kInvPhi;

  if (has_prev && static_cast<int>(std::copysign(1.0, st.last_sign)) != res.sign) {
    const double prev_mag = std::clamp(st.last_mag, 0.05, 1.0);
    const double incoming_mag = mag;
    const double incoming_ratio = incoming_mag / std::max(prev_mag, 1.0e-6);
    const double amp_ratio = std::max(1.0e-6, amp) / std::max(1.0e-6, st.last_amp);
    const double dur_ratio = (st.last_dur > 0.0)
                               ? std::max(1.0e-6, dur_samples) / std::max(1.0e-6, st.last_dur)
                               : incoming_ratio;
    const double gap_ratio = (st.last_gap > 0.0)
                               ? std::max(1.0e-6, gap_samples) / std::max(1.0e-6, st.last_gap)
                               : incoming_ratio;

    const double pan_score = ratio_score(incoming_ratio, fib_ratio);
    const double amp_score = ratio_score(amp_ratio, fib_ratio);
    const double dur_score = ratio_score(dur_ratio, fib_ratio);
    const double gap_score = ratio_score(gap_ratio, fib_ratio);

    res.handshake_score = 0.55 * pan_score + 0.15 * amp_score + 0.15 * dur_score + 0.15 * gap_score;
    const double handshake_threshold = std::clamp(0.42 - 0.18 * strength, 0.22, 0.42);
    res.handshake = (res.handshake_score >= handshake_threshold);

    const double target_mag = aureo::clamp(
        (1.0 - strength) * incoming_mag +
            strength * (phi_mix * std::clamp(prev_mag * fib_ratio, 0.05, 1.0) + (1.0 - phi_mix) * incoming_mag),
        0.05,
        1.0);
    res.pan = static_cast<double>(res.sign) * target_mag;

    const double amp_target = std::max(1.0e-6, st.last_amp) * fib_ratio;
    const double amp_scale_target =
        aureo::clamp(amp_target / std::max(1.0e-6, res.base_amp), 0.15, 4.0);
    const double amp_interp = phi_mix * amp_scale_target + (1.0 - phi_mix);
    res.amp_scale = aureo::clamp((1.0 - strength) * 1.0 + strength * amp_interp, 0.15, 3.0);

    if (st.last_dur > 0.0) {
      const double dur_target = st.last_dur * fib_ratio;
      const double dur_scale_target =
          aureo::clamp(dur_target / std::max(1.0, res.base_dur), 0.25, 4.0);
      const double dur_interp = phi_mix * dur_scale_target + (1.0 - phi_mix);
      res.dur_scale = aureo::clamp((1.0 - strength) * 1.0 + strength * dur_interp, 0.25, 4.0);
    }

    coherence_target = std::clamp(fib_ratio, 0.6, 1.8);

    if (res.handshake) {
      res.fib_index = candidate_index;
      res.fib_direction = candidate_direction;
    }
  } else if (has_prev) {
    const double relaxed_target = std::clamp(st.last_mag * st.fib_ratio_target, 0.05, 1.0);
    const double relaxed = (1.0 - 0.5 * strength) * mag + 0.5 * strength * relaxed_target;
    res.pan = static_cast<double>(res.sign) * aureo::clamp(relaxed, 0.05, 1.0);
    coherence_target = 1.0;
  } else {
    coherence_target = 1.0;
  }

  if (st.pan_samples > 0) {
    const double bias = st.pan_sum / static_cast<double>(st.pan_samples);
    const double corr = strength * 0.4;
    res.pan = aureo::clamp(res.pan - corr * bias, -1.0, 1.0);
  }
  res.sign = (res.pan >= 0.0) ? 1 : -1;

  res.next_coherence = (1.0 - memory) * st.coherence + memory * coherence_target;
  return res;
}

inline void commit(t_aureonoise* x, const Result& res, bool accepted)
{
  if (!x || !x->p_dialogue_on) return;
  auto& st = x->dialogue;
  const double prev_gap = st.last_gap;
  if (!accepted) {
    ++st.failures;
    return;
  }

  st.last_sign = static_cast<double>(res.sign);
  st.last_mag = std::abs(res.pan);
  st.last_amp = std::max(1.0e-6, res.base_amp * res.amp_scale);
  st.last_dur = res.base_dur * res.dur_scale;
  st.last_gap = res.gap_samples;
  st.last_timestamp = res.timestamp;
  st.coherence = res.next_coherence;
  st.coherence_accum += st.coherence;
  ++st.coherence_samples;
  ++st.utterances;
  st.fib_index = std::clamp(res.fib_index, 0, kFibonacciMaxIndex);
  st.fib_direction = (res.fib_direction >= 0) ? 1 : -1;
  if (res.handshake) {
    ++st.handshakes;
    st.fib_ratio_target = res.handshake_ratio_target;
    st.fib_gap_queue.fill(0.0);
    st.fib_gap_queue_pos = 0;
    st.fib_gap_queue_size = 0;
    const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
    const double gap_cap = std::max(1.0, sr * 16.0);
    auto push_gap = [&](double gap) {
      if (st.fib_gap_queue_size >= st.fib_gap_queue.size()) return;
      const double clamped = aureo::clamp(gap, 1.0, gap_cap);
      st.fib_gap_queue[st.fib_gap_queue_size++] = clamped;
    };
    const double base_gap = (prev_gap > 1.0) ? prev_gap : std::max(1.0, res.gap_samples);
    push_gap(base_gap);
    const double phi_lo = std::max(1.0e-6, res.handshake_ratio_target / aureo::kPhi);
    const double phi_hi = std::max(phi_lo * 1.000001, res.handshake_ratio_target * aureo::kPhi);
    const double phi_blend = aureo::map_phi_range(phi_lo, phi_hi, aureo::clamp01(res.handshake_score));
    const double plastic_bias = std::pow(aureo::kPlastic, 0.05 * (res.handshake_score - 0.5));
    const double extend_ratio = aureo::clamp(phi_blend * plastic_bias, 0.25, 8.0);
    push_gap(base_gap * extend_ratio);
  } else {
    if (st.fib_gap_queue_pos >= st.fib_gap_queue_size) {
      st.fib_gap_queue_size = 0;
      st.fib_gap_queue_pos = 0;
      st.fib_gap_queue.fill(0.0);
    }
  }
  st.fib_score = res.handshake_score;
  st.pan_sum += res.pan;
  st.pan_abs_sum += std::abs(res.pan);
  ++st.pan_samples;
}

inline double handshake_ratio(const t_aureonoise* x)
{
  if (!x || x->dialogue.utterances <= 1) return 0.0;
  return static_cast<double>(x->dialogue.handshakes) / static_cast<double>(x->dialogue.utterances - 1);
}

inline double coherence_mean(const t_aureonoise* x)
{
  if (!x || x->dialogue.coherence_samples == 0) return 1.0;
  return x->dialogue.coherence_accum / static_cast<double>(x->dialogue.coherence_samples);
}

} // namespace aureonoise::dialogue
