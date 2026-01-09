#pragma once

#include <array>
#include <vector>

#include "aureo_core/aureo_noise.hpp"

namespace aureonoise {

struct HarmonyConfig {
  bool prime_comb = false;
  bool time_quant = false;
  bool time_strict = false;
  double reference_gap = 0.0;
  double reference_length = 0.0;
  int active_partials = aureo::NoiseAureoState::kMaxPartials;
};

struct HarmonyState {
  HarmonyConfig config{};
  std::array<bool, aureo::NoiseAureoState::kMaxPartials> partial_mask{};
  std::vector<double> time_ratios;
  int active_count = 0;

  bool primes_only() const { return config.prime_comb; }
  bool quantize_time() const { return config.time_quant; }
  bool strict_grid() const { return config.time_strict; }
};

HarmonyState make_harmony_state(const HarmonyConfig& cfg);
double quantize_time_samples(const HarmonyState& state, double samples, double reference);
inline double quantize_gap_samples(const HarmonyState& state, double samples)
{
  return quantize_time_samples(state, samples, state.config.reference_gap);
}
inline double quantize_length_samples(const HarmonyState& state, double samples)
{
  return quantize_time_samples(state, samples, state.config.reference_length);
}

} // namespace aureonoise

