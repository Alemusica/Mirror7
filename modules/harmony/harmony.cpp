#include "harmony.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace aureonoise {

namespace {

constexpr double kLooseTolerance = 0.12;
constexpr int kLooseMaxMultiplier = 8;
constexpr int kStrictMaxMultiplier = 4;

std::vector<double> loose_time_ratios()
{
  return {
    1.0,
    aureo::kInvPhi,
    1.0 / (aureo::kPhi * aureo::kPhi),
    aureo::kPhi,
    aureo::kPhi * aureo::kPhi,
    2.0 / 3.0,
    3.0 / 5.0,
    5.0 / 8.0,
    8.0 / 13.0,
    3.0 / 2.0,
    5.0 / 3.0,
    8.0 / 5.0,
    13.0 / 8.0
  };
}

std::vector<double> strict_time_ratios()
{
  return {
    aureo::kInvPhi,
    1.0 / (aureo::kPhi * aureo::kPhi),
    1.0,
    aureo::kPhi,
    aureo::kPhi * aureo::kPhi
  };
}

std::array<bool, aureo::NoiseAureoState::kMaxPartials>
make_partial_mask(const HarmonyConfig& cfg, int& active_count)
{
  std::array<bool, aureo::NoiseAureoState::kMaxPartials> mask{};
  active_count = 0;
  const int limit = std::clamp(cfg.active_partials, 1, aureo::NoiseAureoState::kMaxPartials);
  for (int i = 0; i < aureo::NoiseAureoState::kMaxPartials; ++i) {
    const bool within = (i < limit);
    const bool prime_ok = (!cfg.prime_comb) || (i == 0) || aureo::NoiseAureoState::is_prime_index(i + 1);
    const bool on = within && prime_ok;
    mask[static_cast<size_t>(i)] = on;
    if (on) ++active_count;
  }
  if (active_count == 0) {
    mask[0] = true;
    active_count = 1;
  }
  return mask;
}

} // namespace

HarmonyState make_harmony_state(const HarmonyConfig& cfg)
{
  HarmonyState state;
  state.config = cfg;
  state.partial_mask = make_partial_mask(cfg, state.active_count);
  if (cfg.time_quant) {
    state.time_ratios = cfg.time_strict ? strict_time_ratios() : loose_time_ratios();
  }
  return state;
}

double quantize_time_samples(const HarmonyState& state, double samples, double reference)
{
  const double fallback_base = (state.config.reference_length > 1.0) ? state.config.reference_length
                                : state.config.reference_gap;
  double base = (reference > 1.0) ? reference : fallback_base;

  if (!state.quantize_time() || samples <= 0.0 || base <= 1.0 || state.time_ratios.empty()) {
    return samples;
  }

  const double normalized = samples / base;
  if (normalized <= 1.0e-6) return samples;

  const int maxMult = state.strict_grid() ? kStrictMaxMultiplier : kLooseMaxMultiplier;

  double best_norm = normalized;
  double best_err = std::numeric_limits<double>::max();

  for (double ratio : state.time_ratios) {
    for (int mult = 1; mult <= maxMult; ++mult) {
      const double candidate = ratio * static_cast<double>(mult);
      const double diff = std::abs(candidate - normalized);
      if (diff < best_err) {
        best_err = diff;
        best_norm = candidate;
      }
      const double inv_candidate = ratio / static_cast<double>(mult);
      if (inv_candidate > 1.0e-6) {
        const double diff_inv = std::abs(inv_candidate - normalized);
        if (diff_inv < best_err) {
          best_err = diff_inv;
          best_norm = inv_candidate;
        }
      }
    }
  }

  if (!state.strict_grid()) {
    const double rel_err = best_err / normalized;
    if (rel_err > kLooseTolerance) return samples;
  }

  double quantized = best_norm * base;
  const double max_allowed = base * 16.0;
  if (quantized < 1.0) quantized = 1.0;
  if (quantized > max_allowed) quantized = max_allowed;
  return quantized;
}

} // namespace aureonoise
