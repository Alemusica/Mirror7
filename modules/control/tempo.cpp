#include "aureonoise3beta/control/tempo.hpp"

#include "aureonoise3beta/core/state.hpp"
#include "aureo_core/aureo_math.hpp"

#include <cmath>
namespace aureonoise {
namespace control {

namespace {

constexpr double kMinRateHz = 0.01;
constexpr double kMaxRateHz = aureo::kMaxEventRateHz;
constexpr double kMinDivision = 0.0625;
constexpr double kMaxDivision = 32.0;
constexpr double kMinBpm = 10.0;
constexpr double kMaxBpm = 400.0;

inline double clamp_rate(double hz)
{
  return std::clamp(hz, kMinRateHz, kMaxRateHz);
}

inline double clamp_bpm(double bpm)
{
  return std::clamp(bpm, kMinBpm, kMaxBpm);
}

} // namespace

void tempo_reset(t_aureonoise* x)
{
  if (!x) return;
  auto& st = x->tempo;
  const double fallback = (x->p_rate > 1.0e-6) ? clamp_rate(x->p_rate) : clamp_rate((st.fallback_rate > 1.0e-6) ? st.fallback_rate : 4.0);
  st.fallback_rate = fallback;
  st.manual_rate = fallback;
  st.sync_rate = 0.0;
  st.current_rate = fallback;
  st.bpm_target = clamp_bpm(st.bpm_target);
  st.bpm_current = st.bpm_target;
}

void tempo_set_manual(t_aureonoise* x, double rate_hz)
{
  if (!x) return;
  auto& st = x->tempo;
  const double clamped = clamp_rate(rate_hz);
  if (clamped > 1.0e-6) {
    st.fallback_rate = clamped;
    st.manual_rate = clamped;
    st.current_rate = clamped;
  }
}

void tempo_set_bpm(t_aureonoise* x, double bpm)
{
  if (!x) return;
  auto& st = x->tempo;
  st.bpm_target = clamp_bpm(bpm);
  if (st.bpm_current <= 0.0) st.bpm_current = st.bpm_target;
}

void tempo_update(t_aureonoise* x, double block_dt)
{
  if (!x) return;
  auto& st = x->tempo;

  const bool sync_enabled = (x->p_sync_enable != 0);
  const bool use_sync = sync_enabled && (x->p_rate <= 1.0e-6);

  // Manual branch
  const double manual = (x->p_rate > 1.0e-6) ? clamp_rate(x->p_rate) : st.fallback_rate;
  st.manual_rate = clamp_rate(manual);

  if (!use_sync) {
    st.sync_rate = 0.0;
    st.current_rate = st.manual_rate;
    if (!sync_enabled)
      st.bpm_current = st.bpm_target = clamp_bpm(st.bpm_target);
    return;
  }

  st.bpm_target = clamp_bpm(st.bpm_target);
  const double slew = std::clamp(x->p_sync_slew, 0.0, 5.0);
  if (block_dt > 0.0 && slew > 0.0) {
    const double tau = std::max(0.01, slew);
    const double alpha = std::exp(-block_dt / tau);
    st.bpm_current = st.bpm_target + (st.bpm_current - st.bpm_target) * alpha;
  } else {
    st.bpm_current = st.bpm_target;
  }

  const double division = std::clamp(x->p_sync_division, kMinDivision, kMaxDivision);
  st.sync_rate = clamp_rate((st.bpm_current / 60.0) * division);
  st.current_rate = st.sync_rate;
}

double tempo_manual_rate(const t_aureonoise* x)
{
  if (!x) return kMinRateHz;
  return x->tempo.manual_rate;
}

double tempo_current_rate(const t_aureonoise* x)
{
  if (!x) return kMinRateHz;
  return x->tempo.current_rate;
}

bool tempo_sync_active(const t_aureonoise* x)
{
  if (!x) return false;
  return (x->p_sync_enable != 0) && (x->p_rate <= 1.0e-6) && (x->tempo.sync_rate > 1.0e-6);
}

} // namespace control
} // namespace aureonoise

