#pragma once

#include <algorithm>

struct t_aureonoise;

namespace aureonoise {
namespace control {

struct TempoState {
  double fallback_rate = 4.0;   // Hz
  double manual_rate = 4.0;     // Hz (p_rate or fallback)
  double sync_rate = 0.0;       // Hz when host sync active
  double current_rate = 4.0;    // Final rate used by scheduler
  double bpm_target = 120.0;
  double bpm_current = 120.0;
};

void tempo_reset(t_aureonoise* x);
void tempo_set_manual(t_aureonoise* x, double rate_hz);
void tempo_set_bpm(t_aureonoise* x, double bpm);
void tempo_update(t_aureonoise* x, double block_dt);

double tempo_manual_rate(const t_aureonoise* x);
double tempo_current_rate(const t_aureonoise* x);
bool tempo_sync_active(const t_aureonoise* x);

} // namespace control
} // namespace aureonoise

