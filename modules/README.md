# Advanced DSP Modules

These modules are extracted from `vellutoblu~` (aureonoise beta5) and contain advanced features not yet fully integrated into the Mirror7 JUCE plugin.

## Module Overview

### `harmony/`
**Status:** 🟡 Not integrated

Harmonic system with:
- **Prime Comb** (`prime_comb`) — Limits partials to prime indices only
- **Time Quantization** (`time_quant`) — Snaps gap/duration to Fibonacci ratios
- **Strict Grid** (`time_strict`) — Forces exact φ, φ², 1/φ ratios only

```cpp
// Example usage
HarmonyConfig cfg;
cfg.prime_comb = true;
cfg.time_quant = true;
auto state = make_harmony_state(cfg);
double quantized = quantize_gap_samples(state, raw_gap);
```

### `control/`
**Status:** ✅ Partially integrated (tempo logic in Mirror7Engine)

Tempo state machine for host sync:
- Manual rate vs sync rate
- BPM slewing
- Transport edge detection

### `core/`
**Status:** ✅ Integrated

- `modal_engine.hpp` — Modal resonator (Wood/Metal/Glass presets)
- `state.hpp` — Max external state definitions (reference)

### `dsp/`
**Status:** 🟡 Partially integrated

- `dialogue.hpp` — Reference dialogue implementation (Mirror7 has its own)
- `scheduler.hpp` — Grain scheduling with Hawkes process burst
- `spatial_utils.hpp` — Binaural utilities
- `dsp.cpp` — Full Max DSP implementation (reference)

## Integration TODO

1. **Harmony → Mirror7Engine**
   - Add `prime_comb`, `time_quant`, `time_strict` parameters
   - Wire into grain scheduling

2. **Burst (Hawkes) → Mirror7Engine**
   - Currently `burst_w = 0.0` is hardcoded
   - Port scheduler.hpp burst logic

3. **Grain Reporting → PluginEditor**
   - Add visual feedback for grain activity
   - Port report_log system from legacy
