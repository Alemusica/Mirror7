# Mirror7 — Φ-Based Binaural Texture Synth

> **Aureonoise Suite** — Unified repository for the Mirror7 JUCE plugin and aureonoise DSP ecosystem.

Mirror7 is a **granular texture synthesizer** with unique spatial characteristics based on the **Golden Ratio (φ)**. It creates "conversational" stereo textures where grains alternate between hemispheres following Fibonacci proportions.

## ✨ Key Features

- **Φ-Based Spatial Dialogue** — Grains alternate L/R with magnitudes following golden ratio proportions
- **Binaural ITD/ILD** — Physics-based interaural time/level differences with head geometry model
- **CEFG Early Reflections** — Convolution-free early reflections from spatial profiles
- **Multiple Noise Engines** — White, Pink, Brown, Aureo (φ-harmonic), Quantum, Velvet
- **DialogueSystem** — Fibonacci-based coherence tracking for organic grain pacing
- **Host Sync** — Lock grain rate to DAW tempo with division and slew
- **Oversampling** — 1×/2×/4× for clean high frequencies

## 📁 Project Structure

```
Mirror7_build/
├── source/                    # JUCE Plugin Source
│   ├── engine/               # DSP Core
│   │   ├── Mirror7Engine     # Main granular engine
│   │   ├── DialogueSystem    # Φ-based alternation logic
│   │   ├── Spatializer       # Binaural + CEFG
│   │   └── NoiseController   # Noise generators
│   ├── plugin/               # JUCE AudioProcessor
│   └── gui/                  # Editor UI
│
├── modules/                   # Advanced DSP Modules (from vellutoblu~)
│   ├── harmony/              # Prime comb, time quantization (TODO: integrate)
│   ├── control/              # Tempo state machine
│   ├── core/                 # Modal engine, state definitions
│   └── dsp/                  # Dialogue, scheduler, spatial utils
│
├── legacy/                    # Historical Reference Code
│   ├── aureo_core_v1/        # Original aureonoise headers (GitHub)
│   └── max_external/         # Original Max external source
│
├── resources/
│   └── spatial_assets/       # CEFG profiles, room meshes
│
├── tests/                     # CTest suite
│   ├── EngineSmokeTest       # Basic render sanity
│   ├── DialogueTest          # Alternation heuristics
│   └── SpatializerTest       # ITD/ILD symmetry
│
└── docs/
    ├── DSP_ANALYSIS.md       # Mathematical breakdown
    ├── ECOSYSTEM_INVENTORY.md # Full project inventory
    └── TESTING.md            # Test guide
```

## 🛠️ Building

### Prerequisites

- CMake ≥ 3.22
- JUCE 8.x with CMake exports
- C++20 compiler (Xcode CLT on macOS)
- Python ≥ 3.8 (for post-build scripts)

### Dependencies

The engine depends on `aureonoise_tilde` headers. Set `AUREONOISE_ROOT`:

```bash
# Clone aureonoise_tilde if not already present
git clone https://github.com/Alemusica/aureonoise_tilde.git ../aureonoise_tilde

# Configure
cmake -B build -S . \
  -DJUCE_DIR=/path/to/JUCE/lib/cmake/JUCE-8.x \
  -DAUREONOISE_ROOT=../aureonoise_tilde

# Build
cmake --build build --config Release
```

### Testing

```bash
cmake --build build --target mirror7_engine_smoke mirror7_dialogue_test mirror7_spatial_test
ctest --output-on-failure --test-dir build
```

## 🎛️ Parameters (~70)

| Group | Key Parameters |
|-------|----------------|
| **Timing** | `rate_hz`, `base_ms`, `len_phi`, `hemis_coupling` |
| **Spatial** | `itd_us`, `ild_db`, `spat_ipd`, `spat_shadow`, `phi_mode` |
| **Noise** | `noise_mode`, `aureo_mix`, `quantum_detail`, `velvet_amt` |
| **Dialogue** | `dialogue_on`, `dialogue_strength`, `dialogue_memory` |
| **Glitch** | `glitch_mix`, `sr_crush`, `bit_crush`, `vhs_wow` |
| **Modal** | `modal_on`, `modal_preset`, `modal_mix` |
| **Sync** | `sync_enable`, `sync_division`, `sync_slew` |

## 📋 Roadmap

### ✅ Completed
- [x] JUCE plugin with full parameter set
- [x] Modular engine (Dialogue, Spatializer, NoiseController)
- [x] Preset save/load (.mir7preset)
- [x] Test suite
- [x] Integration of legacy code

### 🔴 Phase 1 — Bug Fix
- [ ] Enable `burst` parameter (currently hardcoded to 0)
- [ ] Document magic numbers

### 🟠 Phase 2 — Feature Port
- [ ] Integrate `HarmonySystem` from `modules/harmony/`
- [ ] Add `prime_comb` (prime-only harmonics)
- [ ] Add `time_quant` / `time_strict` (Fibonacci grid)

### 🟡 Phase 3 — UX
- [ ] Collapsible parameter sections
- [ ] Grain activity visualizer
- [ ] Factory presets

### 🟢 Phase 4 — CI/CD
- [ ] GitHub Actions build
- [ ] Cross-platform testing

## 🔗 Related Projects

| Project | Description |
|---------|-------------|
| [aureonoise_tilde](https://github.com/Alemusica/aureonoise_tilde) | Parent monorepo (Max externals) — *archived* |
| [phiverb](https://github.com/Alemusica/phiverb) | Φ-based algorithmic reverb |
| aureo-factory | DSP core library (Python/C++) |

## 📜 License

MIT License — © 2025 Alemusica

---

*"Everything is driven by non-periodic relationships: φ, √2, plastic constant, primes."*
