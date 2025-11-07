# Agent README

This document is meant for automation agents (Codex, CI bots, etc.) that need to work on Mirror7 without a lengthy hand-off.

## Quick Facts

- **Project**: Mirror7 texture synth (JUCE + custom aureonoise DSP).
- **Entry Points**: `source/plugin` (AudioProcessor), `source/gui` (Editor), `source/engine` (DSP core).
- **External deps**:
  - `AUREONOISE_ROOT` → clone of `aureonoise_tilde`. Default points to `../aureonoise_tilde`.
  - `JUCE_DIR` → directory containing `JUCEConfig.cmake` (e.g. `.../JUCE/install/lib/cmake/JUCE-8.0.10`).
- **Build dir**: use `build/` (already ignored). Never commit artefacts under `Mirror7_artefacts/`.

## Standard Workflow

1. **Configure** (only once per session, sets cache variables):
   ```bash
   cmake -S . -B build \
     -DAUREONOISE_ROOT=/path/to/aureonoise_tilde \
     -DJUCE_DIR=/path/to/JUCE/install/lib/cmake/JUCE-8.0.10
   ```
2. **Build plug-in + engine**:
   ```bash
   cmake --build build --target Mirror7
   ```
3. **Run tests** (see `docs/TESTING.md` for more detail):
   ```bash
   cmake --build build --target mirror7_engine_smoke
   ctest --output-on-failure --test-dir build
   ```

## Common Tasks

| Task | Where to edit | Notes |
| --- | --- | --- |
| Add / tweak parameters | `source/config/ParamIDs.h`, `PluginProcessor.cpp` (layout + mapping), `PluginEditor.cpp` | Keep layout + editor in sync; apply defaults through `Mirror7Engine::Params`. |
| DSP changes | `source/engine/Mirror7Engine.*` | Split out helpers whenever possible; add focused tests for new scheduling/spatial logic. |
| Build / toolchain updates | `CMakeLists.txt`, `cmake/` (future) | Prefer cache variables instead of hard-coding paths; gate expensive options behind `BUILD_TESTING` etc. |
| Asset handling | `resources/` | Spatial profiles live under `resources/spatial_assets`; keep binary data out of git. |

## Style & Conventions

- C++20, JUCE naming (camel case for members, `snake_case` for locals if already used in the file).
- Comment only when intent is non-obvious (e.g., DSP tricks, phi geometry math).
- Prefer deterministic paths and `CACHE` variables for everything that lives outside the repo.
- Before sending patches, ensure:
  - `cmake --build build --target mirror7_engine_smoke`
  - `ctest --output-on-failure --test-dir build`

## Gotchas

- The aureonoise mono-repo recently renamed `aureonoise3beta` to `vellutoblu~`. `CMakeLists.txt` auto-detects this, but you can override `AUREONOISE_LEGACY_DIR` manually if needed.
- JUCE is not vendored; point `JUCE_DIR` at any existing installation (the repo reuses installations from other projects).
- Generated plug-in bundles land in `build/Mirror7_artefacts`. They are ignored; do not clean them with destructive commands—some developers inspect them between builds.
