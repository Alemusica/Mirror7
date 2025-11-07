# Mirror7 Build

Standalone JUCE-based build of the Mirror7 texture synth, extracted from the `aureonoise_tilde` mono-repo and reorganised into macro-modules for easier maintenance.

## Layout

- `source/engine` &ndash; DSP core (`Mirror7Engine`) plus shared aureonoise dependencies.
- `source/plugin` &ndash; JUCE audio processor entry point and parameter wiring.
- `source/gui` &ndash; Editor components and UI bindings.
- `source/config` &ndash; Parameter IDs and other header-only configuration.
- `resources/spatial_assets` &ndash; Spatial impulse/geometry data required by the engine.
- `scripts/post_build_fix_moduleinfo.py` &ndash; Cleans JUCE's generated `moduleinfo.json` files after builds.

## Prerequisites

- CMake ≥ 3.22
- JUCE with CMake exports available (install via package manager or add as submodule and run `cmake --install`)
- Python ≥ 3.8 (for the post-build manifest fixer)
- Xcode command-line tools (macOS) or equivalent toolchain on other platforms

The project reuses core DSP headers from the original mono-repo. Set `AUREONOISE_ROOT` to the root of `aureonoise_tilde` if it is not located next to this folder:

```bash
cmake -B build -S . -DAUREONOISE_ROOT=/path/to/aureonoise_tilde
```

## Building

Pick a JUCE installation that already exported its CMake files (see `docs/AGENT_README.md` for examples) and point CMake at it:

```bash
cmake -B build -S . \
  -DJUCE_DIR=/path/to/JUCE/install/lib/cmake/JUCE-8.0.10 \
  -DAUREONOISE_ROOT=/path/to/aureonoise_tilde   # optional if it already lives next to this folder
cmake --build build --config Release
```

The plugin artefacts are written to `build/Mirror7_artefacts`. After each build the helper script sanitises JUCE's `moduleinfo.json` files.

## Testing

Mirror7 now ships with a JUCE-based smoke test that instantiates `Mirror7Engine` both at 1× and 4× oversampling and asserts the render is finite:

```bash
cmake --build build --target mirror7_engine_smoke
ctest --output-on-failure --test-dir build
```

See `docs/TESTING.md` for details and guidance on adding more coverage.

## Documentation

- `docs/AGENT_README.md` – quick-start for automation agents and contributors (paths, commands, conventions).
- `docs/TESTING.md` – how to configure, run, and extend the test suite.

## Next Steps

- Integrate automated tests/CI via CTest.
- Package dependencies (e.g., beta7_tools) as submodules for full standalone builds.
- Add preset management utilities and cross-platform asset handling.
