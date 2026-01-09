# Testing Mirror7

The testing story is intentionally minimal today—a single smoke test that instantiates `Mirror7Engine` with and without oversampling and asserts that it produces finite samples. This document explains how to run it and how to extend the suite.

## Prerequisites

- CMake ≥ 3.22
- JUCE with CMake exports. Set `-DJUCE_DIR=/path/to/JUCE/install/lib/cmake/JUCE-<version>`.
- `AUREONOISE_ROOT` pointing at a checkout of `aureonoise_tilde`. If the legacy layout changed, set `-DAUREONOISE_LEGACY_DIR=...` explicitly.

## Running the Smoke Test

1. Configure (once per environment):
   ```bash
   cmake -S . -B build \
     -DJUCE_DIR=/path/to/JUCE/install/lib/cmake/JUCE-8.0.10
   ```
   Use `-DAUREONOISE_ROOT=...` if your mono-repo lives somewhere else.

2. Build and run:
   ```bash
   cmake --build build --target mirror7_engine_smoke
   ctest --output-on-failure --test-dir build
   ```

The `mirror7_engine_smoke` test lives in `tests/EngineSmokeTest.cpp`. It exercises `Mirror7Engine::prepare`, `process`, and the oversampling wrapper; failures are printed to `stderr` for easy diagnosis.

## Additional Targeted Tests

- `mirror7_dialogue_test` (`tests/DialogueTest.cpp`) instantiates the standalone `DialogueSystem`, drives a golden-ratio alternation, and asserts that the coherence score stays healthy and the state flips hemispheres as expected.
- `mirror7_spatial_test` (`tests/SpatializerTest.cpp`) builds a Φ-geometry via `Spatializer`, probes symmetric azimuths, and verifies ITD/ILD are mirrored.

Build and run them the same way:

```bash
cmake --build build --target mirror7_dialogue_test mirror7_spatial_test
ctest --output-on-failure --test-dir build
```

## Extending the Suite

- **Add new tests** next to the smoke test (use straightforward `main()` programs or plug Catch2/Doctest if the dependency footprint stays light).
- **Register** each test with CMake:
  ```cmake
  add_executable(mirror7_dialogue_test tests/DialogueTest.cpp)
  target_link_libraries(mirror7_dialogue_test PRIVATE mirror7_engine)
  add_test(NAME mirror7_dialogue_test COMMAND mirror7_dialogue_test)
  ```
- **Keep runtime low** (< 3 seconds). Aim for deterministic offline renders.
- **Document** each addition in this file so future contributors can discover and run them quickly.
