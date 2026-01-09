# Legacy Code Archive

This folder contains historical reference code from the original `aureonoise_tilde` GitHub repository.

## Contents

### `aureo_core_v1/`
Original DSP headers from the first public release. These are simpler versions of the current `aureo_core` headers used by Mirror7.

**Notable differences from current:**
- `aureo_noise.hpp`: 86 lines vs 334 lines (current has more noise modes)
- Other headers are mostly identical

### `max_external/`
Original Max/MSP external source code (`aureonoise~.mxo`).

**Files:**
- `aureonoise_tilde.cpp` — Main DSP and Max object implementation
- `aureonoise_state.hpp` — State structures and grain definitions
- `aureonoise_attributes.cpp` — Max attribute definitions

**Useful features to potentially port:**
- Grain reporting system (`GrainReport`, `report_log`)
- Poisson enforcement with relaxation (`aureonoise_poisson_enforce`)

## Usage

This code is kept for reference and potential feature extraction. It is **not** compiled as part of the main build.

To integrate features from legacy code:
1. Identify the feature in legacy source
2. Port to appropriate module in `source/engine/` or `modules/`
3. Add tests in `tests/`
4. Update documentation
