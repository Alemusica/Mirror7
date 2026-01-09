# Mirror7 DSP Analysis

This note distils the mathematics hiding in `Mirror7Engine` so it is easier to reason about parameter iteration before refactoring.

## 1. Event Scheduling

### 1.1 Tempo / Host Sync (`Mirror7Engine::updateTempo`)

- Manual rate: `tempo.manual_rate ← clamp(rateHz, 0.01, kMaxEventRateHz)` whenever `rateHz > 0`.
- Sync mode kicks in when `syncEnable = true` **and** `rateHz ≈ 0`. In that branch:
  - Host BPM is slewed with `α = exp(-Δt / τ)` where `τ = max(0.01, syncSlew)` and `Δt` is the current block duration.
  - Event rate becomes `tempo.sync_rate = clamp( (bpm_current / 60) * syncDivision, 0.01, kMaxEventRateHz )`.
  - Transport start edges zero the phase accumulator (`samples_to_next = 1; gap_elapsed = 0`).

### 1.2 Gap Distribution (`scheduleGapSamples`)

Two regimes:
1. **Sync-enabled**: deterministic spacing `gap_samples = round(fs / tempo.current_rate)`.
2. **Free-running**: exponential inter-arrival with a slow φ-modulated rate
   \[
   \lambda(t) = \text{baseRate} \cdot \left[1 + 0.2 \sin\left(2\pi t \cdot \varphi^{-1}\right)\right]
   \]
   and `gap = min(cap, -ln(U)/λ)` with `cap = min(30s, 4/baseRate)` to avoid unbounded tails.

### 1.3 Grain Duration Mapping (`mapLenSamples`)

\[
L = \text{baseLenMs} \cdot 10^{-3} \cdot fs \cdot \varphi^{(2u-1)\cdot lenPhi}
\]
clamped to `[kMinGrainSamples, 4·fs]`. `lenPhi` controls the amount of golden-ratio stretching.

### 1.4 Hemispheric Coupling

The coupling parameter steers new pan targets towards the previous grain:
\[
\text{hemi} = \text{hemisCoupling} \cdot \left[0.35 + 0.65 \cdot \frac{r}{r+1}\right], \quad r = \frac{\text{gap}}{\text{prevDur}}
\]
and the raw pan sample is mixed as `(1-hemi) * pan - hemi * prev_pan`.

### 1.5 Low-Discrepancy Sources

`w_phi`, `w_s2`, `w_pl` step by `φ⁻¹`, `√2⁻¹`, and `plastic⁻¹` respectively, yielding quasi-random sequences with minimal correlation. Additional jitter draws come from the RNG (`u4…u6`) to randomise noise kinds, IPD, etc.

## 2. Dialogue & Parameter Iteration

The dialogue system keeps state (pan sign, magnitude, Fibonacci indices) and adjusts new grains via `dialogueEvaluate`:
- Alternation: forcing pan sign flips by comparing `dialogue.last_sign`.
- Fibonacci handshake: target ratios are chosen from `{1,1,2,3,5,8,13,21}`. The score blends into gap predictions through `dialogue.fib_gap_queue`.
- Resulting scaling: `amp_scale`, `dur_scale`, `pan`, `coherence`, etc. feed back via `dialogueCommit`.

Notable iterative parameters:
- `phi_last_mag/sign`: updated in `commitPhiState` to keep φ-pan alternation.
- `pinnaMix`: smoothed with `pinnaMix += (target - current) / (fs·0.02)`.
- `modalEnergy`: single-pole tracker `E ← 0.995·E + 0.005·e`.
- `panSmoothed`: 20 ms linear ramp to absorb sudden host automation changes.

## 3. Spatial Decoding

### 3.1 Binaural Modes

- **Legacy** (`phiMode = false` or invalid geometry): uses `aureo::compute_binaural`, parameterised by `width`, `itd_us`, `ild_db`, random phase perturbations.
- **Φ Mode**: `aureonoise::spatial::phi::compute_head_result` returns `itd_samples`, `lateral`, etc. Temporal hysteresis uses `hemi` to mix current head result with `prev_itd/prev_ild`.

### 3.2 IPD / Shadow

- IPD coefficient: `ipd_coeff = clamp(0.18 + 0.55·IPD + spread·(shape-0.5), 0, 0.95)` and scaled by `0.7 + 0.3·focus`. Implemented as a pair of allpasses.
- Shadow filter: first-order single-pole low-pass with `α = exp(-2π fc / fs)` and `fc = map_phi_range(800, fc_max, 1 - |pan|·shadowAmt)`.

### 3.3 CEFG (Early Reflections)

When `spatialProfile = 1` the code loads `uteroid_phi.cefg` and precomputes up to 16 taps per grain. Each tap stores:
`delay`, ITD offsets, per-ear gains, and a simple LP integrator with coefficient `α = exp(-2π fc / fs)` where `fc` decays as the listener distance grows.

## 4. Noise / Modulation Path

### 4.1 Noise Families

Per block, the Aureonoise engine is reconfigured:
- **Quantum**: number of active partials = `round( minP + detail·(maxP - minP) )`; decay mapped via `map_phi_range(0.12, 0.52, 1-detail)`.
- **Aureo**: harmonics clamped to `[1, kMaxPartials]`, base frequency tied to `fs`.
- **Velvet**: forced when `setVelvetAmount` is called from the processor UI.

### 4.2 VHS Modulation

`vhsWow` / `vhsFlutter` map to LFO increments:
\[
\Delta f_\text{wow} = \text{map\_phi\_range}(0.1, 1.5, vhsWow) / fs,\quad
\Delta f_\text{flt} = \text{map\_phi\_range}(7, 12, vhsFlutter) / fs
\]
and the modulation factor is `0.5·sin(2π·lfo_wow) + 0.5·sin(2π·lfo_flt)`.

### 4.3 Glitch Parameters

- Grain kind = `aureo::choose_kind(glitchMix, u4)`.
- Sample-rate hold length = `map_sr_hold_base(srCrush, u1)` (log-style mapping).
- Bit depth = `map_bits(bitCrush)` (maps `[0,1]` to `[4,16]` bits, then `levels = 2^{bits-1}-1`).

## 5. Parameter Iteration Summary

| Parameter | Iteration Mechanism | Effect |
|-----------|--------------------|--------|
| `rateHz`, `sync*` | Updates `tempo.current_rate`, influences `scheduleGapSamples`. | Controls mean inter-arrival. |
| `hemisCoupling` | Mixes current random pan with previous pan & ITD/ILD. | Encourages alternating hemispheres. |
| `phiPan`, `phiRatio`, `phi*` | Stored state `phi_last_mag/sign` ensures golden-ratio pan magnitudes; geometry built once per `prepare`. | Enforces φ-based spatial choreography. |
| `pinnaDepth`, `pinnaOn` | Smooth mix parameter `pinnaMix` + filter update flag. | Avoids zipper noise when toggling pinna notch. |
| `modal*` | `modalEnergy` integrator and dry/wet crossfade. | Keeps tail energy smooth between grains. |
| `dialogue*` | Maintains last pan sign/magnitude, Fibonacci queues, coherence counters. | Gives conversational pacing to grain bursts. |
| `vhsWow/flt` | Two phase accumulators updated sample-by-sample. | Drives VHS grain kind behaviour. |

These relationships identify where future refactors should extract subsystems (tempo, dialogue, spatial, noise) without changing behaviour.

The codebase now contains dedicated helpers (`DialogueSystem`, `NoiseController`, `Spatializer`) that map to the sections above, making it easier to evolve each subsystem independently and to test them in isolation.
