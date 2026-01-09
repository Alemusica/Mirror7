# Mirror7 Knowledge Mapping

## 📊 Database Stats (1670 items)

| Source | Count | Relevance |
|--------|-------|-----------|
| arxiv | 1164 | ⭐⭐⭐ Core research |
| pubmed | 189 | ⭐⭐⭐ Neuroscience |
| dsp.stackexchange | 55 | ⭐⭐⭐ Implementation |
| physics.stackexchange | 61 | ⭐⭐ Acoustics |
| sound.stackexchange | 10 | ⭐⭐ Audio production |
| stackoverflow | 82 | ⭐ Code snippets |

---

## ⭐ TOP 10 PAPERS FOR MIRROR7

### 1. Hemispheric Two-Channel Code (MUST READ)
**ID:** `knowledge:arxiv_2111_04637v4`  
**Title:** A hemispheric two-channel code accounts for binaural unmasking in humans

> *"The neurocomputational process that underlies binaural unmasking is based on **two hemispheric channels** that encode interaural differences in their relative neuronal activity."*
>
> Key insight: **complex-valued correlation coefficient** - mathematical model accounting for 98% of variance!

**Application to Mirror7:**
- `DialogueSystem` models this two-channel code
- The "handshake" mechanism mimics interhemispheric correlation
- **TODO:** Implement complex-valued correlation metric

---

### 2. Binaural Beats Interhemispheric Coherence (MUST READ)
**ID:** `knowledge:pubmed_26541421`  
**Title:** Binaural beats increase interhemispheric alpha-band coherence between auditory cortices

> *"BBs enhanced interhemispheric coherence between the auditory cortices. Beat frequencies in the **alpha (10 Hz) and theta (4 Hz)** range both increased interhemispheric coherence selectively at alpha frequencies."*
>
> Key insight: **Not entrainment but binaural integration!**

**Application to Mirror7:**
- Target beat frequencies: **4 Hz (theta), 10 Hz (alpha)**
- Grain scheduling should create implicit beat in this range
- Effect is "binaural integration" not brainwave entrainment

---

### 3. Corpus Callosum Transfer Asymmetry
**ID:** `knowledge:pubmed_29687282`  
**Title:** Early asymmetric inter-hemispheric transfer in the auditory network

> *"Asymmetrical callosal connectivity favoring the **right-to-left hemisphere direction**. This asymmetry might contribute to enhancement of left lateralization for language processing."*

**Application to Mirror7:**
- R→L transfer is faster than L→R
- DialogueSystem should favor R→L direction for "handshake"
- ~60ms latency difference between hemispheres

---

### 4. HRTF Individualization Survey
**ID:** `knowledge:arxiv_HRTF_Individualization`  
**Title:** HRTF Individualization: A Survey

> Comprehensive review of HRTF personalization techniques for binaural synthesis.

**Application to Mirror7:**
- `Spatializer` Φ-mode geometry approximates HRTF
- Consider adding HRTF preset selection
- Pinna filtering in CEFG mode

---

### 5. ITD-ILD Transformation
**ID:** Multiple papers  
**Key finding:** Low-frequency ITD can be transformed to ILD for processing

**Application to Mirror7:**
- Current `Spatializer` uses both ITD and ILD
- Paper suggests they can be unified mathematically
- Consider `Evaluation of an ITD-to-ILD Transformation` approach

---

### 6-10. Additional Key Papers

| # | Title | Key Insight |
|---|-------|-------------|
| 6 | Natural statistics of binaural sounds | Real-world ITD/ILD distributions |
| 7 | Dichotic listening hemispheric lateralization | Ear advantage patterns |
| 8 | EEG coherence interhemispheric sleep/wake | Coherence varies by state |
| 9 | Binaural beats anxiety treatment | Clinical validation at 10 Hz |
| 10 | Sound localization perturbed binaural | Adaptation mechanisms |

---

## 💻 IMPLEMENTATION GUIDES (StackExchange)

### Fractional Delay Line (Critical for ITD!)

**Q:** How to prevent "zipping" effect on modulated fractional delay?

**Problem:** Clicks/artifacts when delay length increases (pitch decreasing)

**Solution:**
```cpp
// Use all-pass interpolation instead of linear
// Or: Hermite interpolation for smoother transitions

float hermite_interp(float *buf, float pos, int size) {
    int i = (int)pos;
    float frac = pos - i;
    float y0 = buf[(i-1+size) % size];
    float y1 = buf[i];
    float y2 = buf[(i+1) % size];
    float y3 = buf[(i+2) % size];
    
    float c0 = y1;
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f*y1 + 2*y2 - 0.5f*y3;
    float c3 = 0.5f*(y3-y0) + 1.5f*(y1-y2);
    
    return ((c3*frac + c2)*frac + c1)*frac + c0;
}
```

**Apply to Mirror7:** `Spatializer::processItd()` delay line

---

### Phase vs Group Delay

**Key insight:** For ITD simulation, **group delay** matters more than phase delay for perception of transients.

---

### RBJ Audio EQ Cookbook

Referenced for biquad filter design - useful for:
- Pink noise filtering
- Pinna resonance modeling
- Tone shaping in `NoiseController`

---

## 🗺️ ROADMAP MAPPING

### Phase 1: Bug Fix & Cleanup

| Task | Knowledge Source | Action |
|------|------------------|--------|
| **Burst parameter hardcoded** | DSP analysis | Enable `burst_w` control |
| **Magic numbers documentation** | Papers #1, #2 | Document 0.42, 0.55, 0.618 origins |

**Relevant findings:**
- 0.618 → φ (golden ratio) - aesthetic/natural timing
- 10 Hz → alpha band coherence target
- 60 ms → interhemispheric transfer latency

---

### Phase 2: Feature Port

| Task | Knowledge Source | Action |
|------|------------------|--------|
| **HarmonySystem integration** | modules/harmony/ | Port from vellutoblu~ |
| **prime_comb** | Fibonacci papers | Prime-only harmonics |
| **time_quant/time_strict** | Phi papers | Fibonacci grid quantization |

**Relevant findings:**
- Fibonacci sequence in music perception (arxiv)
- Golden ratio aesthetic perception (pubmed)
- Irrational relationships avoid beating

---

### Phase 3: GUI & UX

| Task | Knowledge Source | Action |
|------|------------------|--------|
| **Coherence visualizer** | EEG coherence papers | Show L/R correlation |
| **Handshake indicator** | Paper #1 | Two-channel visualization |
| **Factory presets** | Clinical papers | "Relaxation", "Focus", "Meditation" |

**Relevant findings:**
- Alpha (10 Hz) → relaxation
- Theta (4 Hz) → meditation  
- Gamma (40 Hz) → focus/cognition

---

### Phase 4: Quality & CI

| Task | Knowledge Source | Action |
|------|------------------|--------|
| **Test coverage** | All papers | Validate against published data |
| **API documentation** | Implementation guides | Clear parameter explanations |

---

## 🧪 EXPERIMENTAL IDEAS FROM RESEARCH

### 1. Complex-Valued Correlation Metric
From paper #1: Implement `ComplexCorrelation` class that computes:
```
C = <s_L(t) · s_R*(t)> / sqrt(<|s_L|²> · <|s_R|²>)
```
Where `*` is complex conjugate. This could replace simple `handshake_score`.

### 2. Target Frequency Bands
From paper #2 and clinical studies:

| Band | Frequency | Target Effect |
|------|-----------|---------------|
| Delta | 0.5-4 Hz | Deep sleep |
| Theta | 4-8 Hz | Meditation, creativity |
| Alpha | 8-12 Hz | Relaxation, calm focus |
| Beta | 12-30 Hz | Alertness |
| Gamma | 30-100 Hz | Cognitive processing |

Mirror7 grain rate could target these bands!

### 3. R→L Hemisphere Priority
From paper #3: Right-to-left transfer is ~60ms faster.
`DialogueSystem` could weight R→L transitions higher in handshake scoring.

### 4. HRV Coherence (0.1 Hz)
From HeartMath research: Heart rate variability peaks at 0.1 Hz in coherent states.
Consider very slow LFO at 0.1 Hz for "body coherence" mode.

---

## 📚 PAPERS TO READ IN FULL

1. **arxiv:2111.04637v4** - Hemispheric two-channel code (PRIORITY 1)
2. **pubmed:26541421** - Binaural beats coherence (PRIORITY 1)
3. **pubmed:29687282** - Corpus callosum transfer (PRIORITY 2)
4. **arxiv:1402.4648v2** - Natural statistics of binaural sounds
5. **pubmed:28325532** - Binaural beats dental anxiety (clinical validation)

---

## 🔗 LINKS

- ArXiv paper: `https://arxiv.org/abs/2111.04637`
- PubMed: `https://pubmed.ncbi.nlm.nih.gov/26541421/`
- DSP SE: `https://dsp.stackexchange.com/`
