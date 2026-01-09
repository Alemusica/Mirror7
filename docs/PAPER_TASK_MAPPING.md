# Mirror7 - Paper → Task Mapping

## 📊 Knowledge Base Status

| Metric | Value |
|--------|-------|
| **Total Items** | 1858 |
| **With Embeddings** | 721 (39%) |
| **GraphRAG** | ✅ Active |
| **Sources** | arxiv (1255), pubmed (286), dsp.SE (55), sound.SE (10) |

---

## ⭐ PAPER FONDAMENTALI PER UNIONE EMISFERI

### 1. Corpus Callosum & Transcallosal Transfer

| Paper | Key Finding | Applicazione Mirror7 |
|-------|-------------|---------------------|
| **The Callosal Relay Model** | R→L directed flow via CC in gamma-band (30-100 Hz) | DialogueSystem: favorire R→L |
| **Transcallosal connectivity PT** | Musicians show increased transcallosal white matter | Target: musicisti/produttori |
| **Interhemispheric auditory connectivity** | Theta-band coherence between temporal areas | DialogueSystem: theta 4-8 Hz |
| **Agenesis of CC studies** | CC cruciale per integrazione binaurale | Valida approccio Mirror7 |

### 2. Binaural Integration & Coherence

| Paper | Key Finding | Applicazione Mirror7 |
|-------|-------------|---------------------|
| **Binaural beats alpha coherence** (pubmed:26541421) | 10Hz/4Hz BBs aumentano coerenza alpha | Grain rate → 4-10 Hz target |
| **Hemispheric two-channel code** (arxiv:2111.04637) | Complex-valued correlation, 98% variance | Nuova metrica handshake |
| **EEG coherence sleep/wake** | Coherence varia con stato mentale | Presets per diversi stati |

### 3. 40 Hz Gamma Stimulation (MIT)

| Paper | Key Finding | Applicazione Mirror7 |
|-------|-------------|---------------------|
| **Gamma Entrainment Binds Higher-Order Regions** | GENUS riduce amiloide, protegge neuroni | Preset "Cognitive" 40 Hz |
| **Multisensory gamma promotes glymphatic clearance** | 40 Hz attiva pulizia cerebrale | Feature terapeutica |
| **40Hz visual/auditory stimulation** | Gamma entrainment nel hippocampus | Validazione scientifica |
| **40Hz tACS speech perception** | 40 Hz modula percezione speech | Audio processing benefit |

### 4. HeartMath & Heart-Brain Coherence

| Paper | Key Finding | Applicazione Mirror7 |
|-------|-------------|---------------------|
| **HeartMath HRV Biofeedback** | 0.1 Hz respiration → coherence | LFO a 0.1 Hz per body sync |
| **From Dysregulation to Coherence** | Autonomic signaling heart-brain | Slow modulation feature |
| **Global Coherence Initiative** | Collective coherence effects | Community feature? |

---

## 🗺️ MAPPING → ROADMAP TASKS

### Phase 1: Bug Fix & Cleanup

| Task | Paper Evidence | Priority |
|------|----------------|----------|
| **Fix burst_w hardcoded** | - | 🔴 HIGH |
| **Document magic numbers** | Two-channel code (0.55, 0.42) | 🔴 HIGH |
| | Callosal relay (gamma 30-100 Hz) | |
| | Alpha coherence (10 Hz) | |

### Phase 2: Feature Port

| Task | Paper Evidence | Priority |
|------|----------------|----------|
| **Integrate HarmonySystem** | Fibonacci papers, phi research | 🟠 MED |
| **Add time_quant/time_strict** | Non-periodic relationships | 🟠 MED |
| **Add 40 Hz gamma preset** | MIT GENUS research | 🔴 HIGH |
| **Add 10 Hz alpha preset** | Binaural coherence papers | 🔴 HIGH |

### Phase 3: GUI & UX

| Task | Paper Evidence | Priority |
|------|----------------|----------|
| **Coherence visualizer** | EEG coherence papers | 🟠 MED |
| **Handshake indicator** | Two-channel code | 🟠 MED |
| **Factory presets** | Clinical trial papers | 🟡 LOW |
| - Deep Meditation (4-8 Hz) | Theta coherence | |
| - Calm Focus (10 Hz) | Alpha coherence | |
| - Peak Performance (40 Hz) | Gamma GENUS | |
| - Heart Sync (0.1 Hz LFO) | HeartMath | |

### Phase 4: Quality & CI

| Task | Paper Evidence | Priority |
|------|----------------|----------|
| **EEG validation study** | All coherence papers | 🟢 FUTURE |
| **A/B testing φ timing** | Fibonacci aesthetics | 🟢 FUTURE |

---

## 🧪 EXPERIMENTAL FEATURES FROM RESEARCH

### 1. Complex-Valued Correlation (from Two-Channel Code)

```cpp
// Implementazione proposta
struct ComplexCorrelation {
    // C = <s_L(t) · s_R*(t)> / sqrt(<|s_L|²> · <|s_R|²>)
    std::complex<float> compute(const float* left, const float* right, int N);
    float getCoherenceScore() const;  // |C|
    float getPhaseAngle() const;      // arg(C)
};
```

### 2. Frequency Band Presets

```cpp
enum class BrainState {
    DeepMeditation,  // theta: 4-8 Hz grain rate
    CalmFocus,       // alpha: 8-12 Hz grain rate  
    HighPerformance, // gamma: 40 Hz grain rate
    HeartSync        // LFO: 0.1 Hz amplitude mod
};

void setPreset(BrainState state) {
    switch(state) {
        case BrainState::DeepMeditation:
            setGrainRate(6.0f);  // theta center
            break;
        case BrainState::CalmFocus:
            setGrainRate(10.0f); // alpha peak
            break;
        case BrainState::HighPerformance:
            setGrainRate(40.0f); // gamma GENUS
            break;
        case BrainState::HeartSync:
            setLfoRate(0.1f);    // HRV coherence
            break;
    }
}
```

### 3. R→L Hemisphere Priority (from Callosal Relay)

```cpp
// In DialogueSystem
float computeHandshake(float leftPan, float rightPan) {
    // R→L transfer ~60ms faster
    float rl_weight = 1.2f;  // Favor R→L transitions
    float lr_weight = 0.8f;
    
    if (previousPan > 0 && currentPan < 0) {
        // R→L transition
        return baseScore * rl_weight;
    } else if (previousPan < 0 && currentPan > 0) {
        // L→R transition
        return baseScore * lr_weight;
    }
    return baseScore;
}
```

---

## 📚 PAPERS DA LEGGERE (Priority Order)

1. **The Callosal Relay Model** - Meccanismo R→L via CC
2. **Binaural beats alpha coherence** - 10 Hz/4 Hz target
3. **Gamma Entrainment Binds Regions** - 40 Hz neuroprotection
4. **Hemispheric two-channel code** - Mathematical model
5. **HeartMath HRV Biofeedback** - 0.1 Hz resonance

---

## 🔗 Database Queries Utili

```sql
-- Tutti i paper su corpus callosum
SELECT title, content FROM knowledge 
WHERE content CONTAINS 'corpus callosum'
ORDER BY source;

-- Paper con embeddings per ricerca semantica
SELECT * FROM knowledge 
WHERE embedding != NONE 
AND array::len(embedding) > 0;

-- Usa GraphRAG function
fn::get_entity_context('corpus callosum');
fn::search_knowledge($embedding, 10);
```

---

## ✅ ACTION ITEMS

### Immediati (questa sessione)
- [ ] Leggere "Callosal Relay Model" completo
- [ ] Implementare preset 40 Hz gamma
- [ ] Aggiungere R→L weighting al DialogueSystem

### Prossima sessione
- [ ] Implementare ComplexCorrelation metric
- [ ] Aggiungere LFO 0.1 Hz per HeartSync
- [ ] Creare GUI coherence visualizer

### Futuro
- [ ] Studio EEG con Mirror7
- [ ] Clinical validation con HeartMath
