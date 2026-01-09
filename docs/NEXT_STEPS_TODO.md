# Mirror7 - Next Steps & Implementation Plan

## 🚨 ISSUE CRITICO: Embeddings Mancanti

### Problema
- **1858 items** nel database
- **Solo 721** con embeddings (39%)
- **1137 items NON ricercabili** semanticamente!

### Causa
Lo script `knowledge_ingest_mirror7.py` NON genera embeddings.
Serve eseguire lo script separato.

### Soluzione Immediata
```bash
# Genera embeddings per tutti i record mancanti
cd ~/.config/surrealdb/scripts
python3 generate_missing_embeddings.py
```

### Soluzione Permanente
Aggiornare `knowledge_ingest_mirror7.py` per generare embeddings inline:
```python
from sentence_transformers import SentenceTransformer
model = SentenceTransformer('all-MiniLM-L6-v2')

def generate_embedding(text: str) -> List[float]:
    return model.encode(text, normalize_embeddings=True).tolist()
```

---

## 📊 AUDIT TECNICO - Stato Consolidamento

### Codebase Structure
```
Mirror7_build/
├── source/
│   ├── engine/
│   │   ├── Mirror7Engine.cpp/h    # Core DSP - MONOLITHIC
│   │   ├── DialogueSystem.cpp/h   # Hemispheric alternation
│   │   ├── Spatializer.cpp/h      # ITD/ILD/CEFG
│   │   └── NoiseController.cpp/h  # Noise generation
│   ├── plugin/
│   │   └── PluginProcessor.cpp/h  # JUCE wrapper
│   └── gui/
│       └── PluginEditor.cpp/h     # Basic UI
├── modules/                        # From vellutoblu~
│   ├── harmony/                   # HarmonySystem (TO PORT)
│   ├── control/                   # Tempo control
│   ├── core/                      # Modal engine
│   └── dsp/                       # Scheduler, dialogue
├── legacy/                        # From aureonoise_tilde
│   ├── aureo_core_v1/            # Historical C code
│   └── max_external/             # Max/MSP source
└── docs/                          # Documentation
```

### Stato dei Componenti

| Component | Status | Issues | Priority |
|-----------|--------|--------|----------|
| **Mirror7Engine** | ⚠️ Works | Monolithic, burst_w=0 | 🔴 HIGH |
| **DialogueSystem** | ⚠️ Works | Magic numbers undocumented | 🔴 HIGH |
| **Spatializer** | ✅ Good | Minor improvements | 🟡 LOW |
| **NoiseController** | ✅ Good | - | 🟢 OK |
| **HarmonySystem** | ❌ Not integrated | In modules/ | 🟠 MED |
| **GUI** | ⚠️ Basic | Needs coherence viz | 🟠 MED |
| **Tests** | ⚠️ Minimal | Low coverage | 🟠 MED |

### Technical Debt

1. **burst_w hardcoded to 0** - Disables Hawkes burst feature
2. **Magic numbers** - 0.42, 0.55, 0.618 undocumented
3. **Monolithic process()** - 200+ lines single method
4. **No frequency presets** - Missing alpha/theta/gamma modes
5. **No R→L priority** - DialogueSystem symmetric

---

## 📋 TODO LIST - Prioritized

### 🔴 Phase 0: Infrastructure (BEFORE CODING)

- [ ] **0.1** Run `generate_missing_embeddings.py` (~1137 items)
- [ ] **0.2** Verify all 1858 items searchable
- [ ] **0.3** Create GitHub milestone for Phase 1-4
- [ ] **0.4** Set up CI/CD basics (CMake test)

### 🔴 Phase 1: Bug Fix & Core (Week 1)

- [ ] **1.1** Fix `burst_w` - expose as parameter
- [ ] **1.2** Document magic numbers with paper references
- [ ] **1.3** Add R→L priority to DialogueSystem (from Callosal Relay paper)
- [ ] **1.4** Add frequency presets:
  - [ ] `DeepMeditation` - 4-8 Hz theta
  - [ ] `CalmFocus` - 10 Hz alpha
  - [ ] `PeakPerformance` - 40 Hz gamma
  - [ ] `HeartSync` - 0.1 Hz LFO

### 🟠 Phase 2: Feature Port (Week 2-3)

- [ ] **2.1** Port HarmonySystem from modules/harmony/
- [ ] **2.2** Add prime_comb filter
- [ ] **2.3** Add time_quant/time_strict (Fibonacci grid)
- [ ] **2.4** Implement ComplexCorrelation metric

### 🟡 Phase 3: GUI & UX (Week 4)

- [ ] **3.1** Collapsible parameter sections
- [ ] **3.2** Coherence visualizer (L/R correlation)
- [ ] **3.3** Handshake indicator
- [ ] **3.4** Factory presets based on brain states

### 🟢 Phase 4: Quality (Ongoing)

- [ ] **4.1** Increase test coverage to 80%
- [ ] **4.2** GitHub Actions CI/CD
- [ ] **4.3** API documentation

---

## 🤖 MULTI-AGENT ARCHITECTURE

### Objective
Usare agenti multipli con prompt atomici per:
- Risparmiare token (Sonnet vs Opus)
- Parallelizzare task indipendenti
- Aumentare precisione con focus specifico

### Agent Stack

```
┌─────────────────────────────────────────────────────────────┐
│                    ORCHESTRATOR (Opus)                       │
│  - Task decomposition                                        │
│  - Decision making                                           │
│  - Complex reasoning                                         │
└─────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
┌───────────────┐    ┌───────────────┐    ┌───────────────┐
│ CODE_AGENT    │    │ RESEARCH_AGENT│    │ TEST_AGENT    │
│ (Sonnet)      │    │ (Sonnet)      │    │ (Sonnet)      │
│               │    │               │    │               │
│ - Write code  │    │ - Query KB    │    │ - Write tests │
│ - Fix bugs    │    │ - Find papers │    │ - Run tests   │
│ - Refactor    │    │ - Summarize   │    │ - Validate    │
└───────────────┘    └───────────────┘    └───────────────┘
```

### Atomic Prompts (Examples)

#### CODE_AGENT Prompt
```markdown
# Task: Implement frequency preset

## Context
File: source/engine/Mirror7Engine.h
Current: No preset system

## Requirement
Add enum `BrainState` with 4 values:
- DeepMeditation (6 Hz)
- CalmFocus (10 Hz)
- PeakPerformance (40 Hz)
- HeartSync (0.1 Hz LFO)

Add method `void setPreset(BrainState state)`

## Output
Only the code changes, no explanation.
```

#### RESEARCH_AGENT Prompt
```markdown
# Task: Find ITD implementation details

## Query
Search knowledge base for:
- "ITD implementation"
- "interaural time difference algorithm"
- "binaural delay line"

## Output Format
| Paper | Key Formula | Code Reference |
|-------|-------------|----------------|
```

#### TEST_AGENT Prompt
```markdown
# Task: Write test for setPreset

## Target
Mirror7Engine::setPreset(BrainState)

## Test Cases
1. DeepMeditation sets grainRate to 6 Hz
2. CalmFocus sets grainRate to 10 Hz
3. PeakPerformance sets grainRate to 40 Hz
4. HeartSync sets lfoRate to 0.1 Hz

## Framework
Catch2
```

### Token Optimization

| Model | Use Case | Cost/1M tokens |
|-------|----------|----------------|
| **Opus** | Orchestration, complex decisions | ~$15 |
| **Sonnet** | Code generation, simple tasks | ~$3 |
| **Haiku** | Formatting, simple queries | ~$0.25 |

**Strategy:**
1. Opus decompone task in atomic units
2. Sonnet esegue ogni unit
3. Opus valida e integra

---

## 🗂️ FILE STRUCTURE FOR AGENTS

```
Mirror7_build/
├── .agents/
│   ├── prompts/
│   │   ├── code_agent.md
│   │   ├── research_agent.md
│   │   ├── test_agent.md
│   │   └── orchestrator.md
│   ├── tasks/
│   │   ├── pending/           # Task queue
│   │   ├── in_progress/       # Currently assigned
│   │   └── completed/         # Done tasks
│   └── context/
│       ├── codebase_summary.md
│       ├── paper_index.md
│       └── current_state.md
```

---

## 📅 TIMELINE

| Week | Focus | Deliverable |
|------|-------|-------------|
| **0** | Infrastructure | Embeddings complete, CI setup |
| **1** | Phase 1 | Frequency presets, R→L priority |
| **2** | Phase 2a | HarmonySystem integrated |
| **3** | Phase 2b | ComplexCorrelation metric |
| **4** | Phase 3 | GUI improvements |
| **5+** | Phase 4 | Tests, docs, polish |

---

## 🎯 IMMEDIATE ACTIONS

### Right Now
```bash
# 1. Generate missing embeddings
cd ~/.config/surrealdb/scripts
python3 generate_missing_embeddings.py

# 2. Verify completion
curl -s -X POST http://localhost:8000/sql \
  -H "Authorization: Basic cm9vdDpyb290" \
  -H "surreal-ns: research" \
  -H "surreal-db: knowledge" \
  -d "SELECT count() FROM knowledge WHERE embedding != NONE GROUP ALL;"
```

### This Session
1. ✅ Create this TODO document
2. ⬜ Push to GitHub
3. ⬜ Generate embeddings
4. ⬜ Start Phase 1.4 (frequency presets)

### Next Session
1. Implement R→L DialogueSystem priority
2. Document magic numbers
3. Fix burst_w parameter
