# Aureonoise / Mirror7 Ecosystem Inventory

Ultimo aggiornamento: 2026-01-09

## 🎯 Progetti Principali (Attivi)

### Mirror7 JUCE Plugin
| Location | Stato | Note |
|----------|-------|------|
| `/Users/alessioivoycazzaniga/Mirror7_build` | ✅ MAIN | Copia di lavoro principale |
| `/Users/alessioivoycazzaniga/Documents/GitHub/Mirror7_build` | 📦 Backup | Copia originale |
| `Library/Audio/Plug-Ins/VST3/Mirror7.vst3` | 🔌 Installato | Plugin VST3 |
| `Library/Audio/Plug-Ins/Components/Mirror7.component` | 🔌 Installato | Plugin AU |

### Aureonoise Monorepo
| Location | Descrizione |
|----------|-------------|
| `/Users/alessioivoycazzaniga/aureonoise_tilde` | Monorepo principale (home) |
| `/Users/alessioivoycazzaniga/Documents/GitHub/aureonoise_tilde` | Copia GitHub |
| `GitHub: Alemusica/aureonoise_tilde` | 🌐 Repo online (public) |

### Aureo Factory (DSP Core + ML + Plugin)
| Location | Descrizione |
|----------|-------------|
| `/Users/alessioivoycazzaniga/alessioivoycazzaniga/aureo-factory` | Libreria DSP Python/C++, ML lab, plugin |

---

## 🔌 Max Externals Installati

### Aureonoise Family
```
/Users/alessioivoycazzaniga/externals/
├── aureonoise~.mxo          # Versione base
├── aureonoise2beta~.mxo     # Beta 2
├── aureonoise3beta~.mxo     # Beta 3
└── aureonoise_beta4~.mxo    # Beta 4
```

### Mirror Family
```
/Users/alessioivoycazzaniga/externals/
├── beta4_mirror~.mxo        # Mirror base
├── beta4_mirror^2~.mxo      # Mirror v2
├── beta4_mirror^3~.mxo      # Mirror v3
├── mirror7~.mxo             # Mirror 7
├── mirror8~.mxo             # Mirror 8
├── mirror8_1~.mxo           # Mirror 8.1
├── mirror8_2~.mxo           # Mirror 8.2
└── mirror8_3~.mxo           # Mirror 8.3
```

### Max Packages
```
/Users/alessioivoycazzaniga/Documents/Max 9/Packages/
├── aureonoise-beta2/
├── aureonoise2beta/
├── aureonoise3beta/
├── aureonoise4beta/
├── aureonoise_beta4/
└── aureonoise_beta5/
```

---

## 📁 Progetti GitHub Locali (Correlati)

### Phi/Golden Family
| Progetto | Path | GitHub Online |
|----------|------|---------------|
| PhiVerb | `Documents/GitHub/PhiVerb` | ✅ `Alemusica/phiverb` |
| PhiVerb-dwm | `Documents/GitHub/PhiVerb-dwm` | ❓ |
| phi-branch-synth | `Documents/GitHub/phi-branch-synth` | ❓ |
| golden-helix-studio | `Documents/GitHub/golden-helix-studio` | ✅ Private |
| GoldenSynth | `Documents/GitHub/GoldenSynth` | ❓ |
| golden-ratio-toolkit-v3 | `Documents/GitHub/golden-ratio-toolkit-v3` | ❓ |

### Mirror/Aureo Family
| Progetto | Path | Note |
|----------|------|------|
| aureonoise_tilde | `Documents/GitHub/aureonoise_tilde` | Monorepo |
| beta4_mirror_repo | `Documents/GitHub/beta4_mirror_repo` | Legacy |
| mirror_backup | `Documents/GitHub/mirror_backup` | Backup |
| waterchorus/mirror7 | `Documents/GitHub/waterchorus/mirror7` | Legacy in altro repo |

---

## 🔄 Struttura Monorepo aureonoise_tilde

```
aureonoise_tilde/
├── source/
│   ├── aureo_core/          # Header DSP condivisi (✅ usato da Mirror7)
│   ├── uterine/             # Modulo uterino
│   └── vellutoblu~/         # Versione rinominata di aureonoise3beta
│       ├── beta7_tools/     # Spatial assets, phi_model (✅ usato da Mirror7)
│       ├── harmony/         # Sistema armonico (❌ NON in Mirror7)
│       ├── dsp/             # DSP core
│       └── core/            # Modal engine
├── beta4_mirror/            # Versione Max di Mirror (legacy)
├── aureonoise_py/           # Binding Python
└── docs/
```

---

## 📊 Feature Matrix

| Feature | aureonoise3beta~ | Mirror7 JUCE | Note |
|---------|-----------------|--------------|------|
| Granular Engine | ✅ | ✅ | Identico |
| DialogueSystem | ✅ | ✅ | Portato |
| Spatializer (Φ-mode) | ✅ | ✅ | Portato |
| NoiseController | ✅ | ✅ | Portato |
| **HarmonySystem** | ✅ | ❌ | Da portare |
| **Prime Comb** | ✅ | ❌ | Da portare |
| **Time Quantization** | ✅ | ❌ | Da portare |
| **Burst (Hawkes)** | ✅ | ⚠️ Disabilitato | burst_w=0 |
| Modal Engine | ✅ | ✅ | Portato |
| Pinna Notch | ✅ | ✅ | Portato |
| CEFG Early Reflections | ✅ | ✅ | Portato |
| Host Sync | ✅ | ✅ | Portato |
| Preset System | N/A | ✅ | Solo JUCE |
| GUI | N/A | ⚠️ Basic | Da migliorare |

---

## 🎯 Roadmap Integrazione

### Fase 1: Consolidamento
- [ ] Push Mirror7_build su GitHub
- [ ] Sincronizzare con aureonoise_tilde
- [ ] Documentare dipendenze

### Fase 2: Port Feature Mancanti
- [ ] HarmonySystem (prime_comb, time_quant)
- [ ] Attivare burst (Hawkes process)
- [ ] Time strict quantization

### Fase 3: Miglioramenti
- [ ] GUI moderna con sezioni
- [ ] Factory presets
- [ ] CI/CD

---

## 📝 Note

- **aureo-factory** contiene librerie DSP avanzate (Python/C++) e modelli ML per acustica uterina
- **PhiVerb** è un reverb basato su φ, separato ma con filosofia simile
- **golden-helix-studio** potrebbe contenere GUI/UX patterns riutilizzabili
