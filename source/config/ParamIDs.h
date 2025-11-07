#pragma once

namespace ParamIDs {
// Core output
static constexpr const char* outGain      = "out_gain";

// Panning / stereo
static constexpr const char* pan          = "pan";
static constexpr const char* width        = "width";
static constexpr const char* autopan      = "autopan";
static constexpr const char* autopanSp    = "autopan_speed";

// ITD/ILD and spatial guard
static constexpr const char* itdUs        = "itd_us";
static constexpr const char* ildDb        = "ild_db";
static constexpr const char* spatMinDeg   = "spat_min_deg";
static constexpr const char* spatMinMs    = "spat_min_ms";
static constexpr const char* spatIPD      = "spat_ipd";
static constexpr const char* spatShadow   = "spat_shadow";
static constexpr const char* spatMirror   = "spat_mirror";

// Grain timing
static constexpr const char* rateHz       = "rate_hz";
static constexpr const char* baseMs       = "base_ms";
static constexpr const char* lenPhi       = "len_phi";
static constexpr const char* hemisCoupling= "hemis_coupling";

// Envelope
static constexpr const char* envA         = "env_attack";
static constexpr const char* envD         = "env_decay";
static constexpr const char* envS         = "env_sustain";
static constexpr const char* envR         = "env_release";

// Dialogue
static constexpr const char* dialogOn     = "dialogue_on";
static constexpr const char* dialogStr    = "dialogue_strength";
static constexpr const char* dialogMem    = "dialogue_memory";
static constexpr const char* dialogPhi    = "dialogue_phi_mix";

// Glitch/Crush
static constexpr const char* glitchMix    = "glitch_mix";
static constexpr const char* srCrush      = "sr_crush";
static constexpr const char* bitCrush     = "bit_crush";

// Noise
static constexpr const char* noiseMode    = "noise_mode";     // choice
static constexpr const char* noiseColor   = "noise_color";    // choice
static constexpr const char* colorAmt     = "color_amt";
static constexpr const char* velvetAmt    = "velvet_amt";
// Aureo
static constexpr const char* aureoMix     = "aureo_mix";
static constexpr const char* aureoDecay   = "aureo_decay";
static constexpr const char* aureoStride  = "aureo_stride";
static constexpr const char* aureoHarmo   = "aureo_harmonics";
static constexpr const char* aureoVelvet  = "aureo_velvet";   // bool
// Quantum
static constexpr const char* qMix         = "quantum_mix";
static constexpr const char* qDetail      = "quantum_detail";
static constexpr const char* qStride      = "quantum_stride";
static constexpr const char* qBase        = "quantum_base";
static constexpr const char* qVelvet      = "quantum_velvet";

// VHS
static constexpr const char* vhsWow       = "vhs_wow";
static constexpr const char* vhsFlutter   = "vhs_flutter";

// Phi
static constexpr const char* phiPan       = "phi_pan";        // bool
static constexpr const char* phiMode      = "phi_mode";       // bool
static constexpr const char* phiRatio     = "phi_ratio";
static constexpr const char* phiHeadB     = "phi_head_b";
static constexpr const char* phiDistance  = "phi_distance";
static constexpr const char* phiElev      = "phi_elev";
static constexpr const char* phiElevNotch = "phi_elev_notch";
static constexpr const char* phiTorsoMix  = "phi_torso_mix";
static constexpr const char* phiTorsoMs   = "phi_torso_ms";
static constexpr const char* phiTorsoHp   = "phi_torso_hp_hz";

// Modal
static constexpr const char* modalOn      = "modal_on";
static constexpr const char* modalMirror  = "modal_mirror";
static constexpr const char* modalMix     = "modal_mix";
static constexpr const char* modalDecay   = "modal_decay";
static constexpr const char* modalPreset  = "modal_preset";   // choice

// Pinna
static constexpr const char* pinnaOn      = "pinna_on";
static constexpr const char* pinnaDepth   = "pinna_depth";

// Spatial v2 (CEFG)
static constexpr const char* spatialProfile = "spatial_profile"; // choice
static constexpr const char* cefgGain    = "cefg_gain";
static constexpr const char* cefgMix     = "cefg_mix";

// Sync
static constexpr const char* syncEnable  = "sync_enable";
static constexpr const char* syncDivision= "sync_division";
static constexpr const char* syncSlew    = "sync_slew";

// Quality
static constexpr const char* osFactor    = "oversample";
} // namespace ParamIDs
