#pragma once

extern "C" {
#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include "ext_systhread.h"
}

#include <array>
#include <atomic>
#include <cstdint>

#include "aureo_core/aureo_binaural.hpp"
#include "aureo_core/aureo_config.hpp"
#include "aureo_core/aureo_field_adsr.hpp"
#include "aureo_core/aureo_math.hpp"
#include "aureo_core/aureo_noise.hpp"
#include "aureo_core/aureo_pinna.hpp"
#include "aureo_core/aureo_ring.hpp"
#include "aureo_core/aureo_rng.hpp"
#include "aureo_core/aureo_stoch.hpp"
#include "aureo_core/aureo_weyl.hpp"
#include "aureonoise/spatial/cefg.hpp"
#include "aureonoise3beta/core/modal_engine.hpp"
#include "../harmony/harmony.hpp"
#include "../control/tempo.hpp"

static constexpr int kSpatialProfileLegacy = 0;
static constexpr int kSpatialProfileUteroidPhi = 1;
static constexpr int kSpatialCefgMaxTaps = 16;

struct SpatialCefgTapState {
  double base_delay = 0.0;
  double itd_offset_l = 0.0;
  double itd_offset_r = 0.0;
  double base_gain_l = 0.0;
  double base_gain_r = 0.0;
  double gain_l = 0.0;
  double gain_r = 0.0;
  double lp_alpha = 0.0;
  double lp_state_l = 0.0;
  double lp_state_r = 0.0;
  double mix = 0.0;
};

struct SpatialCefgGrainState {
  int tap_count = 0;
  std::array<SpatialCefgTapState, kSpatialCefgMaxTaps> taps{};
};

struct t_aureonoise {
  t_pxobject ob;

  struct GrainReport {
    uint64_t index = 0;
    double   timestamp_sec = 0.0;
    double   gap_samples = 0.0;
    double   dur_samples = 0.0;
    double   amp = 0.0;
    double   pan = 0.0;
    double   itd_samples = 0.0;
    double   ild_db = 0.0;
    double   binaural_azimuth_deg = 0.0;
    double   binaural_focus = 0.0;
    double   binaural_crossfeed = 0.0;
    double   phi_u = 0.0;
    double   s2_u = 0.0;
    double   pl_u = 0.0;
    double   rng4 = 0.0;
    double   rng5 = 0.0;
    double   rng6 = 0.0;
    double   lattice_u = 0.5;
    double   lattice_v = 0.0;
    double   hawkes_lambda = 0.0;
    bool     thermo_on = false;
    bool     lattice_on = false;
    bool     burst_on = false;
    double   dialogue_coherence = 0.0;
    double   dialogue_reply = 0.0;
  };

  double p_rate = 8.0;
  double p_baselen_ms = 120.0;
  double p_len_phi = 0.8;
  double p_width = 1.0;
  double p_itd_us = 700.0;
  double p_ild_db = 20.0;
  double p_spat_min_deg = 12.0;
  double p_spat_min_ms = 35.0;
  double p_spat_ipd = 0.6;
  double p_spat_shadow = 0.7;
  double p_spat_external = 0.0;
  long   p_spat_mirror = 0;   // new: full spatial pan mirror (flip hemisphere)
  double p_env_attack = 0.18;
  double p_env_decay = 0.28;
  double p_env_sustain = 0.55;
  double p_env_release = 0.30;
  double p_hemis_coupling = 0.6;
  long   p_noise_color = static_cast<long>(aureo::NoiseColor::Pink);
  double p_color_amt = 0.65;
  double p_noise_mode = static_cast<double>(aureo::NoiseColor::Pink);
  double p_aureo_mix = 0.5;
  double p_aureo_decay = 0.3;
  double p_aureo_stride = 1.0;
  double p_aureo_harmonics = 12.0;
  long   p_aureo_velvet = 1; // on/off: include velvet residual in Aureo mode
  double p_quantum_mix = 0.55;
  double p_quantum_detail = 0.7;
  double p_quantum_stride = 1.0;
  double p_quantum_base = 220.0;
  double p_quantum_velvet = 0.12;
  double p_vhs_wow = 0.35;
  double p_vhs_flutter = 0.25;
  double p_glitch_mix = 0.5;
  double p_srcrush_amt = 0.2;
  double p_bitcrush_amt = 0.15;
  double p_burst_floor = 0.35;
  double p_burst_phi_mix = 0.6;
  long   p_dialogue_on = 1;
  double p_dialogue_strength = 0.6;
  double p_dialogue_memory = 0.5;
  double p_dialogue_phi_mix = 0.75;
  long   p_modal_on = 0;
  long   p_modal_mirror = 0;
  double p_modal_mix = 0.35;
  double p_modal_decay = 0.5;
  long   p_modal_preset = 1;
  long   p_seed = 20251010;
  long   p_sync_enable = 0;
  double p_sync_division = 1.0;
  double p_sync_slew = 0.5;
  long   p_pinna_on = 0;
  long   p_prime_comb = 0;
  long   p_time_quant = 0;
  long   p_time_strict = 0;
  long   p_phi_pan = 0;
  double p_pinna_depth = 12.0;
  long   p_phi_mode = 1;
  double p_phi_ratio = 1.6180339887;
  double p_phi_head_b = 0.0875;
  double p_phi_distance = 0.12;
  double p_phi_elev = 0.6;
  double p_phi_elev_notch = 0.6;
  double p_phi_torso_mix = 0.12;
  double p_phi_torso_ms = 1.15;
  double p_phi_torso_hp_hz = 500.0;
  long   p_spatial_profile = kSpatialProfileLegacy;
  double p_cefg_gain = 1.1;
  double p_cefg_mix = 3.4;
  long   p_cefg_debug = 0;

#if AUREO_THERMO_LATTICE
  long p_thermo = 1;
  long p_lattice = 1;
  long p_burst = AUREO_BURST_HAWKES ? 1 : 0;
  double p_T = 0.45;
  double p_lat_rate = 250.0;
  double p_lat_eps = 0.18;
  double p_lat_gamma = 1.40;
  double p_lat_sigma = 0.06;
  long   p_lat_x = 8;
  long   p_lat_y = 8;
  long   p_lat_z = 4;
#endif

  double sr = 44100.0;
  aureo::RNG rng;
  aureo::Weyl w_phi;
  aureo::Weyl w_s2;
  aureo::Weyl w_pl;
  aureo::AureoNoiseEngine noise;
  aureonoise::HarmonyState harmony;
  aureonoise::control::TempoState tempo;
  aureo::FieldState field;
  aureo::RingState ring;
  aureo::Pinna pinna;
  aureo::Biquad tone;
  aureo::Biquad pinna_notchL;
  aureo::Biquad pinna_notchR;

  aureonoise::spatial::cefg::Profile cefg_profile;
  bool cefg_profile_loaded = false;
  std::array<SpatialCefgGrainState, aureo::kMaxGrains> cefg_grains{};

  double pinna_mix = 0.0;
  double pinna_mix_target = 0.0;
  double pinna_depth = 12.0;
  double pinna_freq_left = 7200.0;
  double pinna_freq_right = 8200.0;
  double pinna_q = 5.0;
  double pinna_min_freq = 400.0;
  double pinna_notch_base = 7800.0;
  double pinna_notch_spread = 2200.0;
  bool   pinna_enabled = false;
  bool   pinna_filters_dirty = true;

  double phi_last_mag = 0.0;
  int    phi_last_sign = 1;

  double lfo_wow_phase = 0.0;
  double lfo_flut_phase = 0.0;

  uint64_t sample_counter = 0;
  int samples_to_next = 0;
  int gap_elapsed = 0;
  std::array<int, 256> primes{};
  int primes_count = 0;
  std::array<aureo::Grain, aureo::kMaxGrains> grains{};

  void* out_info = nullptr;
  t_systhread_mutex report_mu = nullptr;
  static constexpr size_t kReportCapacity = 512;
  std::array<GrainReport, kReportCapacity> report_log{};
  size_t report_head = 0;
  size_t report_count = 0;
  std::atomic<uint64_t> report_total{0};
  std::atomic<uint64_t> report_dropped{0};

  double prev_pan = 0.0;
  double prev_itd = 0.0;
  double prev_ild = 0.0;
  double last_gap_samples = 0.0;
  double last_dur_samples = 0.0;
  double last_handshake_time = 0.0;

  struct DialogueState {
    double last_sign = 0.0;
    double last_mag = 0.0;
    double last_amp = 1.0;
    double last_dur = 0.0;
    double last_gap = 0.0;
    double last_timestamp = 0.0;
    double coherence = 1.0;
    double coherence_accum = 0.0;
    uint64_t coherence_samples = 0;
    uint64_t utterances = 0;
    uint64_t handshakes = 0;
    uint64_t failures = 0;
    double pan_sum = 0.0;
    double pan_abs_sum = 0.0;
    uint64_t pan_samples = 0;
    int fib_index = 1;
    int fib_direction = 1;
    double fib_ratio_target = 1.0;
    double fib_score = 0.0;
    std::array<double, 3> fib_gap_queue{};
    uint32_t fib_gap_queue_size = 0;
    uint32_t fib_gap_queue_pos = 0;
  } dialogue;

  aureonoise::modal::Engine modal;
  double modal_energy = 0.0;

#if AUREO_THERMO_LATTICE
  aureo::OU ou_pan;
  aureo::OU ou_itd;
  aureo::OU ou_amp;
  aureo::OU ou_rate;
  aureo::Lattice lat;
  double lat_phase = 0.0;
  t_systhread_mutex lat_mu = nullptr;
  double lat_last_v = 0.0;
#if AUREO_BURST_HAWKES
  aureo::Hawkes hawkes;
#endif
#endif

  static constexpr int kExternalDelaySize = 2048;
  static constexpr int kExternalDelayMask = kExternalDelaySize - 1;
  std::array<double, kExternalDelaySize> ext_delayL{};
  std::array<double, kExternalDelaySize> ext_delayR{};
  uint32_t ext_write_idx = 0;
  double ext_lpL = 0.0;
  double ext_lpR = 0.0;
};

#if AUREO_THERMO_LATTICE
void aureonoise_lattice_safe_resize(t_aureonoise* x, int X, int Y, int Z);
#endif

inline void aureonoise_update_pinna_state(t_aureonoise* x)
{
  x->pinna_depth = aureo::clamp(x->p_pinna_depth, 0.0, 24.0);
  x->pinna_mix_target = (x->p_pinna_on != 0) ? 1.0 : 0.0;
}

inline double aureonoise_pinna_weight(const t_aureonoise* x)
{
  return aureo::clamp01(x->pinna_mix) * aureo::clamp01(x->pinna_depth / 24.0);
}

inline void aureonoise_update_pinna_filters(t_aureonoise* x)
{
  double sr = x->sr;
  if (sr <= 0.0) sr = 44100.0;
  if (sr < 4000.0) sr = 4000.0;

  const double width = aureo::clamp(x->p_width, 0.0, 1.0);
  const double base = x->pinna_notch_base;
  const double spread = x->pinna_notch_spread;
  x->pinna_freq_left = aureo::clamp(base - 0.5 * spread * width, 1200.0, 16000.0);
  x->pinna_freq_right = aureo::clamp(base + 0.5 * spread * width, 1200.0, 16000.0);

  x->pinna_notchL.setNotch(sr, x->pinna_freq_left, x->pinna_q, x->pinna_min_freq);
  x->pinna_notchR.setNotch(sr, x->pinna_freq_right, x->pinna_q, x->pinna_min_freq);
  x->pinna_filters_dirty = false;
}

inline void aureonoise_mark_pinna_dirty(t_aureonoise* x)
{
  x->pinna_filters_dirty = true;
}

void aureonoise_setup_attributes(t_class* c);
int  pick_prime_in_range(t_aureonoise* x, int lo, int hi, double u);
void make_small_primes(t_aureonoise* x);
void aureonoise_report(t_aureonoise* x);

void aureonoise_refresh_harmony(t_aureonoise* x);
void aureonoise_update_sync(t_aureonoise* x, double block_dt);
void aureonoise_set_tempo(t_aureonoise* x, double bpm);
void aureonoise_update_modal(t_aureonoise* x);

t_max_err set_rate(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_baselen(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lenphi(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_width(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_prime_comb(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_time_quant(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_time_strict(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_sync_enable(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_sync_division(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_sync_slew(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_pan(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_itd(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_ild(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_min_deg(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_min_ms(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_ipd(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_shadow(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_external(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spat_mirror(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_attack(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_decay(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_sustain(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_env_release(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_hemis_coupling(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_noise_color(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_color_amt(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_noise_mode(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_aureo_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_aureo_decay(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_aureo_stride(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_aureo_harmonics(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_quantum_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_quantum_detail(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_quantum_stride(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_quantum_base(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_quantum_velvet(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_wow(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_vhs_flutter(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_glitch_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_srcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_bitcrush(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_burst_floor(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_burst_phi_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_dialogue_on(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_dialogue_strength(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_dialogue_memory(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_dialogue_phi_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_modal_on(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_modal_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_modal_decay(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_modal_preset(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_modal_mirror(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_seed(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_pinna_on(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_pinna_depth(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_mode(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_ratio(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_head_b(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_distance(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_elev(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_elev_notch(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_torso_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_torso_ms(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_phi_torso_hp(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_spatial_profile(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_cefg_gain(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_cefg_mix(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_cefg_debug(t_aureonoise* x, void*, long ac, t_atom* av);
#if AUREO_THERMO_LATTICE
t_max_err set_thermo(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lattice(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_burst(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_T(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_rate(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_eps(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_gamma(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_sigma(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_x(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_y(t_aureonoise* x, void*, long ac, t_atom* av);
t_max_err set_lat_z(t_aureonoise* x, void*, long ac, t_atom* av);
#endif
