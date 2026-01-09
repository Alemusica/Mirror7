#include "aureonoise3beta/dsp/dsp_api.hpp"
#include "aureonoise3beta/dsp/scheduler.hpp"
#include "aureonoise3beta/control/tempo.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <cstdint>
#include <random>
#ifdef __APPLE__
#include <Security/SecRandom.h>
#include <mach/mach_time.h>
#endif

#ifdef __SSE2__
#include <xmmintrin.h>
#endif

#include "aureonoise3beta/dsp/spatial_utils.hpp"
#include "aureonoise3beta/dsp/dialogue.hpp"
#include "aureonoise/spatial/phi_model.hpp"

namespace phi = aureonoise::spatial::phi;

namespace {

#ifdef AUREONOISE_SPATIAL_ASSETS_DIR
static std::string default_cefg_profile_path()
{
  return std::string(AUREONOISE_SPATIAL_ASSETS_DIR) + "/cefg/uteroid_phi.cefg";
}
#else
static std::string default_cefg_profile_path()
{
  return std::string("spatial_assets/cefg/uteroid_phi.cefg");
}
#endif

inline std::array<double, 3> direction_from_angles(double azimuth_deg,
                                                   double elevation_deg)
{
  const double az = azimuth_deg * (aureo::kPi / 180.0);
  const double el = elevation_deg * (aureo::kPi / 180.0);
  const double cosEl = std::cos(el);
  return {cosEl * std::sin(az),
          std::sin(el),
          cosEl * std::cos(az)};
}

inline double ring_read_delay(const aureo::RingState& ring,
                              uint32_t writeIndex,
                              double delaySamples)
{
  const double clampedDelay = std::max(0.0, delaySamples);
  const double pos = static_cast<double>(writeIndex) - clampedDelay;
  return aureo::lagrange3(ring.data.data(), pos);
}

inline void cefg_clear_grain_state(SpatialCefgGrainState& state)
{
  state.tap_count = 0;
  for (auto& tap : state.taps) {
    tap.lp_state_l = 0.0;
    tap.lp_state_r = 0.0;
    tap.mix = 0.0;
    tap.gain_l = 0.0;
    tap.gain_r = 0.0;
  }
}

inline void cefg_clear_all(t_aureonoise* x)
{
  if (!x) return;
  for (auto& state : x->cefg_grains) {
    cefg_clear_grain_state(state);
  }
}

// Simple entropy pool (SecRandomCopyBytes + jitter fallback)
struct EntropyPool {
  std::array<uint8_t, 8192> buf{};
  size_t pos{0};
  size_t len{0};
  std::random_device rd;

  void refill()
  {
#ifdef __APPLE__
    size_t want = buf.size();
    if (SecRandomCopyBytes(kSecRandomDefault, want, buf.data()) == errSecSuccess) {
      // Lightly stir with time jitter to decorrelate across calls
      uint64_t t = mach_absolute_time();
      for (size_t i = 0; i < want; ++i) buf[i] ^= static_cast<uint8_t>((t >> (i & 15)) & 0xFFu);
      pos = 0; len = want; return;
    }
#endif
    // Fallback: fill with std::random_device + simple xorshift
    uint32_t s = rd();
    auto xorshift32 = [&s](){ s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; };
    for (size_t i = 0; i < buf.size(); i += 4) {
      uint32_t r = xorshift32();
      buf[i+0] = static_cast<uint8_t>((r) & 0xFFu);
      buf[i+1] = static_cast<uint8_t>((r >> 8) & 0xFFu);
      buf[i+2] = static_cast<uint8_t>((r >> 16) & 0xFFu);
      buf[i+3] = static_cast<uint8_t>((r >> 24) & 0xFFu);
    }
    pos = 0; len = buf.size();
  }

  inline uint8_t next_byte()
  {
    if (pos >= len) refill();
    return buf[pos++];
  }

  inline double unit_pm1()
  {
    // Map two bytes -> [-1,1]
    uint16_t v = static_cast<uint16_t>(next_byte()) | (static_cast<uint16_t>(next_byte()) << 8);
    const double u = static_cast<double>(v) / 65535.0; // [0,1]
    return 2.0 * u - 1.0;
  }
};

static EntropyPool g_entropy;

inline void cefg_prepare_grain(t_aureonoise* x,
                               size_t grainIndex,
                               const phi::Geometry& geom,
                               double sampleRate,
                               double phi_distance,
                               double pitch_rad,
                               double ux_pitch,
                               double uy_pitch,
                               double uz_pitch)
{
  if (!x) return;
  auto& state = x->cefg_grains[grainIndex];
  cefg_clear_grain_state(state);
  if (!x->cefg_profile_loaded) return;
  if (x->p_spatial_profile != kSpatialProfileUteroidPhi) return;

  const double mixScale = aureo::clamp(x->p_cefg_mix, 0.0, 4.0);
  if (mixScale <= 1.0e-6) return;

  const double az_deg = std::atan2(ux_pitch, uz_pitch) * (180.0 / aureo::kPi);
  const double el_deg = std::atan2(uy_pitch, std::hypot(ux_pitch, uz_pitch)) * (180.0 / aureo::kPi);
  const double distance_m = 0.8 + 3.2 * aureo::clamp01(phi_distance);
  const double log_distance = std::log10(std::max(0.12, distance_m));
  const auto* cell = x->cefg_profile.find(az_deg, el_deg, log_distance);
  if (!cell) return;

  const double gainScale = aureo::clamp(x->p_cefg_gain, 0.0, 4.0);
  const double fc_base = 2200.0 + 1600.0 * (1.0 - aureo::clamp01(phi_distance));
  const double alpha_base = std::exp(-2.0 * aureo::kPi * fc_base / sampleRate);

  int tapCount = 0;
  for (const auto& ref : cell->reflections) {
    if (tapCount >= kSpatialCefgMaxTaps) break;
    const double baseGain = ref.gain;
    const double scaledGain = baseGain * gainScale;
    if (std::abs(scaledGain) < 1.0e-5) continue;

    const auto dirRef = direction_from_angles(ref.azimuth_deg, ref.elevation_deg);
    const auto refHead = phi::compute_head_result_from_dir(geom, sampleRate, dirRef, pitch_rad, phi_distance);
    const double refPan = aureo::clamp(refHead.lateral, -1.0, 1.0);
    double gainL, gainR;
    aureo::pan_equal_power(refPan, 1.0, gainL, gainR);
    const double ildPush = 1.0 + 0.18 * refHead.lateral;

    auto& tap = state.taps[tapCount++];
    tap.base_delay = std::max(0.0, ref.delay_sec * sampleRate);
    tap.itd_offset_l = (refHead.itd_samples < 0.0) ? -refHead.itd_samples : 0.0;
    tap.itd_offset_r = (refHead.itd_samples > 0.0) ? refHead.itd_samples : 0.0;
    tap.base_gain_l = baseGain * gainL / ildPush;
    tap.base_gain_r = baseGain * gainR * ildPush;
    tap.gain_l = tap.base_gain_l * gainScale;
    tap.gain_r = tap.base_gain_r * gainScale;
    tap.lp_alpha = aureo::clamp(alpha_base, 0.0, 0.9995);
    tap.lp_state_l = 0.0;
    tap.lp_state_r = 0.0;
    tap.mix = mixScale;
  }
  state.tap_count = tapCount;
  if (x->p_cefg_debug && tapCount > 0) {
    object_post(reinterpret_cast<t_object*>(x),
                "[CEFG] grain %zu taps=%d mix=%.2f gain=%.2f dist=%.2f",
                grainIndex, tapCount, x->p_cefg_mix, x->p_cefg_gain, distance_m);
  }
}

static void aureonoise_log_grain_event(t_aureonoise* x, const t_aureonoise::GrainReport& report)
{
  if (!x) return;
  bool stored = false;
  if (x->report_mu) {
    if (systhread_mutex_trylock(x->report_mu) == 0) {
      const size_t cap = x->report_log.size();
      if (cap > 0) {
        x->report_log[x->report_head] = report;
        x->report_head = (x->report_head + 1) % cap;
        if (x->report_count < cap) ++x->report_count;
      }
      stored = true;
      systhread_mutex_unlock(x->report_mu);
    }
  } else {
    const size_t cap = x->report_log.size();
    if (cap > 0) {
      x->report_log[x->report_head] = report;
      x->report_head = (x->report_head + 1) % cap;
      if (x->report_count < cap) ++x->report_count;
    }
    stored = true;
  }

  if (!stored) {
    x->report_dropped.fetch_add(1, std::memory_order_relaxed);
  }
}

static inline std::array<double, 3> rotate_x(const std::array<double, 3>& v, double angle)
{
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  return {v[0],
          c * v[1] - s * v[2],
          s * v[1] + c * v[2]};
}

static inline double aureonoise_allpass(double in, double a, double& z)
{
  const double y = -a * in + z;
  z = in + a * y;
  return y;
}

static inline double aureonoise_shadow_lp(double in, double a, double& z)
{
  const double y = (1.0 - a) * in + a * z;
  z = y;
  return y;
}

static bool aureonoise_poisson_enforce(const t_aureonoise* x, double& pan, double minPan, double minSamples)
{
  pan = aureo::clamp(pan, -1.0, 1.0);
  minPan = std::max(0.0, minPan);
  minSamples = std::max(0.0, minSamples);
  if (minPan <= 0.0 || x->grains.empty()) return true;

  double guard = minPan;
  const double relax = 0.82;
  for (int attempt = 0; attempt < 8; ++attempt) {
    bool conflict = false;
    for (const auto& g : x->grains) {
      if (!g.on) continue;
      if (minSamples > 0.0 && static_cast<double>(g.age) > minSamples) continue;
      const double dist = std::abs(pan - g.pan);
      if (dist < guard) {
        conflict = true;
        const double sign = (pan >= g.pan) ? 1.0 : -1.0;
        pan = aureo::clamp(g.pan + sign * guard, -1.0, 1.0);
        break;
      }
    }
    if (!conflict) return true;
    guard *= relax;
  }

  for (const auto& g : x->grains) {
    if (!g.on) continue;
    if (minSamples > 0.0 && static_cast<double>(g.age) > minSamples) continue;
    if (std::abs(pan - g.pan) < 0.5 * minPan) return false;
  }

  pan = aureo::clamp(pan, -1.0, 1.0);
  return true;
}

} // namespace

void aureonoise_spatial_profile_refresh(t_aureonoise* x)
{
  if (!x) return;
  cefg_clear_all(x);
  x->cefg_profile_loaded = false;

  if (x->p_spatial_profile != kSpatialProfileUteroidPhi) return;

  const std::string path = default_cefg_profile_path();
  if (!x->cefg_profile.load(path)) {
    object_post(reinterpret_cast<t_object*>(x),
                "[aureonoise] impossibile caricare profilo CEFG (%s)", path.c_str());
    return;
  }

  x->cefg_profile_loaded = true;
  aureonoise_spatial_profile_params_changed(x);
  object_post(reinterpret_cast<t_object*>(x),
              "[aureonoise] profilo CEFG attivo: %s", path.c_str());
}

void aureonoise_spatial_profile_params_changed(t_aureonoise* x)
{
  if (!x) return;
  x->p_cefg_gain = aureo::clamp(x->p_cefg_gain, 0.0, 4.0);
  x->p_cefg_mix = aureo::clamp(x->p_cefg_mix, 0.0, 4.0);
  if (!x->cefg_profile_loaded) return;

  for (auto& state : x->cefg_grains) {
    for (int i = 0; i < state.tap_count; ++i) {
      auto& tap = state.taps[i];
      tap.gain_l = tap.base_gain_l * x->p_cefg_gain;
      tap.gain_r = tap.base_gain_r * x->p_cefg_gain;
      tap.mix = x->p_cefg_mix;
    }
  }
}

void aureonoise_set_tempo(t_aureonoise* x, double bpm)
{
  aureonoise::control::tempo_set_bpm(x, bpm);
}

void aureonoise_update_sync(t_aureonoise* x, double block_dt)
{
  aureonoise::control::tempo_update(x, block_dt);
}

void aureonoise_update_modal(t_aureonoise* x)
{
  if (!x) return;
  const double sr = (x->sr > 0.0) ? x->sr : 44100.0;
  x->modal.setSampleRate(sr);
  x->modal.configure(static_cast<int>(x->p_modal_preset), aureo::clamp01(x->p_modal_decay));
  if (!x->p_modal_on) {
    x->modal.reset();
    x->modal_energy = 0.0;
  }
}

void aureonoise_reset_reports(t_aureonoise* x)
{
  if (!x) return;
  if (x->report_mu) systhread_mutex_lock(x->report_mu);
  x->report_head = 0;
  x->report_count = 0;
  x->report_total.store(0, std::memory_order_relaxed);
  x->report_dropped.store(0, std::memory_order_relaxed);
  for (auto& entry : x->report_log) entry = {};
  if (x->report_mu) systhread_mutex_unlock(x->report_mu);
}

void aureonoise_prepare_instance(t_aureonoise* x)
{
  if (!x) return;
  aureonoise::scheduler::reset_sequences(x);
  make_small_primes(x);
  x->ring.clear();
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  x->lfo_wow_phase = 0.0;
  x->lfo_flut_phase = 0.0;
  for (auto& g : x->grains) g = {};
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
  x->phi_last_mag = 0.0;
  x->phi_last_sign = 1;
  x->last_handshake_time = 0.0;
  x->dialogue = {};
  x->dialogue.coherence = 1.0;

  x->modal.reset();
  x->modal_energy = 0.0;
  aureonoise_update_modal(x);

  x->pinna.width = x->p_width;
  x->pinna.itd_us = x->p_itd_us;
  x->pinna.ild_db = x->p_ild_db;
  x->pinna_enabled = (x->p_pinna_on != 0);
  aureonoise_update_pinna_state(x);
  x->pinna_mix = x->pinna_mix_target;
  aureonoise_update_pinna_filters(x);
  x->pinna_notchL.clear();
  x->pinna_notchR.clear();
  x->field.temperature = 0.45;

  cefg_clear_all(x);
  aureonoise_spatial_profile_params_changed(x);

#if AUREO_THERMO_LATTICE
  x->ou_pan.tau = 0.60; x->ou_pan.y = 0.0;
  x->ou_itd.tau = 0.40; x->ou_itd.y = 0.0;
  x->ou_amp.tau = 0.80; x->ou_amp.y = 0.0;
  x->ou_rate.tau = 1.20; x->ou_rate.y = 0.0;
  x->lat.init(static_cast<int>(x->p_lat_x), static_cast<int>(x->p_lat_y), static_cast<int>(x->p_lat_z));
  x->lat.eps = x->p_lat_eps;
  x->lat.gamma = x->p_lat_gamma;
  x->lat.sigma = x->p_lat_sigma;
  x->lat_phase = 0.0;
  systhread_mutex_new(&x->lat_mu, 0);
  x->lat_last_v = 0.0;
#if AUREO_BURST_HAWKES
  x->hawkes.base = 4.0;
  x->hawkes.beta = 30.0;
  x->hawkes.lambda = 0.0;
#endif
#endif
  aureonoise_reset_external_buffers(x);
  aureonoise_refresh_harmony(x);
  aureonoise::control::tempo_reset(x);
  aureonoise_update_sync(x, 0.0);
}

void aureonoise_prepare_dsp(t_aureonoise* x, double sr)
{
  if (!x) return;
  x->sr = (sr > 0 ? sr : 44100.0);
  for (auto& g : x->grains) g = {};
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  aureonoise_reset_reports(x);
  x->pinna.width = x->p_width;
  x->pinna.itd_us = x->p_itd_us;
  x->pinna.ild_db = x->p_ild_db;
  const bool was_on = x->pinna_enabled;
  x->pinna_enabled = (x->p_pinna_on != 0);
  if (was_on != x->pinna_enabled) {
    x->pinna_notchL.clear();
    x->pinna_notchR.clear();
  }
  aureonoise_update_pinna_state(x);
  x->pinna_mix = x->pinna_mix_target;
  aureonoise_update_pinna_filters(x);
#if AUREO_THERMO_LATTICE
  aureonoise_lattice_safe_resize(x,
                                 static_cast<int>(x->p_lat_x),
                                 static_cast<int>(x->p_lat_y),
                                 static_cast<int>(x->p_lat_z));
  if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
  x->lat.eps = x->p_lat_eps;
  x->lat.gamma = x->p_lat_gamma;
  x->lat.sigma = x->p_lat_sigma;
  if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
#endif
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
  x->phi_last_mag = 0.0;
  x->phi_last_sign = 1;
  x->last_handshake_time = 0.0;
  x->dialogue = {};
  x->dialogue.coherence = 1.0;
  x->modal.reset();
  x->modal_energy = 0.0;
  aureonoise_update_modal(x);
  aureonoise_reset_external_buffers(x);
  aureonoise_refresh_harmony(x);
  aureonoise::control::tempo_reset(x);
  aureonoise_update_sync(x, 0.0);
}

void aureonoise_clear(t_aureonoise* x)
{
  if (!x) return;
  x->ring.clear();
  for (auto& g : x->grains) g = {};
  x->samples_to_next = static_cast<int>(std::max(1.0, x->sr * 0.05));
  x->gap_elapsed = 0;
  x->sample_counter = 0;
  x->noise.reset();
  aureonoise_refresh_harmony(x);
  x->pinna_mix = x->pinna_mix_target;
  x->pinna_notchL.clear();
  x->pinna_notchR.clear();
  x->prev_pan = 0.0;
  x->prev_itd = 0.0;
  x->prev_ild = 0.0;
  x->last_gap_samples = 0.0;
  x->last_dur_samples = 0.0;
  x->phi_last_mag = 0.0;
  x->phi_last_sign = 1;
  x->last_handshake_time = 0.0;
  x->dialogue = {};
  x->dialogue.coherence = 1.0;
  x->modal.reset();
  x->modal_energy = 0.0;
  aureonoise_update_modal(x);
  aureonoise_reset_reports(x);
#if AUREO_THERMO_LATTICE
  x->lat_phase = 0.0;
  x->lat_last_v = 0.0;
#endif
  aureonoise_reset_external_buffers(x);
  aureonoise::control::tempo_reset(x);
  aureonoise_update_sync(x, 0.0);
}

void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short* count, double sr, long n, long flags)
{
  (void)count;
  (void)n;
  (void)flags;
  aureonoise_prepare_dsp(x, sr);
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-function-type-mismatch"
#endif
  object_method(dsp64, gensym("dsp_add64"), x, reinterpret_cast<method>(aureonoise_perform64), 0, NULL);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
}

void aureonoise_perform64(t_aureonoise* x, t_object*, double** ins, long numins,
                          double** outs, long numouts, long sampleframes, long, void*)
{
  if (!x || x->ob.z_disabled) return;
  (void)ins;
  (void)numins;
  // Main outs
  double* outL = (numouts >= 1 && outs[0]) ? outs[0] : nullptr;
  double* outR = (numouts >= 2 && outs[1]) ? outs[1] : nullptr;
  if (!outL || !outR) return;

  const double block_dt = (x->sr > 0.0) ? (static_cast<double>(sampleframes) / x->sr) : 0.0;
  aureonoise_update_sync(x, block_dt);
  aureonoise_refresh_harmony(x);

#ifdef __SSE2__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  if (x->pinna_filters_dirty) {
    aureonoise_update_pinna_filters(x);
  }
  const bool pinna_now = (x->p_pinna_on != 0);
  if (pinna_now != x->pinna_enabled) {
    x->pinna_enabled = pinna_now;
    x->pinna_notchL.clear();
    x->pinna_notchR.clear();
  }
  aureonoise_update_pinna_state(x);
  double pinna_mix = x->pinna_mix;
  const double pinna_target = x->pinna_mix_target;
  const double pinna_depth_amt = aureo::clamp01(x->p_pinna_depth / 24.0);
  double mix_denom = x->sr * 0.02;
  if (mix_denom < 1.0) mix_denom = 1.0;
  const double pinna_mix_step = 1.0 / mix_denom;
  bool pinna_active = x->pinna_enabled;

  const bool phi_mode = (x->p_phi_mode != 0);
  if (phi_mode) {
    // Disable legacy global pinna filters when using Φ pipeline.
    pinna_mix = 0.0;
    pinna_active = false;
  }

  uint32_t wi = x->ring.writeIndex;
  x->noise.mode = static_cast<aureo::NoiseColor>(static_cast<int>(x->p_noise_mode));
  x->noise.classic.color = static_cast<aureo::NoiseColor>(x->p_noise_color);
  x->noise.classic.amount = aureo::clamp01(x->p_color_amt);
  const bool quantum_mode = (x->noise.mode == aureo::NoiseColor::Quantum);
  if (quantum_mode) {
    const double mix_q = aureo::clamp01(x->p_quantum_mix);
    const double detail_q = aureo::clamp01(x->p_quantum_detail);
    const double stride_q = std::max(0.1, x->p_quantum_stride);
    const double base_q = aureo::clamp(x->p_quantum_base, 20.0, 2000.0);
    const int maxP = aureo::NoiseAureoState::kMaxPartials;
    const int minP = 6;
    int target = static_cast<int>(std::round(minP + detail_q * (maxP - minP)));
    target = std::clamp(target, 4, maxP);
    x->noise.aureo.phi_pi_mix = mix_q;
    x->noise.aureo.planck_decay = aureo::map_phi_range(0.12, 0.52, 1.0 - detail_q);
    x->noise.aureo.prime_stride = stride_q;
    x->noise.aureo.active_partials = target;
    x->noise.aureo.base_freq = base_q;
    x->noise.aureo.velvet_on = (x->p_quantum_velvet > 1.0e-4);
    x->noise.aureo.velvet_amount = aureo::clamp(x->p_quantum_velvet, 0.0, 1.0) * 0.3;
  } else {
    x->noise.aureo.phi_pi_mix = aureo::clamp01(x->p_aureo_mix);
    x->noise.aureo.planck_decay = aureo::clamp01(x->p_aureo_decay);
    x->noise.aureo.prime_stride = std::max(0.1, x->p_aureo_stride);
    x->noise.aureo.active_partials = static_cast<int>(std::clamp(x->p_aureo_harmonics, 1.0, static_cast<double>(aureo::NoiseAureoState::kMaxPartials)));
    x->noise.aureo.base_freq = std::max(10.0, x->sr * 0.25 * aureo::kInvSqrt2);
    x->noise.aureo.velvet_on = (x->p_aureo_velvet != 0);
    x->noise.aureo.velvet_amount = 0.10;
  }

  const double wowHz = aureo::map_phi_range(0.1, 1.5, aureo::clamp01(x->p_vhs_wow));
  const double fltHz = aureo::map_phi_range(7.0, 12.0, aureo::clamp01(x->p_vhs_flutter));
  const double incWow = wowHz / x->sr;
  const double incFlt = fltHz / x->sr;
#if AUREO_THERMO_LATTICE
  const double lat_inc = aureo::clamp(x->p_lat_rate, 1.0, 2000.0) / x->sr;
#endif
  const double itd_scale = x->p_itd_us * 1.0e-6 * x->sr;
  const double min_pan_norm = aureo::clamp(x->p_spat_min_deg, 0.0, 180.0) / 90.0;
  const double min_time_samples = aureo::clamp(x->p_spat_min_ms, 0.0, 500.0) * 0.001 * x->sr;
  const double ipd_amt_global = aureo::clamp01(x->p_spat_ipd);
  const double shadow_amt_global = aureo::clamp01(x->p_spat_shadow);
  const AureonoiseExternalState extState = aureonoise_prepare_external_state(x);
  const double external_amt = extState.amount;
  uint32_t ext_write = x->ext_write_idx;

  phi::Geometry phi_geom{};
  phi::DistanceResponse phi_distance_resp{};
  const double phi_distance = aureo::clamp01(x->p_phi_distance);
  const double phi_elev = aureo::clamp01(x->p_phi_elev);
  const double phi_elev_notch = aureo::clamp01(x->p_phi_elev_notch);
  const double phi_pitch_rad = 12.0 * (aureo::kPi / 180.0);
  if (phi_mode) {
    phi::GeometryConfig cfg;
    cfg.head_b = aureo::clamp(x->p_phi_head_b, 0.06, 0.11);
    cfg.phi_ratio = std::max(1.0, x->p_phi_ratio);
    phi_geom = phi::build_geometry(cfg);
    phi_distance_resp = phi::compute_distance_response(x->sr, phi_distance);
  }

  for (long n = 0; n < sampleframes; ++n) {
    ++x->gap_elapsed;
    x->lfo_wow_phase += incWow; if (x->lfo_wow_phase >= 1.0) x->lfo_wow_phase -= 1.0;
    x->lfo_flut_phase += incFlt; if (x->lfo_flut_phase >= 1.0) x->lfo_flut_phase -= 1.0;
    const double wow = std::sin(2.0 * aureo::kPi * x->lfo_wow_phase);
    const double flt = std::sin(2.0 * aureo::kPi * x->lfo_flut_phase);
    const double vhs_mod = 0.5 * wow + 0.5 * flt;

    pinna_mix += (pinna_target - pinna_mix) * pinna_mix_step;
    pinna_mix = aureo::clamp01(pinna_mix);

#if AUREO_THERMO_LATTICE
    x->lat_phase += lat_inc;
    if (x->lat_phase >= 1.0) {
      const int k = static_cast<int>(std::floor(x->lat_phase));
      const double dt = static_cast<double>(k) / std::max(1.0, x->p_lat_rate);
      bool locked = false;
      if (x->lat_mu) locked = (systhread_mutex_trylock(x->lat_mu) == 0);
      if (!x->lat_mu || locked) {
        for (int i = 0; i < k; ++i) x->lat.step(x->rng);
      }
      if (locked && x->lat_mu) systhread_mutex_unlock(x->lat_mu);
      const double T = aureo::clamp01(x->p_T);
      x->ou_pan.sigma = 0.40 * T;
      x->ou_itd.sigma = 0.35 * T;
      x->ou_amp.sigma = 0.30 * T;
      x->ou_rate.sigma = 0.25 * T;
      x->ou_pan.step(dt, 0.0, x->rng);
      x->ou_itd.step(dt, 0.0, x->rng);
      x->ou_amp.step(dt, 0.0, x->rng);
      x->ou_rate.step(dt, 0.0, x->rng);
#if AUREO_BURST_HAWKES
      (void)x->hawkes.tick(dt, x->rng);
#endif
      x->lat_phase -= static_cast<double>(k);
    }
#endif

    double nz = x->noise.process(x->rng, x->sr);
#ifdef __aarch64__
    nz += 1.0e-18 * x->rng.uniPM1();
#endif
    // Mix in physical entropy: stronger for Quantum/Velvet modes, light for Aureo
    if (x->noise.mode == aureo::NoiseColor::Quantum || x->noise.mode == aureo::NoiseColor::Velvet) {
      nz = 0.92 * nz + 0.08 * g_entropy.unit_pm1();
    } else if (x->noise.mode == aureo::NoiseColor::Aureo) {
      nz = 0.97 * nz + 0.03 * g_entropy.unit_pm1();
    }
    nz = aureo::soft_tanh(nz * 1.2);
    x->ring.data[wi] = nz;

    if (--x->samples_to_next <= 0) {
      const int gi = aureonoise::scheduler::find_free_grain(x);
      bool spawned = false;
      if (gi >= 0) {
        auto& g = x->grains[gi];
        g = {};
        auto& cefg_state = x->cefg_grains[gi];
        cefg_clear_grain_state(cefg_state);

        const double u1 = x->w_phi.next();
        const double u2 = x->w_s2.next();
        const double u3 = x->w_pl.next();
        const double u4 = x->rng.uni01();
        const double u5 = x->rng.uni01();
        const double u6 = x->rng.uni01();

        const double gap_samples = static_cast<double>(x->gap_elapsed);
        const double prev_dur = std::max(1.0, x->last_dur_samples);
        const double ratio_gap = gap_samples / prev_dur;
        const double coupling = aureo::clamp01(x->p_hemis_coupling);
        const double time_weight = aureo::clamp(0.35 + 0.65 * (ratio_gap / (ratio_gap + 1.0)), 0.35, 1.0);
        const double hemi = coupling * time_weight;

#if AUREO_THERMO_LATTICE
        double lat_u = 0.5;
        if (x->p_lattice) {
          double v = x->lat_last_v;
          bool locked_lat = false;
          if (x->lat_mu) locked_lat = (systhread_mutex_trylock(x->lat_mu) == 0);
          if (!x->lat_mu || locked_lat) {
            v = x->lat.probe(x->w_phi.next());
            x->lat_last_v = v;
          }
          if (locked_lat && x->lat_mu) systhread_mutex_unlock(x->lat_mu);
          lat_u = 0.5 + 0.5 * std::tanh(v);
        }
        const double oup = x->p_thermo ? aureo::clamp(x->ou_pan.y, -1.0, 1.0) : 0.0;
        const double oua = x->p_thermo ? aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, 0.5 + 0.5 * std::tanh(x->ou_amp.y)) : 1.0;
        const double oui = x->p_thermo ? x->ou_itd.y : 0.0;
#else
        const double lat_u = 0.5;
        const double oup = 0.0;
        const double oua = 1.0;
        const double oui = 0.0;
#endif

        const double amp_shape = std::pow(std::max(1e-9, u2), 0.35);
        const double amp_lat = aureo::map_phi_range(1.0 / aureo::kPhi, aureo::kPhi, lat_u);
        const double amp = aureo::kAmpNorm * amp_shape * amp_lat * oua;
        const int base_dur_samples = aureonoise::scheduler::map_len_samples(x, u1);
        const double base_dur = std::max(1.0, static_cast<double>(base_dur_samples));

        double pan = 2.0 * u3 - 1.0;
#if AUREO_THERMO_LATTICE
        pan += 0.25 * (2.0 * lat_u - 1.0) + 0.35 * oup;
#endif
        pan = aureo::clamp((1.0 - hemi) * pan - hemi * x->prev_pan, -1.0, 1.0);
        pan = aureonoise_apply_phi_pan(x, pan);

        const double timestamp_sec = (x->sr > 0.0)
                                       ? (static_cast<double>(x->sample_counter) / x->sr)
                                       : 0.0;
        auto dialogue_result = aureonoise::dialogue::evaluate(
            x, timestamp_sec, pan, amp, base_dur, gap_samples);
        pan = dialogue_result.pan;
        double amp_dialogue = amp * dialogue_result.amp_scale;
        double dur_dialogue = std::max(1.0, base_dur * dialogue_result.dur_scale);

        bool ok = aureonoise_poisson_enforce(x, pan, min_pan_norm, min_time_samples);
        aureonoise_commit_phi_state(x, pan, ok);
        if (!ok) {
          aureonoise::dialogue::commit(x, dialogue_result, false);
        } else {
          const uint64_t report_idx = x->report_total.fetch_add(1, std::memory_order_relaxed) + 1;
          const int dur_samples = std::max<int>(1, static_cast<int>(std::round(dur_dialogue)));
          g.dur = dur_samples;
          const double sr_ref = (x->sr > 0.0) ? x->sr : 44100.0;
          const double dur_norm = aureo::clamp(static_cast<double>(g.dur) / (sr_ref * 0.35), 0.0, 1.0);
          const double gap_norm = aureo::clamp(gap_samples / (sr_ref * 0.5), 0.0, 1.0);
          double burst_amp_scale = 1.0;
          pan = aureonoise_apply_burst_position(x, pan, dur_norm, gap_norm, burst_amp_scale);
          double amp_final = amp_dialogue * burst_amp_scale;
          if (x->p_spat_mirror) {
            pan = -pan;
          }
          g.amp = amp_final;
          g.pan = pan;

          double ild_db = 0.0;
          aureo::BinauralCoefficients binaural{};
          bool have_binaural = false;
          if (phi_mode) {
            double panGainL = 1.0;
            double panGainR = 1.0;
            aureo::pan_equal_power(pan, 1.0, panGainL, panGainR);
            g.panL = panGainL;
            g.panR = panGainR;
            const double ext_gain = 1.0 + 0.18 * external_amt * std::abs(pan);
            g.panL *= ext_gain;
            g.panR *= ext_gain;
            g.crossfeed = 0.0;
            g.focus = 0.0;

            const double theta = aureo::clamp(pan, -1.0, 1.0) * (aureo::kPi * 0.5);
            std::array<double, 3> dir{std::sin(theta), 0.0, std::cos(theta)};
            std::array<double, 3> dir_pitch = rotate_x(dir, -phi_pitch_rad);
            const double ux_pitch = dir_pitch[0];
            const double uy_pitch = dir_pitch[1];
            const double uz_pitch = dir_pitch[2];
            cefg_prepare_grain(x, static_cast<size_t>(gi), phi_geom, x->sr,
                               phi_distance, phi_pitch_rad,
                               ux_pitch, uy_pitch, uz_pitch);
            const auto headRes = phi::compute_head_result(phi_geom, x->sr, pan, phi_pitch_rad, phi_distance);
            double itd = headRes.itd_samples;
#if AUREO_THERMO_LATTICE
            itd += ((2.0 * lat_u - 1.0) * 0.33 + 0.33 * oui) * (x->p_itd_us * 1.0e-6 * x->sr);
#endif
            g.itd = aureo::clamp((1.0 - hemi) * itd - hemi * x->prev_itd, -itd_scale, itd_scale);

            const double ild_max = aureo::clamp(x->p_ild_db, 0.0, 24.0);
            const double ild_shape = std::pow(std::abs(headRes.lateral), 1.3);
            double ild = ild_max * ild_shape * (0.70 + 0.30 * (1.0 - phi_distance));
            const bool toRight = (ux_pitch >= 0.0);
            ild_db = toRight ? ild : -ild;
            const double gIL = aureo::db_to_lin(toRight ? -ild : +ild);
            const double gIR = aureo::db_to_lin(toRight ? +ild : -ild);
            g.gL = gIL;
            g.gR = gIR;

            const double ext_shadow = std::abs(headRes.lateral);
            const double reff_dir = headRes.direction_radius;
            const double fc_open = std::min(9000.0, 0.35 * 343.0 / std::max(1.0e-4, reff_dir));
            double fc_occ = std::max(1200.0, fc_open * (0.55 + 0.45 * (1.0 - ext_shadow)));
            const double alpha_occ = std::exp(-2.0 * aureo::kPi * fc_occ / x->sr);
            g.shadow_a = alpha_occ;
            g.shadow_left = (pan > 0.0);
            g.shadow_right = (pan < 0.0);

            double amp_scale_final = 1.0;
            amp_scale_final *= phi_distance_resp.direct_gain;
            g.amp *= amp_scale_final;
          } else {
            binaural = x->pinna.coefficients(x->sr, pan, u2, u5);
            have_binaural = true;
            g.panL = binaural.gainL;
            g.panR = binaural.gainR;
            const double ext_gain = 1.0 + 0.18 * external_amt * std::abs(pan);
            g.panL *= ext_gain;
            g.panR *= ext_gain;
            g.crossfeed = binaural.crossfeed;
            g.focus = aureo::clamp01(binaural.focus + 0.25 * external_amt);

            double itd = binaural.itd_samples;
#if AUREO_THERMO_LATTICE
            itd += ((2.0 * lat_u - 1.0) * 0.33 + 0.33 * oui) * (x->p_itd_us * 1.0e-6 * x->sr);
#endif
            g.itd = (1.0 - hemi) * itd - hemi * x->prev_itd;

            const double ild_max = aureo::clamp(x->p_ild_db, 0.0, 24.0);
            ild_db = binaural.ild_db;
            ild_db = aureo::clamp((1.0 - hemi) * ild_db - hemi * x->prev_ild, -ild_max, ild_max);
            ild_db = aureo::clamp(ild_db * (1.0 + 0.45 * external_amt), -ild_max, ild_max);
            g.crossfeed *= (1.0 - 0.45 * external_amt);
            g.gL = aureo::db_to_lin(ild_db);
            g.gR = aureo::db_to_lin(-ild_db);
          }

          double ipd_coeff = 0.0;
          if (ipd_amt_global > 1.0e-6) {
            const double ipd_shape = std::abs(2.0 * u6 - 1.0);
            const double ipd_base = 0.18 + 0.55 * ipd_amt_global;
            const double ipd_spread = 0.25 * ipd_amt_global;
            ipd_coeff = aureo::clamp(ipd_base + ipd_spread * (ipd_shape - 0.5), 0.0, 0.95);
          }
          ipd_coeff *= (0.7 + 0.3 * g.focus);
          g.ipd_coeff = ipd_coeff;

          const double shadow_factor = shadow_amt_global * std::abs(pan);
          if (shadow_factor > 1.0e-6) {
            const double fc_min = 800.0;
            const double fc_max = std::max(fc_min, std::min(8000.0, 0.45 * x->sr));
            const double fc = aureo::map_phi_range(fc_min, fc_max, 1.0 - shadow_factor);
            const double alpha = std::exp(-2.0 * aureo::kPi * fc / x->sr);
            g.shadow_a = aureo::clamp(alpha, 0.0, 0.9999);
            g.shadow_left = (pan > 0.0);
            g.shadow_right = (pan < 0.0);
          }

          g.kind = aureo::choose_kind(x->p_glitch_mix, u4);

          int srN = aureo::map_sr_hold_base(x->p_srcrush_amt, u1);
#if AUREO_SR_PRIME_SNAP
          g.sr_holdN = pick_prime_in_range(x, std::max(1, srN - 7), srN + 7, u2);
#else
          g.sr_holdN = srN;
#endif
          g.sr_holdCnt = g.sr_holdN;
          const int bits = aureo::map_bits(x->p_bitcrush_amt);
          g.q_levels = (1 << (bits - 1)) - 1;

          const auto envShape = x->field.make_envelope(x->p_env_attack,
                                                       x->p_env_decay,
                                                       x->p_env_sustain,
                                                       x->p_env_release,
                                                       gap_samples,
                                                       static_cast<double>(g.dur),
                                                       std::abs(pan));
          g.env = envShape;
          g.on = true;
          g.age = 0;
          x->prev_pan = pan;
          x->prev_itd = g.itd;
          x->prev_ild = ild_db;
          x->last_gap_samples = gap_samples;
          x->last_dur_samples = static_cast<double>(g.dur);
          x->gap_elapsed = 0;
          t_aureonoise::GrainReport report;
          report.index = report_idx;
          report.timestamp_sec = timestamp_sec;
          report.gap_samples = gap_samples;
          report.dur_samples = static_cast<double>(g.dur);
          report.amp = g.amp;
          report.pan = g.pan;
          report.itd_samples = g.itd;
          report.ild_db = ild_db;
          if (phi_mode) {
            report.binaural_azimuth_deg = pan * 90.0;
            report.binaural_focus = 0.0;
            report.binaural_crossfeed = 0.0;
          } else if (have_binaural) {
            report.binaural_azimuth_deg = binaural.azimuth_rad * (180.0 / aureo::kPi);
            report.binaural_focus = binaural.focus;
            report.binaural_crossfeed = binaural.crossfeed;
          } else {
            report.binaural_azimuth_deg = pan * 90.0;
            report.binaural_focus = 0.0;
            report.binaural_crossfeed = 0.0;
          }
          report.phi_u = u1;
          report.s2_u = u2;
          report.pl_u = u3;
          report.rng4 = u4;
          report.rng5 = u5;
          report.rng6 = u6;
#if AUREO_THERMO_LATTICE
          report.lattice_u = lat_u;
          report.lattice_v = x->lat_last_v;
          report.thermo_on = (x->p_thermo != 0);
          report.lattice_on = (x->p_lattice != 0);
#if AUREO_BURST_HAWKES
          report.hawkes_lambda = x->hawkes.lambda;
          report.burst_on = (x->p_burst != 0) && (x->hawkes.lambda > x->hawkes.base + 1.0);
#else
          report.hawkes_lambda = 0.0;
          report.burst_on = false;
#endif
#else
          report.lattice_u = 0.5;
          report.lattice_v = 0.0;
          report.thermo_on = false;
          report.lattice_on = false;
          report.hawkes_lambda = 0.0;
          report.burst_on = false;
#endif
          auto commit_result = dialogue_result;
          commit_result.pan = pan;
          const double base_amp_safe = std::max(1.0e-9, commit_result.base_amp);
          commit_result.amp_scale = amp_final / base_amp_safe;
          const double base_dur_safe = std::max(1.0, commit_result.base_dur);
          commit_result.dur_scale = static_cast<double>(g.dur) / base_dur_safe;
          aureonoise::dialogue::commit(x, commit_result, true);
          report.dialogue_coherence = x->dialogue.coherence;
          report.dialogue_reply = commit_result.handshake ? 1.0 : 0.0;
          aureonoise_log_grain_event(x, report);
          spawned = true;
        }
      }
      x->samples_to_next = aureonoise::scheduler::schedule_gap_samples(x);
      if (!spawned && gi >= 0) {
        x->grains[gi].on = false;
      }
    }

    double yL = 0.0, yR = 0.0;
    for (size_t gi = 0; gi < x->grains.size(); ++gi) {
      auto& g = x->grains[gi];
      auto& cefgState = x->cefg_grains[gi];
      if (!g.on) { cefgState.tap_count = 0; continue; }
      if (g.age >= g.dur) { g.on = false; cefgState.tap_count = 0; continue; }

      const double phase = static_cast<double>(g.age) / static_cast<double>(g.dur);
      const double env = x->field.envelope(phase, g.env);

      double itd = g.itd + (vhs_mod * 0.25 * itd_scale);
#if AUREO_THERMO_LATTICE
      if (x->p_lattice) itd += 0.18 * std::tanh(x->lat_last_v) * itd_scale;
#endif
      double sL, sR;
      aureo::ring_read_stereo_itd_frac(x->ring, wi, itd, sL, sR);
      sL *= g.panL * g.gL;
      sR *= g.panR * g.gR;

      if (std::abs(g.crossfeed) > 1.0e-6) {
        const double baseL = sL;
        const double baseR = sR;
        sL = baseL + g.crossfeed * baseR;
        sR = baseR + g.crossfeed * baseL;
      }

      if (--g.sr_holdCnt <= 0) {
        g.heldL = sL;
        g.heldR = sR;
        g.sr_holdCnt = g.sr_holdN;
      }
      sL = g.heldL;
      sR = g.heldR;

      if (g.q_levels > 0) {
        sL = std::round(sL * g.q_levels) / static_cast<double>(g.q_levels);
        sR = std::round(sR * g.q_levels) / static_cast<double>(g.q_levels);
      }

      if (g.kind == aureo::GrainKind::VhsDrop) {
        const double att = 0.5 + 0.5 * (1.0 - std::abs(vhs_mod));
        sL *= att; sR *= att;
      } else if (g.kind == aureo::GrainKind::Stutter) {
        if ((g.age & 7) == 0) { sL *= 0.2; sR *= 0.2; }
      }

      if (g.ipd_coeff > 1.0e-6) {
        sL = aureonoise_allpass(sL, g.ipd_coeff, g.ipd_zL);
        sR = aureonoise_allpass(sR, -g.ipd_coeff, g.ipd_zR);
      }

      if (g.shadow_a > 1.0e-6) {
        if (g.shadow_left) {
          sL = aureonoise_shadow_lp(sL, g.shadow_a, g.shadow_zL);
        }
        if (g.shadow_right) {
          sR = aureonoise_shadow_lp(sR, g.shadow_a, g.shadow_zR);
        }
      }

      if (phi_mode) {
        if (cefgState.tap_count > 0) {
          double erL = 0.0;
          double erR = 0.0;
          for (int ti = 0; ti < cefgState.tap_count; ++ti) {
            auto& tap = cefgState.taps[ti];
            const double sampleL = ring_read_delay(x->ring, wi, tap.base_delay + tap.itd_offset_l);
            const double sampleR = ring_read_delay(x->ring, wi, tap.base_delay + tap.itd_offset_r);
            if (tap.lp_alpha > 0.0) {
              tap.lp_state_l = tap.lp_alpha * tap.lp_state_l + (1.0 - tap.lp_alpha) * sampleL;
              tap.lp_state_r = tap.lp_alpha * tap.lp_state_r + (1.0 - tap.lp_alpha) * sampleR;
            } else {
              tap.lp_state_l = sampleL;
              tap.lp_state_r = sampleR;
            }
            erL += tap.mix * tap.gain_l * tap.lp_state_l;
            erR += tap.mix * tap.gain_r * tap.lp_state_r;
          }
          sL += erL;
          sR += erR;
        }
      } else {
        cefgState.tap_count = 0;
      }

      const double amp_env = g.amp * env;
      yL += amp_env * sL;
      yR += amp_env * sR;
      ++g.age;
    }

    const double modalMix = aureo::clamp(x->p_modal_mix, 0.0, 1.0);
    if (x->p_modal_on && modalMix > 1.0e-6) {
      double modalOutL = yL;
      double modalOutR = yR;
      x->modal.process(yL, yR, modalOutL, modalOutR);
      if (x->p_modal_mirror) {
        modalOutL = -modalOutL;
        std::swap(modalOutL, modalOutR);
      }
      yL = (1.0 - modalMix) * yL + modalMix * modalOutL;
      yR = (1.0 - modalMix) * yR + modalMix * modalOutR;
      const double modalEnergySample = 0.5 * (modalOutL * modalOutL + modalOutR * modalOutR);
      x->modal_energy = 0.995 * x->modal_energy + 0.005 * modalEnergySample;
    } else {
      x->modal_energy *= 0.995;
    }

    aureonoise_apply_external_mix(x, extState, yL, yR, ext_write);

    if (pinna_active) {
      const double dryL = yL;
      const double dryR = yR;
      const double depthAmt = pinna_depth_amt;
      double wetL = x->pinna_notchL.proc(dryL);
      double wetR = x->pinna_notchR.proc(dryR);
      wetL = dryL + depthAmt * (wetL - dryL);
      wetR = dryR + depthAmt * (wetR - dryR);
      const double wet = pinna_mix;
      const double dry_w = 1.0 - wet;
      yL = dry_w * dryL + wet * wetL;
      yR = dry_w * dryR + wet * wetR;
    }

    yL = std::tanh(yL * aureo::kOutDrive) / aureo::kOutDrive;
    yR = std::tanh(yR * aureo::kOutDrive) / aureo::kOutDrive;

    outL[n] = yL;
    outR[n] = yR;

    wi = (wi + 1u) & aureo::kRingMask;
    ++x->sample_counter;
  }

  x->pinna_mix = pinna_mix;
  x->ring.writeIndex = wi;
  if (extState.active) x->ext_write_idx = ext_write;
}

#if AUREO_THERMO_LATTICE
void aureonoise_lattice_safe_resize(t_aureonoise* x, int X, int Y, int Z)
{
  if (!x) return;
  X = std::max(2, X);
  Y = std::max(2, Y);
  Z = std::max(1, Z);
  const size_t n = static_cast<size_t>(X) * static_cast<size_t>(Y) * static_cast<size_t>(Z);
  std::vector<double> newx(n, 0.0);
  std::vector<double> newtmp(n, 0.0);
  if (x->lat_mu) systhread_mutex_lock(x->lat_mu);
  x->lat.X = X;
  x->lat.Y = Y;
  x->lat.Z = Z;
  x->lat.x.swap(newx);
  x->lat.tmp.swap(newtmp);
  x->lat_phase = 0.0;
  x->lat_last_v = 0.0;
  if (x->lat_mu) systhread_mutex_unlock(x->lat_mu);
}
#endif
