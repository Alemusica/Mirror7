#pragma once

#include <juce_dsp/juce_dsp.h>
#include <memory>
#include <vector>

#include "aureo_noise.hpp"
#include "aureo_binaural.hpp"
#include "aureo_stoch.hpp"
#include "aureo_math.hpp"
#include "aureo_ring.hpp"
#include "aureo_field_adsr.hpp"
#include "aureo_weyl.hpp"
#include "aureo_pinna.hpp"
#include "aureonoise/spatial/phi_model.hpp"

#if __has_include("aureonoise3beta/legacy/core/modal_engine.hpp")
#include "aureonoise3beta/legacy/core/modal_engine.hpp"
#else
#include "core/modal_engine.hpp"
#endif

// Optional: early reflections profile (CEFG). We only declare type; can be disabled.
#include "aureonoise/spatial/cefg.hpp"

class Mirror7Engine
{
public:
    struct Params
    {
        // timing/scheduling
        double rateHz = 8.0;
        double baseLenMs = 120.0;
        double lenPhi = 0.8;
        double width = 1.0;

        // binaural
        double itd_us = 700.0;
        double ild_db = 20.0;
        double spatMinDeg = 12.0;
        double spatMinMs = 35.0;
        double spatIPD = 0.6;
        double spatShadow = 0.7;
        double spatExternal = 0.0; // not used yet
        bool   spatMirror = false;

        // envelope
        double envAttack = 0.18;
        double envDecay = 0.28;
        double envSustain = 0.55;
        double envRelease = 0.30;

        // coupling
        double hemisCoupling = 0.6;

        // noise
        int    noiseColor = (int) aureo::NoiseColor::Pink;
        double colorAmt = 0.65;
        int    noiseMode = (int) aureo::NoiseColor::Pink; // Pink/Aureo/Quantum/Velvet etc.

        // aureo
        double aureoMix = 0.5;
        double aureoDecay = 0.3;
        double aureoStride = 1.0;
        double aureoHarmonics = 12.0;
        bool   aureoVelvet = true;

        // quantum
        double qMix = 0.55;
        double qDetail = 0.7;
        double qStride = 1.0;
        double qBase = 220.0;
        double qVelvet = 0.12;

        // vhs
        double vhsWow = 0.35;
        double vhsFlutter = 0.25;

        // glitch
        double glitchMix = 0.5;
        double srCrush = 0.2;
        double bitCrush = 0.15;

        // burst
        double burstFloor = 0.35;
        double burstPhiMix = 0.6;

        // dialogue
        bool   dialogueOn = true;
        double dialogueStrength = 0.6;
        double dialogueMemory = 0.5;
        double dialoguePhiMix = 0.75;

        // modal
        bool   modalOn = false;
        bool   modalMirror = false;
        double modalMix = 0.35;
        double modalDecay = 0.5;
        int    modalPreset = 1;

        // pinna/global
        bool   pinnaOn = false;
        double pinnaDepth = 12.0; // dB range 0..24

        // phi
        bool   phiPan = false;
        bool   phiMode = true;
        double phiRatio = aureo::kPhi;
        double phiHeadB = 0.0875;
        double phiDistance = 0.12;
        double phiElev = 0.6;
        double phiElevNotch = 0.6;
        double phiTorsoMix = 0.12;
        double phiTorsoMs = 1.15;
        double phiTorsoHpHz = 500.0;

        // spatial profile
        int    spatialProfile = 0; // 0 legacy, 1 uteroid_phi
        double cefgGain = 1.1;
        double cefgMix = 3.4;
        bool   cefgDebug = false;

        // host sync
        bool   syncEnable = false;
        double syncDivision = 1.0; // 1=quarter
        double syncSlew = 0.5;     // seconds
    };

    struct DialogueState
    {
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
    };

    Mirror7Engine() = default;

    void prepare (double fs, int maxBlock);
    void reset();
    void setOversampling (int factor); // param choice index: 0->1x, 1->2x, 2->4x
    int  getOversamplingFactor() const { return osFactor; }

    // realtime param updates
    void setVelvetAmount (double amt);
    void setPanWidth (double pan, double width);
    void setITDILD (double itdUs, double ildDb);
    void setAutopan (bool enable, double speed);
    Params& getParams() { return params; }

    // host transport info
    void updateHostTransport (double bpm /* <0 if unknown */, bool isPlaying);

    void process (juce::AudioBuffer<float>& buffer, float outGainDb);

private:
    // helpers
    struct CefgTap { double base_delay=0, itd_l=0, itd_r=0, base_gain_l=0, base_gain_r=0, gain_l=0, gain_r=0, lp_a=0, lp_l=0, lp_r=0, mix=0; };
    struct CefgGrainState { int tap_count=0; std::array<CefgTap, 16> taps{}; };

    double fs = 44100.0;
    int osFactor = 1;
    std::unique_ptr<juce::dsp::Oversampling<float>> oversampling;

    Params params;

    // core state
    aureo::RNG rng;
    aureo::AureoNoiseEngine noise;
    aureo::RingState ring;
    aureo::FieldState field;
    aureo::Weyl w_phi, w_s2, w_pl;
    std::array<aureo::Grain, aureo::kMaxGrains> grains{};
    std::array<CefgGrainState, aureo::kMaxGrains> cefgGrains{};
    uint32_t samples_to_next = 1;
    uint32_t gap_elapsed = 0;
    uint64_t sample_counter = 0;

    // autopan
    aureo::OU panOU, itdOU, ampOU, rateOU;
    bool autopan = false;
    double autopanSpeed = 0.12; // Hz-ish

    // pinna
    double pinna_mix = 0.0, pinna_target = 0.0;
    double pinna_depth = 12.0;

    // last φ state
    double prev_pan = 0.0, prev_itd = 0.0, prev_ild = 0.0;
    double phi_last_mag = 0.0; int phi_last_sign = 1;

    // wow/flutter
    double lfo_wow = 0.0, lfo_flt = 0.0;

    // dialogue
    DialogueState dialogue;

    // Φ geometry and distance
    aureonoise::spatial::phi::Geometry phiGeom{};
    bool phiGeomValid = false;

    // CEFG profile
    aureonoise::spatial::cefg::Profile cefgProfile;
    bool cefgLoaded = false;

    // fractional delays for ITD
    juce::dsp::DelayLine<float, juce::dsp::DelayLineInterpolationTypes::Lagrange3rd> dL { 4096 }, dR { 4096 };

    // smoothing for static pan param
    juce::SmoothedValue<float, juce::ValueSmoothingTypes::Linear> panSmoothed;

    // util
    static inline double allpass (double in, double a, double& z) { const double y = -a * in + z; z = in + a * y; return y; }
    static inline double shlp    (double in, double a, double& z) { const double y = (1.0 - a) * in + a * z; z = y; return y; }
    static inline double ringReadDelay (const aureo::RingState& ring, uint32_t wi, double delaySamples)
    {
        const double d = std::max (0.0, delaySamples);
        const double pos = (double) wi - d;
        return aureo::lagrange3 (ring.data.data(), pos);
    }

    // scheduling
    int  findFreeGrain();
    int  scheduleGapSamples();
    uint32_t mapLenSamples (double u);

    // dialogue helpers
    void dialogueReset();
    struct DialogueResult { double pan=0, amp_scale=1, dur_scale=1, next_coherence=1, base_amp=1, base_dur=0, gap_samples=0, timestamp=0, handshake_score=0, handshake_ratio_target=aureo::kPhi; bool handshake=false; int sign=1; int fib_index=1; int fib_direction=1; };
    DialogueResult dialogueEvaluate (double timestamp_sec, double pan, double amp, double dur_samples, double gap_samples);
    void dialogueCommit (const DialogueResult& res, bool accepted);

    // φ helpers
    double applyPhiPan (double pan);
    void   commitPhiState (double pan, bool ok);
    double applyBurstPosition (double pan, double dur_norm, double gap_norm, double& amp_scale);

    // modal
    aureonoise::modal::Engine modal;
    double modalEnergy = 0.0;

    // pinna global (legacy-style notch)
    aureo::Biquad pinnaL, pinnaR;
    bool pinnaFiltersDirty = true;
    bool pinnaEnabled = false;
    double pinnaMinHz = 400.0;
    double pinnaBaseHz = 7800.0;
    double pinnaSpread = 2200.0;
    double pinnaQ = 5.0;
    double pinnaFreqL = 7200.0;
    double pinnaFreqR = 8200.0;

    // tempo state
    struct Tempo { double fallback_rate=4.0, manual_rate=4.0, sync_rate=0.0, current_rate=4.0, bpm_target=120.0, bpm_current=120.0; bool transport_running=false, pending_start=false; } tempo;
    double hostBpm = -1.0; bool hostPlaying = false; bool lastHostPlaying = false;

    void updateTempo (double block_dt);
    void updatePinnaFilters();
    void maybeLoadCefg();
    void cefgPrepareGrain (size_t grainIndex, const aureonoise::spatial::phi::Geometry& geom, double sampleRate, double phi_distance, double pitch_rad);
};
