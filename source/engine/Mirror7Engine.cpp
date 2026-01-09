#include "engine/Mirror7Engine.h"
#include <cmath>

using aureo::clamp; using aureo::clamp01; using aureo::map_phi_range; using aureo::db_to_lin;

Mirror7Engine::Mirror7Engine()
    : spatializer (AUREONOISE_SPATIAL_ASSETS_DIR)
{
}

void Mirror7Engine::prepare (double sampleRate, int maxBlock)
{
    fs = (sampleRate > 0.0 ? sampleRate : 44100.0);
    const juce::uint32 blockSize = (juce::uint32) juce::jmax (1, maxBlock);
    juce::dsp::ProcessSpec spec { fs, blockSize, 1 };
    dL.prepare (spec);
    dR.prepare (spec);
    dL.reset();
    dR.reset();

    noiseController.reset();
    ring.data.fill (0.0); ring.writeIndex = 0;
    // FieldState holds only parameters and last stats; no reset needed
    w_phi.x = rng.uni01(); w_phi.set_step (aureo::kInvPhi);
    w_s2.x  = rng.uni01(); w_s2.set_step (aureo::kInvSqrt2);
    w_pl.x  = rng.uni01(); w_pl.set_step (aureo::kInvPlastic);
    panOU.tau = 2.0; panOU.sigma = 0.15; itdOU.tau = 2.0; itdOU.sigma = 0.12; ampOU.tau = 2.0; ampOU.sigma = 0.12; rateOU.tau = 2.0; rateOU.sigma = 0.10;
    lfo_wow = lfo_flt = 0.0;
    pinna_mix = 0.0; pinna_target = params.pinnaOn ? 1.0 : 0.0; pinna_depth = clamp (params.pinnaDepth, 0.0, 24.0);
    for (auto& g : grains) g.on = false;
    samples_to_next = 1; gap_elapsed = 0; sample_counter = 0;
    panSmoothed.reset (fs, 0.02); panSmoothed.setCurrentAndTargetValue (0.0f);

    // Configure noise defaults
    noiseController.configure (params, fs);

    // Φ geometry
    SpatialParamsView spatialParams { params.phiMode, params.phiHeadB, params.phiRatio, params.cefgGain, params.cefgMix, params.spatialProfile };
    spatializer.reset (fs);
    spatializer.updateGeometry (spatialParams);

    // Dialogue
    dialogue.setSampleRate (fs);
    dialogue.reset();

    // Modal engine
    modal.setSampleRate (fs);
    modal.configure (params.modalPreset, params.modalDecay);

    // Tempo init
    tempo.fallback_rate = (params.rateHz > 1.0e-6 ? params.rateHz : 4.0);
    tempo.manual_rate = tempo.fallback_rate;
    tempo.current_rate = tempo.fallback_rate;
    tempo.sync_rate = 0.0;
    tempo.bpm_target = (hostBpm > 0.0 ? hostBpm : 120.0);
    tempo.bpm_current = tempo.bpm_target;
    tempo.transport_running = false;
    tempo.pending_start = false;
    lastHostPlaying = hostPlaying;

    // Pinna filters initial setup
    pinnaFiltersDirty = true;
}

void Mirror7Engine::reset()
{
    for (auto& g : grains) g.on = false;
    ring.data.fill (0.0); ring.writeIndex = 0; samples_to_next = 1; gap_elapsed = 0; sample_counter = 0;
    dL.reset();
    dR.reset();
}

void Mirror7Engine::setOversampling (int factor)
{
    factor = juce::jlimit (0, 2, factor); // 0,1,2
    const int newFactor = (factor == 0 ? 1 : (factor == 1 ? 2 : 4));
    if (newFactor == osFactor) return;
    osFactor = newFactor;
    if (osFactor == 1) { oversampling.reset(); return; }
    const int stages = (osFactor == 2 ? 1 : 2);
    oversampling = std::make_unique<juce::dsp::Oversampling<float>> (2, stages, juce::dsp::Oversampling<float>::filterHalfBandPolyphaseIIR, true);
    oversampling->initProcessing (512);
}

void Mirror7Engine::updateHostTransport (double bpm, bool isPlaying)
{
    hostBpm = bpm;
    hostPlaying = isPlaying;
}

void Mirror7Engine::updateTempo (double block_dt)
{
    const bool syncEnabled = params.syncEnable && (params.rateHz <= 1.0e-6);
    // manual branch
    tempo.fallback_rate = (params.rateHz > 1.0e-6 ? clamp (params.rateHz, 0.01, aureo::kMaxEventRateHz) : tempo.fallback_rate);
    tempo.manual_rate = (params.rateHz > 1.0e-6 ? tempo.fallback_rate : tempo.fallback_rate);

    if (!syncEnabled)
    {
        tempo.sync_rate = 0.0;
        tempo.current_rate = tempo.manual_rate;
        tempo.bpm_target = (hostBpm > 0.0 ? hostBpm : tempo.bpm_target);
        tempo.bpm_current = tempo.bpm_target;
    }
    else
    {
        // sync branch
        if (hostBpm > 0.0)
            tempo.bpm_target = hostBpm;
        const double slew = juce::jlimit (0.0, 5.0, params.syncSlew);
        if (block_dt > 0.0 && slew > 0.0)
        {
            const double tau = std::max (0.01, slew);
            const double alpha = std::exp (-block_dt / tau);
            tempo.bpm_current = tempo.bpm_target + (tempo.bpm_current - tempo.bpm_target) * alpha;
        }
        else
            tempo.bpm_current = tempo.bpm_target;
        const double div = juce::jlimit (0.0625, 32.0, params.syncDivision);
        tempo.sync_rate = juce::jlimit (0.01, aureo::kMaxEventRateHz, (tempo.bpm_current / 60.0) * div);
        tempo.current_rate = tempo.sync_rate;
    }

    // transport edge
    if (syncEnabled)
    {
        if (!lastHostPlaying && hostPlaying)
        {
            tempo.transport_running = true;
            tempo.pending_start = true;
        }
        else if (!hostPlaying)
        {
            tempo.transport_running = false;
        }
        lastHostPlaying = hostPlaying;
    }
}

void Mirror7Engine::updatePinnaFilters()
{
    double width = clamp (params.width, 0.0, 1.0);
    pinnaFreqL = juce::jlimit (1200.0, 16000.0, pinnaBaseHz - 0.5 * pinnaSpread * width);
    pinnaFreqR = juce::jlimit (1200.0, 16000.0, pinnaBaseHz + 0.5 * pinnaSpread * width);
    pinnaL.setNotch (fs, pinnaFreqL, pinnaQ, pinnaMinHz);
    pinnaR.setNotch (fs, pinnaFreqR, pinnaQ, pinnaMinHz);
    pinnaFiltersDirty = false;
}

void Mirror7Engine::setVelvetAmount (double amt)
{
    auto& noise = noiseController.getEngine();
    noise.aureo.velvet_on = true;
    noise.mode = aureo::NoiseColor::Velvet;
    noise.aureo.velvet_amount = clamp (amt, 0.0, 1.0);
}

void Mirror7Engine::setPanWidth (double pan, double width)
{
    panSmoothed.setTargetValue ((float) juce::jlimit (-1.0, 1.0, pan));
    params.width = clamp (width, 0.0, 1.0);
}

void Mirror7Engine::setITDILD (double itdUs, double ildDb)
{
    params.itd_us = clamp (itdUs, 0.0, 1600.0);
    params.ild_db = clamp (ildDb, 0.0, 30.0);
}

void Mirror7Engine::setAutopan (bool enable, double speed)
{
    autopan = enable; autopanSpeed = juce::jlimit (0.02, 2.0, speed);
}

int Mirror7Engine::findFreeGrain()
{
    for (int i = 0; i < aureo::kMaxGrains; ++i) if (!grains[(size_t)i].on) return i; return -1;
}

uint32_t Mirror7Engine::mapLenSamples (double u)
{
    const double base = clamp (params.baseLenMs, aureo::kMinBaseLengthMs, 2000.0) * 0.001 * fs;
    const double kexp = (2.0 * u - 1.0) * clamp01 (params.lenPhi);
    double L = base * std::pow (aureo::kPhi, kexp);
    L = clamp (L, aureo::kMinGrainSamples, fs * 4.0);
    return (uint32_t) L;
}

int Mirror7Engine::scheduleGapSamples()
{
    // Sync branch: quantized to host tempo division
    if (params.syncEnable && params.rateHz <= 1.0e-6 && tempo.current_rate > 1.0e-6)
    {
        const double gap_samples = std::max (1.0, fs / tempo.current_rate);
        return (int) std::round (gap_samples);
    }
    // Poisson around base rate
    double baseRate = clamp (params.rateHz, 0.001, aureo::kMaxEventRateHz);
    const double t = (fs > 0.0 ? (double) sample_counter / fs : 0.0);
    double lambda = baseRate * (1.0 + 0.2 * std::sin (2.0 * aureo::kPi * (t * (1.0 / aureo::kPhi))));
    lambda = std::max (lambda, 1e-3);
    const double U = std::max (1.0e-12, (double) rng.uni01());
    double gap_sec = -std::log (U) / lambda;
    const double cap_sec = std::min (30.0, 4.0 / std::max (baseRate, 1.0e-6));
    gap_sec = std::min (gap_sec, cap_sec);
    double gap_samples = std::max (1.0, gap_sec * fs);
    return (int) std::max (1.0, std::round (gap_samples));
}

double Mirror7Engine::applyBurstPosition (double pan, double dur_norm, double gap_norm, double& amp_scale)
{
    amp_scale = 1.0; // Simple burst mapping (no Hawkes yet)
    const double burst_w = 0.0; if (burst_w <= 1.0e-6) return pan;
    const double clamped_dur = clamp01 (dur_norm); const double clamped_gap = clamp01 (gap_norm);
    const double sign = (pan >= 0.0 ? 1.0 : -1.0); double mag = std::abs (pan);
    const double center_pull = burst_w * (0.35 + 0.45 * clamped_dur); mag = mag * (1.0 - center_pull) + center_pull * 0.12;
    const double edge_push = burst_w * (0.35 + 0.30 * (1.0 - clamped_dur)); mag = clamp (mag + edge_push * (1.0 - clamped_gap), 0.0, 1.0);
    amp_scale = 1.0 + burst_w * (0.6 + 0.4 * clamped_dur); return clamp (sign * mag, -1.0, 1.0);
}

Mirror7Engine::BlockDerivedParams Mirror7Engine::buildBlockDerivedParams() const
{
    BlockDerivedParams ctx;
    const double fsSafe = (fs > 0.0 ? fs : 44100.0);
    ctx.wowInc = map_phi_range (0.1, 1.5, clamp01 (params.vhsWow)) / fsSafe;
    ctx.flutterInc = map_phi_range (7.0, 12.0, clamp01 (params.vhsFlutter)) / fsSafe;
    ctx.itdScale = params.itd_us * 1.0e-6 * fsSafe;
    ctx.minPanNorm = clamp (params.spatMinDeg, 0.0, 180.0) / 90.0;
    ctx.ipdAmt = clamp01 (params.spatIPD);
    ctx.shadowAmt = clamp01 (params.spatShadow);
    ctx.phiDistance = clamp01 (params.phiDistance);
    ctx.phiPitchRad = 12.0 * (aureo::kPi / 180.0);
    return ctx;
}

void Mirror7Engine::process (juce::AudioBuffer<float>& buffer, float outGainDb)
{
    juce::dsp::AudioBlock<float> block (buffer);
    auto doProcess = [&] (juce::dsp::AudioBlock<float>& blk)
    {
        // Tempo update and transport
        const double block_dt = (fs > 0.0 ? (double) blk.getNumSamples() / fs : 0.0);
        updateTempo (block_dt);
        if (tempo.pending_start)
        {
            samples_to_next = 1; gap_elapsed = 0; sample_counter = 0; tempo.pending_start = false;
        }

        if (pinnaFiltersDirty) updatePinnaFilters();
        const bool pinna_now = params.pinnaOn;
        if (pinna_now != pinnaEnabled) { pinnaEnabled = pinna_now; pinnaL.clear(); pinnaR.clear(); }
        double pinnaMix = pinna_mix;
        const double pinnaTarget = (params.pinnaOn ? 1.0 : 0.0);
        const double depthAmt = juce::jlimit (0.0, 1.0, params.pinnaDepth / 24.0);
        double denom = fs * 0.02; if (denom < 1.0) denom = 1.0; const double pinnaStep = 1.0 / denom;
        // Update noise engine configuration per-block
        noiseController.configure (params, fs);
        auto* L = blk.getChannelPointer (0); auto* R = blk.getChannelPointer (1);
        const int N = (int) blk.getNumSamples();
        const auto blockParams = buildBlockDerivedParams();
        DialogueSettings dialogueSettings;
        dialogueSettings.enabled = params.dialogueOn;
        dialogueSettings.strength = params.dialogueStrength;
        dialogueSettings.memory = params.dialogueMemory;
        dialogueSettings.phiMix = params.dialoguePhiMix;
        const bool phiModeActive = params.phiMode && spatializer.hasValidGeometry();

        for (int n = 0; n < N; ++n)
        {
            ++gap_elapsed; lfo_wow += blockParams.wowInc; if (lfo_wow >= 1.0) lfo_wow -= 1.0; lfo_flt += blockParams.flutterInc; if (lfo_flt >= 1.0) lfo_flt -= 1.0;
            const double wow = std::sin (2.0 * aureo::kPi * lfo_wow); const double flt = std::sin (2.0 * aureo::kPi * lfo_flt); const double vhs_mod = 0.5 * wow + 0.5 * flt;

            // noise -> ring
            double nz = noiseController.process (rng, fs); ring.data[ring.writeIndex] = nz;

            if (--samples_to_next <= 0)
            {
                const int gi = findFreeGrain();
                if (gi >= 0)
                {
                    auto& g = grains[(size_t) gi]; g = {}; auto& cefg = cefgGrains[(size_t) gi]; cefg.tap_count = 0;
                    const double u1 = w_phi.next(); const double u2 = w_s2.next(); const double u3 = w_pl.next(); const double u4 = rng.uni01(); const double u5 = rng.uni01(); const double u6 = rng.uni01();
                    const double gap_samples = (double) gap_elapsed; const double prev_dur = std::max (1.0, (double) g.dur); const double ratio_gap = gap_samples / prev_dur;
                    const double coupling = clamp01 (params.hemisCoupling); const double time_weight = clamp (0.35 + 0.65 * (ratio_gap / (ratio_gap + 1.0)), 0.35, 1.0); const double hemi = coupling * time_weight;
                    const double amp_shape = std::pow (std::max (1e-9, u2), 0.35); const double amp = aureo::kAmpNorm * amp_shape;
                    const int base_dur_samples = (int) mapLenSamples (u1); const double base_dur = std::max (1.0, (double) base_dur_samples);
                    double pan = 2.0 * u3 - 1.0;
                    pan = clamp ((1.0 - hemi) * pan - hemi * prev_pan, -1.0, 1.0);
                    pan = spatializer.applyPhiPan (params.phiPan, pan);
                    const double timestamp_sec = (fs > 0.0 ? (double) sample_counter / fs : 0.0);
                    auto dres = dialogue.evaluate (dialogueSettings, timestamp_sec, pan, amp, base_dur, gap_samples);
                    pan = dres.pan; double amp_dialogue = amp * dres.amp_scale; double dur_dialogue = std::max (1.0, base_dur * dres.dur_scale);

                    // Poisson pan/time enforcement (simple)
                    auto enforce = [&](double& p){ p = clamp (p, -1.0, 1.0); for (const auto& og : grains) { if (!og.on) continue; const double dist = std::abs (p - og.pan); if (dist < blockParams.minPanNorm) { const double sign = (p >= og.pan) ? 1.0 : -1.0; p = clamp (og.pan + sign * blockParams.minPanNorm, -1.0, 1.0); } } };
                    enforce (pan);
                    spatializer.commitPhiState (params.phiPan, pan, true);

                    const int dur_samples = std::max (1, (int) std::round (dur_dialogue)); g.dur = (uint32_t) dur_samples;
                    const double sr_ref = (fs > 0.0 ? fs : 44100.0);
                    const double dur_norm = clamp ((double) g.dur / (sr_ref * 0.35), 0.0, 1.0); const double gap_norm = clamp (gap_samples / (sr_ref * 0.5), 0.0, 1.0);
                    double burst_amp_scale = 1.0; pan = applyBurstPosition (pan, dur_norm, gap_norm, burst_amp_scale); double amp_final = amp_dialogue * burst_amp_scale; if (params.spatMirror) pan = -pan;

                    g.amp = amp_final; g.pan = pan; double gL, gR; aureo::pan_equal_power (pan, params.width, gL, gR); g.panL = gL; g.panR = gR; g.gL = 1.0; g.gR = 1.0; g.focus = 1.0; g.crossfeed = 0.0; g.ipd_coeff = 0.0; g.ipd_zL = g.ipd_zR = 0.0; g.shadow_a = 0.0; g.shadow_zL = g.shadow_zR = 0.0; g.shadow_left = (pan > 0.0); g.shadow_right = (pan < 0.0); g.sr_holdN = 1; g.sr_holdCnt = 1; g.q_levels = 0;

                    bool have_binaural = true;
                    double ild_db = 0.0;
                    if (!phiModeActive)
                    {
                        const auto bin = aureo::compute_binaural (fs, pan, params.width, params.itd_us, params.ild_db, rng.uni01(), rng.uni01());
                        g.itd = bin.itd_samples;
                        ild_db = bin.ild_db;
                        g.focus = bin.focus;
                        g.crossfeed = bin.crossfeed;
                        g.gL = db_to_lin (ild_db);
                        g.gR = db_to_lin (-ild_db);
                    }
                    else
                    {
                        const auto& geom = spatializer.getGeometry();
                        auto head = aureonoise::spatial::phi::compute_head_result (geom, fs, pan, blockParams.phiPitchRad, blockParams.phiDistance);
                        g.itd = (1.0 - hemi) * head.itd_samples - hemi * prev_itd;
                        const double ildMax = clamp (params.ild_db, 0.0, 24.0);
                        ild_db = head.lateral * params.ild_db;
                        ild_db = clamp ((1.0 - hemi) * ild_db - hemi * prev_ild, -ildMax, ildMax);
                        g.gL = db_to_lin (ild_db);
                        g.gR = db_to_lin (-ild_db);
                        have_binaural = false;
                    }

                    double ipd_coeff = 0.0; if (blockParams.ipdAmt > 1.0e-6) { const double ipd_shape = std::abs (2.0 * u6 - 1.0); const double ipd_base = 0.18 + 0.55 * blockParams.ipdAmt; const double ipd_spread = 0.25 * blockParams.ipdAmt; ipd_coeff = clamp (ipd_base + ipd_spread * (ipd_shape - 0.5), 0.0, 0.95); }
                    ipd_coeff *= (0.7 + 0.3 * g.focus); g.ipd_coeff = ipd_coeff;

                    const double shadow_factor = blockParams.shadowAmt * std::abs (pan); if (shadow_factor > 1.0e-6) { const double fc_min = 800.0; const double fc_max = std::max (fc_min, std::min (8000.0, 0.45 * fs)); const double fc = map_phi_range (fc_min, fc_max, 1.0 - shadow_factor); const double alpha = std::exp (-2.0 * aureo::kPi * fc / fs); g.shadow_a = clamp (alpha, 0.0, 0.9999); g.shadow_left = (pan > 0.0); g.shadow_right = (pan < 0.0); }

                    // glitch kind
                    g.kind = aureo::choose_kind (params.glitchMix, u4);
                    int srN = aureo::map_sr_hold_base (params.srCrush, u1); g.sr_holdN = std::max (1, srN); g.sr_holdCnt = g.sr_holdN; const int bits = aureo::map_bits (params.bitCrush); g.q_levels = (1 << (bits - 1)) - 1;

                    const auto envShape = field.make_envelope (params.envAttack, params.envDecay, params.envSustain, params.envRelease, gap_samples, (double) g.dur, std::abs (pan)); g.env = envShape;
                    g.on = true; g.age = 0; prev_pan = pan; prev_itd = g.itd; prev_ild = ild_db; gap_elapsed = 0;

                    // Prepare CEFG taps for this grain
                    if (phiModeActive && params.spatialProfile == 1 && spatializer.ensureCefgLoaded())
                    {
                        const auto& geom = spatializer.getGeometry();
                        spatializer.prepareCefgGrain ((size_t) gi, params, geom, fs, blockParams.phiDistance, blockParams.phiPitchRad, cefgGrains, grains);
                    }

                    // commit dialogue
                    auto cres = dres;
                    cres.pan = pan;
                    const double base_amp_safe = std::max (1.0e-9, cres.base_amp);
                    cres.amp_scale = amp_final / base_amp_safe;
                    const double base_dur_safe = std::max (1.0, cres.base_dur);
                    cres.dur_scale = (double) g.dur / base_dur_safe;
                    dialogue.commit (dialogueSettings, cres, true);
                }
                samples_to_next = (uint32_t) scheduleGapSamples();
            }

            double yL = 0.0, yR = 0.0;
            for (auto& g : grains)
            {
                if (!g.on) continue; if (g.age >= g.dur) { g.on = false; continue; }
                const double env = field.envelope ((double) g.age / (double) std::max (1u, g.dur), g.env);
                const double xL = ring.data[ring.writeIndex];
                // simple sample-rate hold and bitcrush
                if (--g.sr_holdCnt <= 0) { g.heldL = xL; g.heldR = xL; g.sr_holdCnt = g.sr_holdN; }
                double sL = g.heldL * g.panL; double sR = g.heldR * g.panR; sL *= g.gL; sR *= g.gR;
                // ITD using per-channel delay
                const float delayL = (float) (g.itd > 0.0 ? g.itd : 0.0);
                const float delayR = (float) (g.itd < 0.0 ? -g.itd : 0.0);
                dL.setDelay (delayL);
                dR.setDelay (delayR);
                const float ydl = dL.popSample (0) + (float) sL;
                const float ydr = dR.popSample (0) + (float) sR;
                dL.pushSample (0, (float) sL);
                dR.pushSample (0, (float) sR);
                sL = ydl; sR = ydr;
                // ipd allpass
                if (g.ipd_coeff > 1.0e-6) { sL = allpass (sL, g.ipd_coeff, g.ipd_zL); sR = allpass (sR, -g.ipd_coeff, g.ipd_zR); }
                // shadow filtering
                if (g.shadow_a > 1.0e-6) { if (g.shadow_left) sL = shlp (sL, g.shadow_a, g.shadow_zL); if (g.shadow_right) sR = shlp (sR, g.shadow_a, g.shadow_zR); }

                // Early reflections (CEFG)
                auto& cefg = cefgGrains[&g - grains.data()];
                if (params.phiMode && cefg.tap_count > 0)
                {
                    double erL = 0.0, erR = 0.0;
                    for (int ti = 0; ti < cefg.tap_count; ++ti)
                    {
                        auto& t = cefg.taps[(size_t) ti];
                        const double smpL = ringReadDelay (ring, ring.writeIndex, t.base_delay + t.itd_l);
                        const double smpR = ringReadDelay (ring, ring.writeIndex, t.base_delay + t.itd_r);
                        t.lp_l = t.lp_a * t.lp_l + (1.0 - t.lp_a) * smpL;
                        t.lp_r = t.lp_a * t.lp_r + (1.0 - t.lp_a) * smpR;
                        erL += t.mix * t.gain_l * t.lp_l;
                        erR += t.mix * t.gain_r * t.lp_r;
                    }
                    sL += erL; sR += erR;
                }
                // glitch special cases
                if (g.kind == aureo::GrainKind::VhsDrop) { const double att = 0.5 + 0.5 * (1.0 - std::abs (vhs_mod)); sL *= att; sR *= att; }
                else if (g.kind == aureo::GrainKind::Stutter) { if ((g.age & 7) == 0) { sL *= 0.2; sR *= 0.2; } }
                // amp env
                const double amp_env = g.amp * env; yL += amp_env * sL; yR += amp_env * sR; ++g.age;
            }

            // modal (optional)
            const double modalMix = clamp (params.modalMix, 0.0, 1.0);
            if (params.modalOn && modalMix > 1.0e-6)
            {
                double mL = yL, mR = yR;
                modal.process (yL, yR, mL, mR);
                if (params.modalMirror) { mL = -mL; std::swap (mL, mR); }
                yL = (1.0 - modalMix) * yL + modalMix * mL;
                yR = (1.0 - modalMix) * yR + modalMix * mR;
                const double e = 0.5 * (mL * mL + mR * mR);
                modalEnergy = 0.995 * modalEnergy + 0.005 * e;
            }
            else { modalEnergy *= 0.995; }

            // pinna (global notch)
            if (pinnaEnabled)
            {
                const double dryL = yL, dryR = yR;
                double wetL = pinnaL.proc (dryL);
                double wetR = pinnaR.proc (dryR);
                const double depthAmt = juce::jlimit (0.0, 1.0, params.pinnaDepth / 24.0);
                wetL = dryL + depthAmt * (wetL - dryL);
                wetR = dryR + depthAmt * (wetR - dryR);
                const double wet = pinnaMix; const double dry = 1.0 - wet;
                yL = dry * dryL + wet * wetL;
                yR = dry * dryR + wet * wetR;
            }
            pinnaMix += (pinnaTarget - pinnaMix) * pinnaStep; pinnaMix = clamp01 (pinnaMix); pinna_mix = pinnaMix;

            // soft drive
            yL = std::tanh (yL * aureo::kOutDrive) / aureo::kOutDrive; yR = std::tanh (yR * aureo::kOutDrive) / aureo::kOutDrive;
            L[n] = (float) yL; R[n] = (float) yR; ring.writeIndex = (ring.writeIndex + 1u) & aureo::kRingMask; ++sample_counter;
        }
    };

    if (oversampling)
    {
        auto os = oversampling->processSamplesUp (block); doProcess (os); oversampling->processSamplesDown (block);
    }
    else
    {
        doProcess (block);
    }

    buffer.applyGain (juce::Decibels::decibelsToGain (outGainDb));
}
