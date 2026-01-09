#include "engine/DialogueSystem.h"

#include <algorithm>
#include <cmath>

using aureo::clamp;
using aureo::clamp01;
using aureo::map_phi_range;

void DialogueSystem::setSampleRate (double newFs)
{
    fs = (newFs > 0.0 ? newFs : 44100.0);
}

void DialogueSystem::reset()
{
    state = {};
    state.coherence = 1.0;
    state.fib_index = 1;
    state.fib_direction = 1;
    state.fib_ratio_target = aureo::kPhi;
}

DialogueSystem::DebugState DialogueSystem::getDebugState() const
{
    return { state.last_sign, state.utterances, state.coherence };
}

static double ratio_score (double ratio, double target)
{
    if (!std::isfinite (ratio) || ratio <= 0.0 || !std::isfinite (target) || target <= 0.0)
        return 0.0;
    const double log_diff = std::abs (std::log (ratio / target));
    constexpr double tol = 0.4;
    return std::clamp (1.0 - log_diff / tol, 0.0, 1.0);
}

DialogueResult DialogueSystem::evaluate (const DialogueSettings& settings,
                                         double timestamp_sec,
                                         double pan,
                                         double amp,
                                         double dur_samples,
                                         double gap_samples)
{
    DialogueResult res;
    res.pan = clamp (pan, -1.0, 1.0);
    res.base_amp = amp;
    res.base_dur = dur_samples;
    res.gap_samples = gap_samples;
    res.timestamp = timestamp_sec;
    res.sign = (res.pan >= 0.0 ? 1 : -1);

    if (!settings.enabled)
        return res;

    static constexpr double fibWalk[] = { 1, 1, 2, 3, 5, 8, 13, 21 };
    const int fibCount = (int) (sizeof (fibWalk) / sizeof (double));
    auto clampIdx = [fibCount](int i){ return std::clamp (i, 0, fibCount - 1); };
    int curIndex = clampIdx (state.fib_index);
    int curDir = (state.fib_direction >= 0 ? 1 : -1);
    int candDir = curDir;
    int rawIdx = curIndex + candDir;
    if (rawIdx < 0 || rawIdx >= fibCount)
    {
        candDir = -candDir;
        rawIdx = curIndex + candDir;
    }
    int candIndex = clampIdx (rawIdx);
    const double fibPrev = std::max (1.0, fibWalk[curIndex]);
    const double fibNext = std::max (1.0, fibWalk[candIndex]);
    double fibRatio = fibNext / fibPrev;
    if (!std::isfinite (fibRatio) || fibRatio <= 0.0)
        fibRatio = aureo::kPhi;
    res.handshake_ratio_target = fibRatio;
    res.fib_index = curIndex;
    res.fib_direction = curDir;

    const double strength = clamp01 (settings.strength);
    const double memory = clamp01 (settings.memory);
    const double phi_mix = clamp01 (settings.phiMix);
    const double mag = std::max (1.0e-6, std::abs (res.pan));
    const bool has_prev = (std::abs (state.last_sign) > 0.5);

    double coherence_target = aureo::kInvPhi;

    if (has_prev && (int) std::copysign (1.0, state.last_sign) != res.sign)
    {
        const double prev_mag = std::clamp (state.last_mag, 0.05, 1.0);
        const double incoming_mag = mag;
        const double incoming_ratio = incoming_mag / std::max (prev_mag, 1.0e-6);
        const double amp_ratio = std::max (1.0e-6, res.base_amp) / std::max (1.0e-6, state.last_amp);
        const double dur_ratio = (state.last_dur > 0.0) ? std::max (1.0e-6, dur_samples) / std::max (1.0e-6, state.last_dur)
                                                       : incoming_ratio;
        const double gap_ratio = (state.last_gap > 0.0) ? std::max (1.0e-6, gap_samples) / std::max (1.0e-6, state.last_gap)
                                                       : incoming_ratio;
        const double pan_score = ratio_score (incoming_ratio, fibRatio);
        const double amp_score = ratio_score (amp_ratio, fibRatio);
        const double dur_score = ratio_score (dur_ratio, fibRatio);
        const double gap_score = ratio_score (gap_ratio, fibRatio);
        res.handshake_score = 0.55 * pan_score + 0.15 * amp_score + 0.15 * dur_score + 0.15 * gap_score;
        const double handshake_threshold = std::clamp (0.42 - 0.18 * strength, 0.22, 0.42);
        res.handshake = (res.handshake_score >= handshake_threshold);
        const double target_mag = clamp ((1.0 - strength) * incoming_mag
                                         + strength * (phi_mix * std::clamp (prev_mag * fibRatio, 0.05, 1.0)
                                                       + (1.0 - phi_mix) * incoming_mag),
                                         0.05, 1.0);
        res.pan = (double) res.sign * target_mag;
        const double amp_target = std::max (1.0e-6, state.last_amp) * fibRatio;
        const double amp_scale_target = clamp (amp_target / std::max (1.0e-6, res.base_amp), 0.15, 4.0);
        const double amp_interp = phi_mix * amp_scale_target + (1.0 - phi_mix);
        res.amp_scale = clamp ((1.0 - strength) * 1.0 + strength * amp_interp, 0.15, 3.0);
        if (state.last_dur > 0.0)
        {
            const double dur_target = state.last_dur * fibRatio;
            const double dur_scale_target = clamp (dur_target / std::max (1.0, res.base_dur), 0.25, 4.0);
            const double dur_interp = phi_mix * dur_scale_target + (1.0 - phi_mix);
            res.dur_scale = clamp ((1.0 - strength) * 1.0 + strength * dur_interp, 0.25, 4.0);
        }
        coherence_target = std::clamp (fibRatio, 0.6, 1.8);
        if (res.handshake)
        {
            res.fib_index = candIndex;
            res.fib_direction = candDir;
        }
    }
    else if (has_prev)
    {
        const double relaxed_target = std::clamp (state.last_mag * state.fib_ratio_target, 0.05, 1.0);
        const double relaxed = (1.0 - 0.5 * strength) * mag + 0.5 * strength * relaxed_target;
        res.pan = (double) res.sign * clamp (relaxed, 0.05, 1.0);
        coherence_target = 1.0;
    }
    else
    {
        coherence_target = 1.0;
    }

    if (state.pan_samples > 0)
    {
        const double bias = state.pan_sum / (double) state.pan_samples;
        const double corr = strength * 0.4;
        res.pan = clamp (res.pan - corr * bias, -1.0, 1.0);
    }
    res.sign = (res.pan >= 0.0 ? 1 : -1);
    res.next_coherence = (1.0 - memory) * state.coherence + memory * coherence_target;
    return res;
}

void DialogueSystem::commit (const DialogueSettings& settings,
                             const DialogueResult& res,
                             bool accepted)
{
    if (!settings.enabled)
        return;
    if (!accepted)
    {
        ++state.failures;
        return;
    }
    const double prev_gap = state.last_gap;
    state.last_sign = (double) res.sign;
    state.last_mag = std::abs (res.pan);
    state.last_amp = std::max (1.0e-6, res.base_amp * res.amp_scale);
    state.last_dur = res.base_dur * res.dur_scale;
    state.last_gap = res.gap_samples;
    state.last_timestamp = res.timestamp;
    state.coherence = res.next_coherence;
    state.coherence_accum += state.coherence;
    ++state.coherence_samples;
    ++state.utterances;
    state.fib_index = std::clamp (res.fib_index, 0, 7);
    state.fib_direction = (res.fib_direction >= 0 ? 1 : -1);
    if (res.handshake)
    {
        ++state.handshakes;
        state.fib_ratio_target = res.handshake_ratio_target;
        state.fib_gap_queue.fill (0.0);
        state.fib_gap_queue_pos = 0;
        state.fib_gap_queue_size = 0;
        const double gap_cap = std::max (1.0, fs * 16.0);
        auto push_gap = [&](double g)
        {
            if (state.fib_gap_queue_size >= state.fib_gap_queue.size())
                return;
            const double clamped = clamp (g, 1.0, gap_cap);
            state.fib_gap_queue[state.fib_gap_queue_size++] = clamped;
        };
        const double base_gap = (prev_gap > 1.0 ? prev_gap : std::max (1.0, res.gap_samples));
        push_gap (base_gap);
        const double phi_lo = std::max (1.0e-6, res.handshake_ratio_target / aureo::kPhi);
        const double phi_hi = std::max (phi_lo * 1.000001, res.handshake_ratio_target * aureo::kPhi);
        const double phi_blend = map_phi_range (phi_lo, phi_hi, clamp01 (res.handshake_score));
        const double plastic_bias = std::pow (aureo::kPlastic, 0.05 * (res.handshake_score - 0.5));
        const double extend_ratio = clamp (phi_blend * plastic_bias, 0.25, 8.0);
        push_gap (base_gap * extend_ratio);
    }
    else
    {
        if (state.fib_gap_queue_pos >= state.fib_gap_queue_size)
        {
            state.fib_gap_queue_size = 0;
            state.fib_gap_queue_pos = 0;
            state.fib_gap_queue.fill (0.0);
        }
    }
    state.fib_score = res.handshake_score;
    state.pan_sum += res.pan;
    state.pan_abs_sum += std::abs (res.pan);
    ++state.pan_samples;
}
