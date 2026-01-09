#pragma once

#include <array>

#include "aureo_math.hpp"

struct DialogueSettings
{
    bool enabled = true;
    double strength = 0.6;
    double memory = 0.5;
    double phiMix = 0.75;
};

struct DialogueResult
{
    double pan = 0.0;
    double amp_scale = 1.0;
    double dur_scale = 1.0;
    double next_coherence = 1.0;
    double base_amp = 1.0;
    double base_dur = 0.0;
    double gap_samples = 0.0;
    double timestamp = 0.0;
    double handshake_score = 0.0;
    double handshake_ratio_target = aureo::kPhi;
    bool handshake = false;
    int sign = 1;
    int fib_index = 1;
    int fib_direction = 1;
};

class DialogueSystem
{
public:
    DialogueSystem() = default;

    void setSampleRate (double newFs);
    void reset();

    DialogueResult evaluate (const DialogueSettings& settings,
                             double timestamp_sec,
                             double pan,
                             double amp,
                             double dur_samples,
                             double gap_samples);

    void commit (const DialogueSettings& settings,
                 const DialogueResult& res,
                 bool accepted);

    struct DebugState
    {
        double lastSign = 0.0;
        uint64_t utterances = 0;
        double coherence = 1.0;
    };

    DebugState getDebugState() const;

private:
    struct State
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
        double fib_ratio_target = aureo::kPhi;
        double fib_score = 0.0;
        std::array<double, 3> fib_gap_queue{};
        uint32_t fib_gap_queue_size = 0;
        uint32_t fib_gap_queue_pos = 0;
    };

    double fs = 44100.0;
    State state;
};
