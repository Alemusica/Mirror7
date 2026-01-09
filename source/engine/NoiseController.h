#pragma once

#include "aureo_noise.hpp"
#include "aureo_math.hpp"

class NoiseController
{
public:
    void reset();

    template <typename Params>
    void configure (const Params& params, double sampleRate)
    {
        using aureo::clamp;
        using aureo::clamp01;
        engine.mode = (aureo::NoiseColor) params.noiseMode;
        engine.classic.color = (aureo::NoiseColor) params.noiseColor;
        engine.classic.amount = clamp01 (params.colorAmt);
        const bool quantum_mode = (engine.mode == aureo::NoiseColor::Quantum);
        if (quantum_mode)
        {
            const double mix_q = clamp01 (params.qMix);
            const double detail_q = clamp01 (params.qDetail);
            const double stride_q = std::max (0.1, params.qStride);
            const double base_q = clamp (params.qBase, 20.0, 2000.0);
            const int maxP = aureo::NoiseAureoState::kMaxPartials;
            const int minP = 6;
            int target = (int) std::round (minP + detail_q * (maxP - minP));
            target = std::clamp (target, 4, maxP);
            engine.aureo.phi_pi_mix = mix_q;
            engine.aureo.planck_decay = aureo::map_phi_range (0.12, 0.52, 1.0 - detail_q);
            engine.aureo.prime_stride = stride_q;
            engine.aureo.active_partials = target;
            engine.aureo.base_freq = base_q;
            engine.aureo.velvet_on = (params.qVelvet > 1.0e-4);
            engine.aureo.velvet_amount = clamp (params.qVelvet, 0.0, 1.0) * 0.3;
        }
        else if (engine.mode == aureo::NoiseColor::Aureo)
        {
            engine.aureo.phi_pi_mix = clamp01 (params.aureoMix);
            engine.aureo.planck_decay = clamp01 (params.aureoDecay);
            engine.aureo.prime_stride = std::max (0.1, params.aureoStride);
            engine.aureo.active_partials = (int) std::clamp (params.aureoHarmonics, 1.0, (double) aureo::NoiseAureoState::kMaxPartials);
            engine.aureo.base_freq = std::max (10.0, sampleRate * 0.25 * aureo::kInvSqrt2);
            engine.aureo.velvet_on = params.aureoVelvet;
            engine.aureo.velvet_amount = 0.10;
        }
    }

    double process (aureo::RNG& rng, double sampleRate);

    aureo::AureoNoiseEngine& getEngine() { return engine; }
    const aureo::AureoNoiseEngine& getEngine() const { return engine; }

private:
    aureo::AureoNoiseEngine engine;
};
