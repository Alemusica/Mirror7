#include "engine/NoiseController.h"

void NoiseController::reset()
{
    engine.reset();
}

double NoiseController::process (aureo::RNG& rng, double sampleRate)
{
    auto sample = engine.process (rng, sampleRate);
    return aureo::soft_tanh (sample * 1.2);
}
