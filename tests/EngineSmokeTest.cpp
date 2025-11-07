#include "engine/Mirror7Engine.h"

#include <juce_dsp/juce_dsp.h>

#include <cmath>
#include <iostream>

namespace
{
bool bufferIsFinite (const juce::AudioBuffer<float>& buffer)
{
    for (int ch = 0; ch < buffer.getNumChannels(); ++ch)
    {
        const float* data = buffer.getReadPointer (ch);
        for (int i = 0; i < buffer.getNumSamples(); ++i)
        {
            if (!std::isfinite (data[i]))
                return false;
        }
    }
    return true;
}
} // namespace

int main()
{
    Mirror7Engine engine;
    engine.setAutopan (false, 0.2);

    juce::AudioBuffer<float> buffer (2, 256);
    buffer.clear();

    engine.prepare (48000.0, 512);
    engine.process (buffer, 0.0f);

    if (!bufferIsFinite (buffer))
    {
        std::cerr << "[mirror7] Engine produced NaN/Inf samples at 1x oversampling.\n";
        return 1;
    }

    engine.setOversampling (2); // request 4x OS
    buffer.clear();
    engine.process (buffer, -3.0f);

    if (!bufferIsFinite (buffer))
    {
        std::cerr << "[mirror7] Engine produced NaN/Inf samples at 4x oversampling.\n";
        return 1;
    }

    return 0;
}
