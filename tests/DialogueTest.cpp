#include "engine/DialogueSystem.h"
#include "aureo_math.hpp"

#include <iostream>

int main()
{
    DialogueSystem system;
    system.setSampleRate (48000.0);
    system.reset();

    DialogueSettings settings;
    settings.enabled = true;
    settings.strength = 1.0;
    settings.memory = 0.5;
    settings.phiMix = 1.0;

    auto first = system.evaluate (settings, 0.0, 0.45, 1.0, 1200.0, 2400.0);
    system.commit (settings, first, true);

    const double phi = aureo::kPhi;
    auto second = system.evaluate (settings, 0.5, -0.45, 1.0 * phi, 1200.0 * phi, 2400.0 * phi);
    if (second.sign >= 0)
    {
        std::cerr << "[mirror7] Dialogue alternation failed to flip sign.\n";
        return 1;
    }
    if (second.handshake_score < 0.2)
    {
        std::cerr << "[mirror7] Dialogue coherence score too low for golden-ratio pairings.\n";
        return 1;
    }
    system.commit (settings, second, true);

    const auto dbg = system.getDebugState();
    if (dbg.lastSign >= 0.0)
    {
        std::cerr << "[mirror7] Dialogue state did not preserve negative sign after commit.\n";
        return 1;
    }

    return 0;
}
