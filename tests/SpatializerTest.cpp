#include "engine/Mirror7Engine.h"
#include "aureo_math.hpp"

#include <cmath>
#include <iostream>

int main()
{
    Spatializer spatializer (AUREONOISE_SPATIAL_ASSETS_DIR);
    spatializer.reset (48000.0);

    Mirror7Engine::Params params;
    params.phiMode = true;
    params.phiRatio = aureo::kPhi;

    SpatialParamsView view { params.phiMode, params.phiHeadB, params.phiRatio, params.cefgGain, params.cefgMix, params.spatialProfile };
    spatializer.updateGeometry (view);

    if (!spatializer.hasValidGeometry())
    {
        std::cerr << "[mirror7] Spatializer failed to build Φ-geometry.\n";
        return 1;
    }

    const auto& geom = spatializer.getGeometry();
    const double pitch = 12.0 * (aureo::kPi / 180.0);
    const double distance = 0.25;
    auto left = aureonoise::spatial::phi::compute_head_result (geom, 48000.0, -0.35, pitch, distance);
    auto right = aureonoise::spatial::phi::compute_head_result (geom, 48000.0, 0.35, pitch, distance);

    if (std::abs (left.itd_samples + right.itd_samples) > 1.0e-6)
    {
        std::cerr << "[mirror7] Φ-mode ITD continuity failed (not mirrored).\n";
        return 1;
    }

    if (std::abs (left.lateral + right.lateral) > 1.0e-3)
    {
        std::cerr << "[mirror7] Φ-mode ILD lateral asymmetry detected.\n";
        return 1;
    }

    return 0;
}
