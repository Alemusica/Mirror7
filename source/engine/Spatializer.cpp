#include "engine/Spatializer.h"

#include <cmath>

using aureo::clamp;

Spatializer::Spatializer (std::string assets)
    : assetsRoot (std::move (assets))
{
}

void Spatializer::reset (double newSampleRate)
{
    sampleRate = (newSampleRate > 0.0 ? newSampleRate : 44100.0);
    phi_last_mag = 0.0;
    phi_last_sign = 1;
    cefgLoaded = false;
}

void Spatializer::updateGeometry (const SpatialParamsView& params)
{
    if (params.phiMode)
    {
        aureonoise::spatial::phi::GeometryConfig cfg;
        cfg.head_b = clamp (params.phiHeadB, 0.06, 0.11);
        cfg.phi_ratio = std::max (1.0, params.phiRatio);
        phiGeom = aureonoise::spatial::phi::build_geometry (cfg);
        phiGeomValid = true;
    }
    else
    {
        phiGeomValid = false;
    }
}

double Spatializer::applyPhiPan (bool enabled, double pan)
{
    if (!enabled)
        return clamp (pan, -1.0, 1.0);
    if (phi_last_mag > 1.0e-6)
    {
        const int next_sign = (phi_last_sign >= 0) ? -1 : 1;
        double target_mag = phi_last_mag * aureo::kPhi;
        if (target_mag > 0.999)
        {
            const double fallback = phi_last_mag / aureo::kPhi;
            if (fallback >= 0.05)
                target_mag = fallback;
            else
                target_mag = std::max (0.05, std::abs (pan));
        }
        target_mag = clamp (target_mag, 0.05, 1.0);
        return clamp ((double) next_sign * target_mag, -1.0, 1.0);
    }
    const double base_mag = std::max (0.05, std::abs (pan));
    const int first_sign = (pan >= 0.0 ? 1 : -1);
    return clamp ((double) first_sign * base_mag, -1.0, 1.0);
}

void Spatializer::commitPhiState (bool enabled, double pan, bool ok)
{
    if (!enabled)
        return;
    const double mag = std::abs (pan);
    if (!ok && mag < 1.0e-6)
    {
        phi_last_mag = 0.0;
        phi_last_sign = 1;
        return;
    }
    phi_last_mag = mag;
    phi_last_sign = (pan >= 0.0 ? 1 : -1);
}

bool Spatializer::ensureCefgLoaded()
{
    if (cefgLoaded)
        return true;
    std::string path = assetsRoot + "/cefg/uteroid_phi.cefg";
    cefgLoaded = cefgProfile.load (path);
    return cefgLoaded;
}
