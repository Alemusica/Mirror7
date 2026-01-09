#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

#include "aureo_math.hpp"
#include "aureonoise/spatial/phi_model.hpp"
#include "aureonoise/spatial/cefg.hpp"
#include "aureo_ring.hpp"
#include "aureo_field_adsr.hpp"
#include "aureo_binaural.hpp"

inline std::array<double, 3> dirFromAngles (double az_deg, double el_deg)
{
    const double az = az_deg * (aureo::kPi / 180.0);
    const double el = el_deg * (aureo::kPi / 180.0);
    const double cosEl = std::cos (el);
    return { cosEl * std::sin (az), std::sin (el), cosEl * std::cos (az) };
}

struct SpatialParamsView
{
    bool phiMode = true;
    double phiHeadB = 0.0875;
    double phiRatio = aureo::kPhi;
    double cefgGain = 1.1;
    double cefgMix = 3.4;
    int spatialProfile = 0;
};

struct CefgTap
{
    double base_delay = 0;
    double itd_l = 0;
    double itd_r = 0;
    double base_gain_l = 0;
    double base_gain_r = 0;
    double gain_l = 0;
    double gain_r = 0;
    double lp_a = 0;
    double lp_l = 0;
    double lp_r = 0;
    double mix = 0;
};

struct CefgGrainState
{
    int tap_count = 0;
    std::array<CefgTap, 16> taps{};
};

class Spatializer
{
public:
    explicit Spatializer (std::string assetsRoot);

    void reset (double sampleRate);
    void updateGeometry (const SpatialParamsView& params);

    double applyPhiPan (bool enabled, double pan);
    void commitPhiState (bool enabled, double pan, bool ok);

    bool hasValidGeometry() const { return phiGeomValid; }
    const aureonoise::spatial::phi::Geometry& getGeometry() const { return phiGeom; }

    bool ensureCefgLoaded();

    template <typename Params, typename GrainArray>
    void prepareCefgGrain (size_t grainIndex,
                           const Params& params,
                           const aureonoise::spatial::phi::Geometry& geom,
                           double sampleRate,
                           double phi_distance,
                           double pitch_rad,
                           std::array<CefgGrainState, aureo::kMaxGrains>& states,
                           const GrainArray& grains)
    {
        if (!cefgLoaded || params.spatialProfile != 1)
            return;

        auto& state = states[grainIndex];
        state.tap_count = 0;
        const double mixScale = std::clamp (params.cefgMix, 0.0, 4.0);
        if (mixScale <= 1.0e-6)
            return;

        const double az_guess = grains[grainIndex].pan * 90.0;
        const double distance_m = 0.8 + 3.2 * std::clamp (phi_distance, 0.0, 1.0);
        const double log_distance = std::log10 (std::max (0.12, distance_m));
        const auto* cell = cefgProfile.find (az_guess, 0.0, log_distance);
        if (!cell)
            return;

        const double gainScale = std::clamp (params.cefgGain, 0.0, 4.0);
        const double fc_base = 2200.0 + 1600.0 * (1.0 - std::clamp (phi_distance, 0.0, 1.0));
        const double alpha_base = std::exp (-2.0 * aureo::kPi * fc_base / sampleRate);
        int tapCount = 0;
        for (const auto& ref : cell->reflections)
        {
            if (tapCount >= (int) state.taps.size())
                break;
            const double baseGain = ref.gain;
            const double scaledGain = baseGain * gainScale;
            if (std::abs (scaledGain) < 1.0e-5)
                continue;
            const auto dir = dirFromAngles (ref.azimuth_deg, ref.elevation_deg);
            const auto refHead = aureonoise::spatial::phi::compute_head_result_from_dir (geom, sampleRate, dir, pitch_rad, phi_distance);
            const double refPan = std::clamp (refHead.lateral, -1.0, 1.0);
            double gL, gR;
            aureo::pan_equal_power (refPan, 1.0, gL, gR);
            const double ildPush = 1.0 + 0.18 * refHead.lateral;
            auto& tap = state.taps[(size_t) tapCount++];
            tap.base_delay = std::max (0.0, ref.delay_sec * sampleRate);
            tap.itd_l = (refHead.itd_samples < 0.0 ? -refHead.itd_samples : 0.0);
            tap.itd_r = (refHead.itd_samples > 0.0 ? refHead.itd_samples : 0.0);
            tap.base_gain_l = baseGain * gL / ildPush;
            tap.base_gain_r = baseGain * gR * ildPush;
            tap.gain_l = tap.base_gain_l * gainScale;
            tap.gain_r = tap.base_gain_r * gainScale;
            tap.lp_a = std::clamp (alpha_base, 0.0, 0.9995);
            tap.lp_l = 0.0;
            tap.lp_r = 0.0;
            tap.mix = mixScale;
        }
        state.tap_count = tapCount;
    }

private:
    std::string assetsRoot;
    double sampleRate = 44100.0;
    aureonoise::spatial::phi::Geometry phiGeom{};
    bool phiGeomValid = false;
    double phi_last_mag = 0.0;
    int phi_last_sign = 1;
    aureonoise::spatial::cefg::Profile cefgProfile;
    bool cefgLoaded = false;
};
