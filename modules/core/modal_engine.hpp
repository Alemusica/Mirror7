#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "aureo_core/aureo_math.hpp"

namespace aureonoise::modal {

struct ModeDesc {
  double freq;  // reference frequency (for ratio derivation)
  double decay;
  double gain;
};

inline std::vector<ModeDesc> preset_modes(int preset)
{
  switch (preset) {
    case 1: // wood
      return {
        {205.0, 1.60, 1.00}, {312.0, 1.45, 0.88}, {421.0, 1.32, 0.78},
        {538.0, 1.18, 0.70}, {660.0, 1.00, 0.60}, {795.0, 0.88, 0.52},
        {932.0, 0.80, 0.44}, {1105.0,0.72, 0.36}
      };
    case 2: // metal
      return {
        {185.0, 2.40, 0.95}, {260.0, 2.20, 0.90}, {340.0, 2.00, 0.85},
        {420.0, 1.85, 0.80}, {520.0, 1.70, 0.74}, {630.0, 1.55, 0.68},
        {760.0, 1.40, 0.62}, {905.0, 1.25, 0.56}, {1085.0,1.15, 0.48},
        {1265.0,1.05, 0.42}
      };
    case 3: // glassy
      return {
        {480.0, 1.80, 0.80}, {720.0, 1.60, 0.72}, {960.0, 1.44, 0.64},
        {1200.0,1.28, 0.58}, {1440.0,1.14, 0.52}, {1680.0,1.02, 0.46},
        {1920.0,0.92, 0.40}, {2160.0,0.84, 0.35}, {2400.0,0.76, 0.30}
      };
    default:
      return {};
  }
}

class Engine {
public:
  void setSampleRate(double sr)
  {
    if (sr <= 0.0) sr = 44100.0;
    if (std::abs(sr - sr_) > 1.0) {
      sr_ = sr;
      rebuildCoefficients();
    }
  }

  void configure(int preset, double decayNorm)
  {
    preset = std::clamp(preset, 0, 3);
    decayNorm = aureo::clamp(decayNorm, 0.0, 1.0);
    bool presetChanged = (preset_ != preset);
    preset_ = preset;
    decayNorm_ = decayNorm;
    if (presetChanged) {
      buildModes();
    } else {
      rebuildCoefficients();
    }
  }

  void setBaseHz(double baseHz)
  {
    baseHz = std::clamp(baseHz, 20.0, 20000.0);
    if (std::abs(baseHz - baseHz_) > 1.0e-6) {
      baseHz_ = baseHz;
      rebuildCoefficients();
    }
  }

  void setTransposeSemitones(double semis)
  {
    if (std::abs(semis - transpose_semitones_) > 1.0e-9) {
      transpose_semitones_ = semis;
      freq_scale_ = std::pow(2.0, transpose_semitones_ / 12.0);
      rebuildCoefficients();
    }
  }

  // Switch between preset table and harmonic bank
  void setBankActive(bool enable)
  {
    if (enable != (source_ == Source::Bank)) {
      source_ = enable ? Source::Bank : Source::Preset;
      buildModes();
    }
  }

  void setHarmonicBank(int count, double stretch, double tilt, double decaySlope, bool oddOnly)
  {
    bank_count_ = std::clamp(count, 1, 128);
    bank_stretch_ = aureo::clamp(stretch, -1.0, 1.0);
    bank_tilt_ = aureo::clamp(tilt, -1.0, 1.0);
    bank_decay_slope_ = std::max(0.0, decaySlope);
    bank_oddOnly_ = oddOnly;
    source_ = Source::Bank;
    buildModes();
  }

  void reset()
  {
    for (auto& m : modes_) {
      m.y1L = m.y2L = 0.0;
      m.y1R = m.y2R = 0.0;
    }
  }

  void process(double inL, double inR, double& outL, double& outR)
  {
    if (modes_.empty()) {
      outL = inL;
      outR = inR;
      return;
    }

    double sumL = 0.0;
    double sumR = 0.0;
    for (auto& m : modes_) {
      double yL = m.b0 * inL + m.a1 * m.y1L + m.a2 * m.y2L;
      double yR = m.b0 * inR + m.a1 * m.y1R + m.a2 * m.y2R;
      m.y2L = m.y1L;
      m.y1L = yL;
      m.y2R = m.y1R;
      m.y1R = yR;
      sumL += m.gain * yL;
      sumR += m.gain * yR;
    }

    outL = sumL;
    outR = sumR;
  }

private:
  enum class Source { Preset, Bank };
  Source source_ = Source::Preset;
  struct ModeState {
    double ratio = 1.0; // relative to baseHz_
    double decay = 1.0;
    double gain = 1.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double b0 = 0.0;
    double y1L = 0.0;
    double y2L = 0.0;
    double y1R = 0.0;
    double y2R = 0.0;
  };

  double sr_ = 44100.0;
  int    preset_ = 0;
  double decayNorm_ = 0.6;
  double baseHz_ = 220.0;              // tunable base frequency
  double transpose_semitones_ = 0.0;   // additional pitch shift
  double freq_scale_ = 1.0;            // = 2^(semitones/12)
  std::vector<ModeState> modes_;

  // Harmonic bank parameters
  int    bank_count_ = 0;
  double bank_stretch_ = 0.0;   // [-1..1]
  double bank_tilt_ = 0.0;      // [-1..1] spectral tilt
  double bank_decay_slope_ = 0.5;// >=0, more -> faster decay for higher modes
  bool   bank_oddOnly_ = false;

  void buildModes()
  {
    if (source_ == Source::Preset) {
      buildPresetModes();
    } else {
      buildBankModes();
    }
  }

  void buildPresetModes()
  {
    modes_.clear();
    if (preset_ <= 0) return;
    auto table = preset_modes(preset_);
    modes_.reserve(table.size());
    double f0 = 0.0;
    if (!table.empty()) f0 = std::max(1.0, table.front().freq);
    for (const auto& desc : table) {
      ModeState st;
      const double ratio = (f0 > 0.0) ? (desc.freq / f0) : 1.0;
      st.ratio = (std::isfinite(ratio) && ratio > 0.0) ? ratio : 1.0;
      st.decay = desc.decay;
      st.gain = desc.gain;
      modes_.push_back(st);
    }
    rebuildCoefficients();
  }

  void buildBankModes()
  {
    modes_.clear();
    if (bank_count_ <= 0) return;
    modes_.reserve(static_cast<size_t>(bank_count_));
    const double s = aureo::clamp(bank_stretch_, -1.0, 1.0);
    const double tilt = aureo::clamp(bank_tilt_, -1.0, 1.0);
    const double powTilt = -0.75 * tilt; // +tilt -> darker
    double sumG = 0.0;
    for (int k = 0; k < bank_count_; ++k) {
      const int n = bank_oddOnly_ ? (2 * k + 1) : (k + 1);
      const double ratio = std::pow(static_cast<double>(n), 1.0 + 0.6 * s);
      ModeState st;
      st.ratio = std::max(1.0e-6, ratio);
      st.decay = 1.0 / std::pow(static_cast<double>(n), std::max(0.0, bank_decay_slope_));
      st.gain = std::pow(static_cast<double>(n), powTilt);
      modes_.push_back(st);
      sumG += std::abs(st.gain);
    }
    if (sumG > 1.0e-12) {
      const double norm = 1.0 / sumG; // keep energy controlled
      for (auto& m : modes_) m.gain *= norm;
    }
    rebuildCoefficients();
  }

  void rebuildCoefficients()
  {
    if (modes_.empty()) return;
    for (auto& m : modes_) {
      const double target = baseHz_ * freq_scale_ * m.ratio;
      const double freq = std::clamp(target, 40.0, 0.45 * sr_);
      double T60 = std::max(0.05, m.decay * (0.4 + 0.6 * decayNorm_));
      const double r = std::exp(-6.907755278982137 / (T60 * sr_));
      const double omega = 2.0 * aureo::kPi * freq / sr_;
      m.b0 = (1.0 - r);
      m.a1 = 2.0 * r * std::cos(omega);
      m.a2 = -r * r;
      m.y1L = m.y2L = 0.0;
      m.y1R = m.y2R = 0.0;
    }
  }
};

} // namespace aureonoise::modal
