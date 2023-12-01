#pragma once

namespace nusyst {
/// N.B. This is bad. It is a hard copy of simb::int_type from
/// LArSoft/nusimdata/SimulationBase/MCNeutrino.h It means that standalone
/// nusystematics has access to this enum as well
enum class simb_mode_copy {
  kUnknownInteraction = -1,
  kQE = 0,
  kRes = 1,
  kDIS = 2,
  kCoh = 3,
  kCohElastic = 4,
  kElectronScattering = 5,
  kIMDAnnihilation = 6,
  kInverseBetaDecay = 7,
  kGlashowResonance = 8,
  kAMNuGamma = 9,
  kMEC = 10,
  kDiffractive = 11,
  kEM = 12,
  kWeakMix = 13,
};

inline std::string tostr(simb_mode_copy const &c) {
  switch (c) {
  case simb_mode_copy::kUnknownInteraction: {
    return "kUnknownInteraction";
  }
  case simb_mode_copy::kQE: {
    return "kQE";
  }
  case simb_mode_copy::kRes: {
    return "kRes";
  }
  case simb_mode_copy::kDIS: {
    return "kDIS";
  }
  case simb_mode_copy::kCoh: {
    return "kCoh";
  }
  case simb_mode_copy::kCohElastic: {
    return "kCohElastic";
  }
  case simb_mode_copy::kElectronScattering: {
    return "kElectronScattering";
  }
  case simb_mode_copy::kIMDAnnihilation: {
    return "kIMDAnnihilation";
  }
  case simb_mode_copy::kInverseBetaDecay: {
    return "kInverseBetaDecay";
  }
  case simb_mode_copy::kGlashowResonance: {
    return "kGlashowResonance";
  }
  case simb_mode_copy::kAMNuGamma: {
    return "kAMNuGamma";
  }
  case simb_mode_copy::kMEC: {
    return "kMEC";
  }
  case simb_mode_copy::kDiffractive: {
    return "kDiffractive";
  }
  case simb_mode_copy::kEM: {
    return "kEM";
  }
  case simb_mode_copy::kWeakMix: {
    return "kWeakMix";
  }
  }
  return "kUnknownInteraction";
}
} // namespace nusyst
