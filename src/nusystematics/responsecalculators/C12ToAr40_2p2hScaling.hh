#pragma once

#include "nusystematics/utility/GENIEUtils.hh"

namespace nusyst {

inline double GetC_Ar2p2hScalingWeight(double parameter_value = 0) {
  static double const central_value = 1;
  static double const uncertainty = 0.2;
  return (central_value + parameter_value * uncertainty);
}

} // namespace nusyst
