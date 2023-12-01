#pragma once

#include "nusystematics/utility/GENIEUtils.hh"

namespace nusyst {

static double const A_OneSig = 0.9; // GeV2
static double const B_OneSig = 0.3; // GeV

static double const A_CV = 0;
static double const B_CV = 0;

// #define MINERVAE2p2h_DEBUG

///\brief Get an ad hoc factor to add energy dependence to 2p2h below the level
/// of MINERvA resolution.
///
/// MINERvA sees no energy dependence at the 10% level between 3 GeV and 6 GeV
/// Want at Enu = 3 GeV, a normalization weight of < 1.1
///
inline double Get_MINERvA2p2h2EnergyDependencyScaling(int mode_simb, int is_CC,
                                                      double nu_Energy_GeV,
                                                      double A_val = 0,
                                                      double B_val = 0) {

  if (!is_CC) {
    return 1;
  }

  if (mode_simb != e2i(simb_mode_copy::kMEC)) {
    return 1;
  }

  double A = A_CV + A_val * A_OneSig;
  double B = B_CV + B_val * B_OneSig;

#ifdef MINERVAE2p2h_DEBUG
  std::cout << "A = " << A << ", A_val = " << A_val << std::endl;
  std::cout << "B = " << B << ", B_val = " << B_val << std::endl;
#endif

  return 1.0 / (1.0 + (A / pow(nu_Energy_GeV, 2)) + (B / nu_Energy_GeV));
}

} // namespace nusyst