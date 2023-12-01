#pragma once

namespace nusyst {

///
///\note Motivated by McFarland-Day calculation Phys. Rev. D 86, 053003 (2012)
/// Fig 6.
inline double GetNueNueBarXSecRatioWeight(int nu_pdg, int is_CC,
                                          double nu_Energy_GeV,
                                          double parameter_value = 1) {

  static double const A_nu = 0.01;
  static double const A_nubar = -0.018;
  static double const B_nu = -1;
  static double const B_nubar = -1;
  static double const E_min = 0.2;

  static double const central_value = 0;
  static double const uncertainty = 1;

  if (!is_CC) {
    return 1;
  }

  if (nu_pdg == 12) {
    return 1.0 / (1.0 + (central_value + parameter_value * uncertainty) *
                            (A_nu * pow(std::max(E_min, nu_Energy_GeV), B_nu)));
  }

  if (nu_pdg == -12) {
    return 1.0 /
           (1.0 + (central_value + parameter_value * uncertainty) *
                      (A_nubar * pow(std::max(E_min, nu_Energy_GeV), B_nubar)));
  }

  return 1;
}
} // namespace nusyst
