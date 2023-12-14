#pragma once

namespace nusyst {

inline std::pair<double, double> Getq3Boundary(double nu_Energy_GeV,
                                               double q0_GeV) {
  static const double mass_mu_GeV = 0.1056;
  return {fabs(nu_Energy_GeV -
               sqrt(nu_Energy_GeV * nu_Energy_GeV - mass_mu_GeV * mass_mu_GeV -
                    2 * nu_Energy_GeV * q0_GeV + q0_GeV * q0_GeV)),
          fabs(nu_Energy_GeV +
               sqrt(nu_Energy_GeV * nu_Energy_GeV - mass_mu_GeV * mass_mu_GeV -
                    2 * nu_Energy_GeV * q0_GeV + q0_GeV * q0_GeV))};
}

inline bool IsInSharedq3q0PhaseSpace(double nu_Energy_GeV, double q0_GeV,
                                     double q3_GeV) {

  std::pair<double, double> q3b = Getq3Boundary(nu_Energy_GeV, q0_GeV);
  return ((q3_GeV > q3b.first) && (q3_GeV < q3b.second));
}

inline double GetNueNumuRatioWeight(int nu_pdg, int is_CC, double nu_Energy_GeV,
                                    double q0_GeV, double q3_GeV,
                                    double parameter_value = 1) {

  static double const central_value = 0;
  static double const uncertainty = 1;
  if (abs(nu_pdg) != 12 || !is_CC) {
    return 1;
  }

  return IsInSharedq3q0PhaseSpace(nu_Energy_GeV, q0_GeV, q3_GeV)
             ? 1
             : (central_value + parameter_value * uncertainty);
}
} // namespace nusyst
