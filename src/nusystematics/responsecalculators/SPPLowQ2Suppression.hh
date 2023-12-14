#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

///
///
///\note From Phys. Rev. D 91, 012005 (2015) Fig. 9
inline double GetMINOSSPPLowQ2SuppressionWeight(int mode_simb, int is_CC,
                                                double Q2_GeV,
                                                double parameter_value = 1) {

  static double const A = 1.01;
  static double const Q2_Max = 0.7;
  static double const Q0 = 0.156;
  static double const central_value = 0;
  static double const uncertainty = 1;

  if ((mode_simb != e2i(simb_mode_copy::kRes)) || !is_CC) {
    return 1;
  }

  if ((Q2_GeV > Q2_Max) || (Q2_GeV < 0)) {
    return 1;
  }

  return (central_value + parameter_value * uncertainty) *
         (A / (1 + exp(1 - (sqrt(Q2_GeV) / Q0))));
}

///
///
///\note See P. Stowells thesis or upcoming MINERvA pion tune publication.
inline double GetMINERvASPPLowQ2SuppressionWeight(int mode_simb, int is_CC,
                                                  double Q2_GeV,
                                                  double parameter_value = 1) {

  // Fit parameters for FrInel+Low-Q2 tune were
  // MA_RES = 0.93 +/- 0.05
  // NormRes = 116 +/- 7
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 120 +/- 32
  // ThetaPi = 1.0 (at limit)
  // FrInel = 132 +/- 27
  // R1 = 0.37 +/- 0.09
  // R2 = 0.60 +/- 0.16
  //
  // Fit parameters for FrAbs+Low-Q2 tune
  // MA_RES = 0.92 +/- 0.02
  // NormRes = 116 +/- 3
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 99 +/- 31
  // ThetaPi = 1.0 (at limit)
  // FrAbs = 48 +/- 21
  // R1 = 0.32 +/- 0.06
  // R2 = 0.5 (limit)
  static double const Q2_Max = 0.7;
  static double const Q2_t1 = 0;
  static double const Q2_t2 = 0.35;
  static double const R1 = 0.3;
  static double const R2 = 0.6;

  if ((mode_simb != e2i(simb_mode_copy::kRes)) || !is_CC) {
    return 1;
  }

  if ((Q2_GeV > Q2_Max) || (Q2_GeV < 0)) {
    return 1;
  }

  double RQ2 = (R2 * ((Q2_GeV - Q2_t1) * (Q2_GeV - Q2_Max)) /
         ((Q2_t2 - Q2_t1) * (Q2_t2 - Q2_Max))) +
        (((Q2_GeV - Q2_t1) * (Q2_GeV - Q2_t2)) /
         ((Q2_Max - Q2_t1) * (Q2_Max - Q2_t2)));
  return 1 - parameter_value * ((1 - R1) * pow((1 - RQ2), 2));
}
} // namespace nusyst
