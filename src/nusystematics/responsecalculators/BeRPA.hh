#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

namespace BeRPA_consts {
/// C. Wilkinson's fitted nominals to C12 calculation
/// A = 0.59 +/- 20%
/// B = 1.05 +/- 20%
/// D = 1.13 +/- 15%
/// E = 0.88 +/- 40%
/// U = 1.2

static size_t const kC12 = 0;
static size_t const kAr40 = 1;
static double const A_central[] = {0.59, 0};
static double const A_frac_uncert[] = {0.2, 0};
static double const B_central[] = {1.05, 0};
static double const B_frac_uncert[] = {0.2, 0};
static double const D_central[] = {1.13, 0};
static double const D_frac_uncert[] = {0.15, 0};
static double const E_central[] = {0.88, 0};
static double const E_frac_uncert[] = {0.4, 0};
static double const U_central[] = {1.2, 0};

} // namespace BeRPA_consts

inline double
EvalBeRPA(double Q2_GeV2,
          double A = BeRPA_consts::A_central[BeRPA_consts::kC12],
          double B = BeRPA_consts::B_central[BeRPA_consts::kC12],
          double D = BeRPA_consts::D_central[BeRPA_consts::kC12],
          double E = BeRPA_consts::E_central[BeRPA_consts::kC12],
          double U = BeRPA_consts::U_central[BeRPA_consts::kC12]) {

  // Q2 transition; if Q2 less than U -> polynominal, otherwise exponential
  // decay
  if (Q2_GeV2 < U) {
    // xprime as prescribed by C. Wilkinson
    const double xprime = Q2_GeV2 / U;
    const double one_minus_xprime = 1. - xprime;
    const double C = D + U * E * (D - 1) / 3.;
    return A * one_minus_xprime * one_minus_xprime * one_minus_xprime +
           3 * B * one_minus_xprime * one_minus_xprime * xprime +
           3 * C * one_minus_xprime * xprime * xprime +
           D * xprime * xprime * xprime;
  } else {
    return 1 + (D - 1) * std::exp(-E * (Q2_GeV2 - U));
  }
}

///\brief An interface to BeRPA that accepts dial 'tweaks' away from nominal
/// instead of absolute values.
///
/// For the nominal parameter value, X_tweak = 0, for +1 sigma X_tweak = 1,
/// etc...
inline double GetBeRPAWeight(int mode_simb, bool is_CC, double Q2_GeV2,
                             double A_tweak = 0, double B_tweak = 0,
                             double D_tweak = 0, double E_tweak = 0,
                             size_t tgtidx = BeRPA_consts::kC12) {

  if ((mode_simb != e2i(simb_mode_copy::kQE)) || !is_CC) {
    return 1;
  }

  return EvalBeRPA(Q2_GeV2,
                   BeRPA_consts::A_central[tgtidx] +
                       A_tweak * BeRPA_consts::A_central[tgtidx] *
                           BeRPA_consts::A_frac_uncert[tgtidx],
                   BeRPA_consts::B_central[tgtidx] +
                       B_tweak * BeRPA_consts::B_central[tgtidx] *
                           BeRPA_consts::B_frac_uncert[tgtidx],
                   BeRPA_consts::D_central[tgtidx] +
                       D_tweak * BeRPA_consts::D_central[tgtidx] *
                           BeRPA_consts::D_frac_uncert[tgtidx],
                   BeRPA_consts::E_central[tgtidx] +
                       E_tweak * BeRPA_consts::E_central[tgtidx] *
                           BeRPA_consts::E_frac_uncert[tgtidx],
                   BeRPA_consts::U_central[tgtidx]);
}

} // namespace nusyst
