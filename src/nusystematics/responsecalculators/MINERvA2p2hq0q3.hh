#pragma once

#include <array>

#include "nusystematics/utility/enumclass2int.hh"

namespace nusyst {
enum class PI {
  kNorm = 0,
  kMeanQ0 = 1,
  kMeanQ3 = 2,
  kSigmaQ0 = 3,
  kSigmaQ3 = 4,
  kCorrelation = 5
};

enum class MECTweak_t { kCV = 0, kNNOnly, kNPOnly, kQEOnly };

double Gaussian2D(double q0, double q3, std::array<double, 6> const &GParams) {
  // Unpack stats values from input array and combine w class members

  double z =
      (q0 - GParams[e2i(PI::kMeanQ0)]) * (q0 - GParams[e2i(PI::kMeanQ0)]) /
          GParams[e2i(PI::kSigmaQ0)] / GParams[e2i(PI::kSigmaQ0)] +
      (q3 - GParams[e2i(PI::kMeanQ3)]) * (q3 - GParams[e2i(PI::kMeanQ3)]) /
          GParams[e2i(PI::kSigmaQ3)] / GParams[e2i(PI::kSigmaQ3)] -
      2 * GParams[e2i(PI::kCorrelation)] * (q0 - GParams[e2i(PI::kMeanQ0)]) *
          (q3 - GParams[e2i(PI::kMeanQ3)]) /
          (GParams[e2i(PI::kSigmaQ0)] * GParams[e2i(PI::kSigmaQ3)]);

  double ret =
      GParams[e2i(PI::kNorm)] * exp(-0.5 * z /
                                    (1 - GParams[e2i(PI::kCorrelation)] *
                                             GParams[e2i(PI::kCorrelation)]));

  // Need to add 1 to the results
  return ret;
}
} // namespace nusyst
