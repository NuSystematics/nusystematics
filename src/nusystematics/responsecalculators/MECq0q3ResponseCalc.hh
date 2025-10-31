#ifndef NUSYST_MEC_Q0Q3_RESPONSECALC_HH
#define NUSYST_MEC_Q0Q3_RESPONSECALC_HH

#include <memory>
#include <TH2D.h>
#include <algorithm>
#include <stdexcept>

namespace nusyst {

/// Lightweight helper that evaluates a single 2‑D weight histogram.
class MECq0q3ResponseCalc {
public:



  /// Constructor takes ownership of the histogram (cloned internally).
  MECq0q3ResponseCalc(TH2D* h, double w_min = 0.0, double w_max = 5.0,bool mapIsQ3xQ0 = false);

  /// Central weight with **bilinear interpolation** inside the map.
  double GetCentralWeight(double q0, double q3) const;

  /// Up / Down side variation (ivar = ±1) – symmetric envelope 2 − w_CV.
  double GetVariation(int ivar, double q0, double q3) const;

private:
  std::unique_ptr<TH2D> fHist;   ///< owned, thread‑safe clone of the histogram
  double fWmin, fWmax;      ///< clamp limits
  bool   fMapIsQ3xQ0{false};           ///< true if map is in (q3, q0) coordinates
};

// ---------------------------------------------------------------------------
inline MECq0q3ResponseCalc::MECq0q3ResponseCalc(TH2D* h,
                                                         double wmin,
                                                         double wmax,
                                                         bool mapIsQ3xQ0)
  : fHist(h ? static_cast<TH2D*>(h->Clone()) : nullptr),
    fWmin(wmin), fWmax(wmax)
{
  if (!fHist)
    throw std::runtime_error("Null histogram passed to ResponseCalc");
  fHist->SetDirectory(nullptr);          // detach from any ROOT file
  // Note: kCanRebin is not available in this ROOT version, binning is fixed by design
}

// ---------------------------------------------------------------------------
inline double clamp(double x, double a, double b) {
  return std::max(a, std::min(b, x));
}

// ---------------------------------------------------------------------------
inline double MECq0q3ResponseCalc::GetCentralWeight(double q0, double q3) const
{

  // Histogram axis order = (x = q3, y = q0)
  const double xMin = fHist->GetXaxis()->GetXmin();
  const double xMax = fHist->GetXaxis()->GetXmax();
  const double yMin = fHist->GetYaxis()->GetXmin();
  const double yMax = fHist->GetYaxis()->GetXmax();

  // If kinematics fall outside trained / provided map domain, return neutral weight
  if (q3 < xMin || q3 > xMax || q0 < yMin || q0 > yMax) {
    return 1.0; // do not extrapolate beyond Valencia map coverage
  }

  const double w = fHist->Interpolate(q3, q0); // bilinear inside domain

  return clamp(w, fWmin, fWmax);
}

// ---------------------------------------------------------------------------
inline double MECq0q3ResponseCalc::GetVariation(int ivar,
                                                 double q0, double q3) const
{
  if (ivar == 0) return GetCentralWeight(q0, q3);
  const double w_cv = GetCentralWeight(q0, q3);
  return clamp(2.0 - w_cv, fWmin, fWmax);  // symmetric envelope
}

} // namespace nusyst

#endif // NUSYST_MEC_Q0Q3_RESPONSECALC_HH
