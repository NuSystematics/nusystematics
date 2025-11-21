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

  /// Up / Down side variation (ivar = ±1) – symmetric envelope 2 − w_CV.
  double GetVariation(int ivar, double q0, double q3) const;

  void   SetUseNearestBin(bool v) { fUseNearestBin = v; }
  void   SetEdgeClamp(bool v)     { fEdgeClamp = v; }
  void   SetOutOfRangeWeight(double w) { fOutOfRangeWeight = w; }

private:
  std::unique_ptr<TH2D> fHist;   ///< owned, thread‑safe clone of the histogram
  double fWmin, fWmax;      ///< clamp limits
  bool   fMapIsQ3xQ0{false};           ///< true if map is in (q3, q0) coordinates
  bool   fUseNearestBin{false};   ///< default off → legacy bilinear
  bool   fEdgeClamp{false};       ///< default off → legacy out-of-range=1
  double fOutOfRangeWeight{1.0};  ///< weight to return when outside histogram bounds (default 1.0 for backward compatibility, set to 0.0 to suppress)
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

  auto in_x = (q3 >= xMin && q3 <= xMax);
  auto in_y = (q0 >= yMin && q0 <= yMax);

  double w = 1.0;
  if (fUseNearestBin) {
    // --- piecewise-constant: nearest-bin content ---
    int ix = fHist->GetXaxis()->FindBin(q3);
    int iy = fHist->GetYaxis()->FindBin(q0);
    // clamp to valid bins (1..Nbins) if requested, else treat OOR as configurable
    if (!in_x || !in_y) {
      if (!fEdgeClamp) return fOutOfRangeWeight;
      ix = std::clamp(ix, 1, fHist->GetNbinsX());
      iy = std::clamp(iy, 1, fHist->GetNbinsY());
    }
    w = fHist->GetBinContent(ix, iy);
  } else {
    // --- legacy: bilinear interpolation inside domain; configurable weight outside ---
    if (!in_x || !in_y) return fOutOfRangeWeight;
    w = fHist->Interpolate(q3, q0);
  }
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
