#ifndef NUSYST_MEC_Q0Q3_RESPONSECALC_HH
#define NUSYST_MEC_Q0Q3_RESPONSECALC_HH

#include <memory>
#include <TH2D.h>
#include <TH3D.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

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

/// Evaluates a q3 x q0 x Enu map while preserving the legacy interpolation
/// order: interpolate q3/q0 within each Enu plane, clamp each plane weight,
/// then linearly blend the two plane weights in Enu.
class MECq0q3ResponseCalc3D {
public:
  MECq0q3ResponseCalc3D(TH3D* h, double w_min = 0.0,
                        double w_max = 5.0);

  double GetCentralWeight(double q0, double q3, double enu,
                          double enuSnapTol) const;

  void SetUseNearestBin(bool v) { fUseNearestBin = v; }
  void SetEdgeClamp(bool v) { fEdgeClamp = v; }
  void SetOutOfRangeWeight(double w) { fOutOfRangeWeight = w; }

private:
  double GetPlaneWeight(double q0, double q3, int iz) const;

  std::unique_ptr<TH3D> fHist;
  std::vector<double> fEnergyGrid;
  double fWmin, fWmax;
  bool fUseNearestBin{false};
  bool fEdgeClamp{false};
  double fOutOfRangeWeight{1.0};
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

// ---------------------------------------------------------------------------
inline MECq0q3ResponseCalc3D::MECq0q3ResponseCalc3D(
    TH3D* h, double wmin, double wmax)
  : fHist(h ? static_cast<TH3D*>(h->Clone()) : nullptr),
    fWmin(wmin), fWmax(wmax)
{
  if (!fHist)
    throw std::runtime_error("Null histogram passed to 3D ResponseCalc");
  fHist->SetDirectory(nullptr);

  fEnergyGrid.reserve(fHist->GetNbinsZ());
  for (int iz = 1; iz <= fHist->GetNbinsZ(); ++iz)
    fEnergyGrid.push_back(fHist->GetZaxis()->GetBinCenter(iz));
}

// ---------------------------------------------------------------------------
inline double MECq0q3ResponseCalc3D::GetPlaneWeight(
    double q0, double q3, int iz) const
{
  const double xMin = fHist->GetXaxis()->GetXmin();
  const double xMax = fHist->GetXaxis()->GetXmax();
  const double yMin = fHist->GetYaxis()->GetXmin();
  const double yMax = fHist->GetYaxis()->GetXmax();
  const bool inX = (q3 >= xMin && q3 <= xMax);
  const bool inY = (q0 >= yMin && q0 <= yMax);

  double w = 1.0;
  if (fUseNearestBin) {
    int ix = fHist->GetXaxis()->FindBin(q3);
    int iy = fHist->GetYaxis()->FindBin(q0);
    if (!inX || !inY) {
      if (!fEdgeClamp) return fOutOfRangeWeight;
      ix = std::clamp(ix, 1, fHist->GetNbinsX());
      iy = std::clamp(iy, 1, fHist->GetNbinsY());
    }
    w = fHist->GetBinContent(ix, iy, iz);
  } else {
    if (!inX || !inY) return fOutOfRangeWeight;
    double enuCenter = fHist->GetZaxis()->GetBinCenter(iz);
    // ROOT's TH3::Interpolate requires a second valid bin on the side of
    // each coordinate.  At the final z-bin center it selects the overflow
    // bin even though the interpolation fraction is zero.  Nudge the query
    // toward the preceding plane so the final stored plane is evaluated.
    if (iz == fHist->GetNbinsZ())
      enuCenter = std::nextafter(enuCenter,
                                 -std::numeric_limits<double>::infinity());
    w = fHist->Interpolate(q3, q0, enuCenter);
  }
  return clamp(w, fWmin, fWmax);
}

// ---------------------------------------------------------------------------
inline double MECq0q3ResponseCalc3D::GetCentralWeight(
    double q0, double q3, double enu, double enuSnapTol) const
{
  auto itHi = std::lower_bound(fEnergyGrid.begin(), fEnergyGrid.end(), enu);
  size_t ih = (itHi == fEnergyGrid.end()) ? fEnergyGrid.size() - 1
                                          : std::distance(fEnergyGrid.begin(), itHi);
  size_t il = (ih == 0) ? 0 : ih - 1;

  const double elo = fEnergyGrid[il];
  const double ehi = fEnergyGrid[ih];
  double t = 0.0;
  if (std::fabs(enu - elo) <= enuSnapTol) {
    ih = il;
  } else if (std::fabs(enu - ehi) <= enuSnapTol) {
    il = ih;
  } else {
    t = (ih == il || ehi <= elo)
        ? 0.0
        : std::clamp((enu - elo) / (ehi - elo), 0.0, 1.0);
  }

  const double wLo = GetPlaneWeight(q0, q3, static_cast<int>(il) + 1);
  const double wHi = GetPlaneWeight(q0, q3, static_cast<int>(ih) + 1);
  return clamp((1.0 - t) * wLo + t * wHi, fWmin, fWmax);
}

} // namespace nusyst

#endif // NUSYST_MEC_Q0Q3_RESPONSECALC_HH
