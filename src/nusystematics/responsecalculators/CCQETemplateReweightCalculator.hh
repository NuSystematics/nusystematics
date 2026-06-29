#ifndef nusystematics_RESPONSE_CALCULATORS_CCQETemplateReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCQETemplateReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"
#include "systematicstools/interpreters/PolyResponse.hh"
#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/string_parsers.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

#include <map>
#include <memory>
#include <string>

NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_Template_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_Template_FILEPATH);
NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_Template_PDG);

namespace nusyst {

  // Return the string flavor label for a neutrino PDG code, or "" if
  // the PDG is not one of the four supported flavors (±12, ±14)
  inline std::string pdg_to_flavor_label(int pdg) {
    switch (pdg) {
      case  12: return "nue";
      case -12: return "nuebar";
      case  14: return "numu";
      case -14: return "numubar";
      default:  return "";
    }
  }

  class CCQETemplateReweightCalculator {

  protected:

    // Baseline and target absolute cross-section histograms per PDG code
    // (outer key) and energy range (inner key = 0 for Enu < ENuBoundary)
    //
    // Histogram naming convention in the combined ROOT file:
    //   {flavor_label}_{tune}   e.g. "numu_AR23_20i_01_000"
    //
    // enu_range 0 → LowE (Enu < ENuBoundary): loaded from file.
    // enu_range 1 → HighE (Enu ≥ ENuBoundary): returns weight 1.0, not loaded.
    std::map<int, std::map<int, std::unique_ptr<TH3D>>> map_PDG_ENuRange_to_BaselineXSec;
    std::map<int, std::map<int, std::unique_ptr<TH3D>>> map_PDG_ENuRange_to_TargetXSec;

    // Per-PDG, per-enu_range axis clamping limits (bin-center coordinates)
    // Stored to avoid repeated calls to GetBinCenter during event processing
    // Both baseline and target should share identical binning
    std::map<int, std::map<int, double>> pdg_x_FirstBinCenter, pdg_x_LastBinCenter;
    std::map<int, std::map<int, double>> pdg_y_FirstBinCenter, pdg_y_LastBinCenter;
    std::map<int, std::map<int, double>> pdg_z_FirstBinCenter, pdg_z_LastBinCenter;

    double ENuBoundary{2.0};
    std::string baseline_tune{};
    std::string target_tune{};

  public:

    explicit CCQETemplateReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~CCQETemplateReweightCalculator() {}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    // Compute the reweight for a single event
    //
    // Enu_GeV         : neutrino energy 
    // bin_kin         : kinematic coordinates {Y, Z} matching the histogram axis layout
    //                   (ex. Y=q0, Z=q3 for q3q0 mode)
    // parameter_value : 0 → no change, 1 → reweight
    // pdg_nu          : neutrino PDG code -- must be one of ±12, ±14

    double GetTemplateReweight(double Enu_GeV, std::array<double, 2> bin_kin,
                               double parameter_value, int pdg_nu);

    std::string GetCalculatorName() const { return "CCQETemplateReweightCalculator"; }
  };


  // ---------------------------------------------------------------------------
  // GetTemplateReweight
  // ---------------------------------------------------------------------------

  inline double CCQETemplateReweightCalculator::GetTemplateReweight(
      double Enu_GeV, std::array<double, 2> bin_kin,
      double parameter_value, int pdg_nu) {

    // ------------------------------------------------------------------
    // Validate PDG
    // ------------------------------------------------------------------
    std::string flavor = pdg_to_flavor_label(pdg_nu);
    if (flavor.empty()) {
      throw invalid_CCQE_Template_PDG()
          << "[ERROR]: Neutrino PDG " << pdg_nu
          << " is not one of the four supported flavors (±12, ±14).";
    }
    if (map_PDG_ENuRange_to_BaselineXSec.find(pdg_nu) ==
        map_PDG_ENuRange_to_BaselineXSec.end()) {
      throw invalid_CCQE_Template_PDG()
          << "[ERROR]: No template loaded for flavor " << flavor
          << " (PDG " << pdg_nu << "). "
          << "Check that make_reweight_templates.py was run for this flavor "
          << "and that combine_flavor_templates.py included it.";
    }

    // ------------------------------------------------------------------
    // Energy-range gate
    // ------------------------------------------------------------------
    if (Enu_GeV >= ENuBoundary) {
      return 1.0;  // No reweighting above ENuBoundary
    }

    // ------------------------------------------------------------------
    // Clamp coordinates to the innermost bin centers to avoid extrapolation
    // ------------------------------------------------------------------
    static const double epsil = 1E-6;
    const int er = 0;  // enu_range for LowE

    double Enu_c = std::max(Enu_GeV,    pdg_x_FirstBinCenter.at(pdg_nu).at(er) + epsil);
    Enu_c        = std::min(Enu_c,      pdg_x_LastBinCenter .at(pdg_nu).at(er) - epsil);

    double Y_c   = std::max(bin_kin[0], pdg_y_FirstBinCenter.at(pdg_nu).at(er) + epsil);
    Y_c          = std::min(Y_c,        pdg_y_LastBinCenter .at(pdg_nu).at(er) - epsil);

    double Z_c   = std::max(bin_kin[1], pdg_z_FirstBinCenter.at(pdg_nu).at(er) + epsil);
    Z_c          = std::min(Z_c,        pdg_z_LastBinCenter .at(pdg_nu).at(er) - epsil);

    // ------------------------------------------------------------------
    // Compute ratio = sigma_target / sigma_baseline 
    // ------------------------------------------------------------------
    double sigma_baseline =
        map_PDG_ENuRange_to_BaselineXSec.at(pdg_nu).at(er)->Interpolate(Enu_c, Y_c, Z_c);
    double sigma_target   =
        map_PDG_ENuRange_to_TargetXSec  .at(pdg_nu).at(er)->Interpolate(Enu_c, Y_c, Z_c);

    if (sigma_baseline <= 0.0) {
      // Baseline has no cross-section here 
      // return weight = 1 (no change)
      return 1.0;
    }

    double ratio = sigma_target / sigma_baseline;

    // Guard against NaN from interpolation at sparse template regions.
    if (ratio != ratio) {
      std::cout << "[CCQETemplateReweightCalculator] WARNING: NaN ratio for "
                << flavor << " (PDG=" << pdg_nu << ") at"
                << " Enu=" << Enu_GeV << " Y=" << bin_kin[0] << " Z=" << bin_kin[1]
                << "; returning weight=1." << std::endl;
      return 1.0;
    }

    // Linear interpolation between no-change (param=0) and full reweight (param=1).
    //   weight = 1 + parameter_value * (ratio - 1)
    return 1.0 + parameter_value * (ratio - 1.0);
  }


  // ---------------------------------------------------------------------------
  // LoadInputHistograms
  // ---------------------------------------------------------------------------

  inline void CCQETemplateReweightCalculator::LoadInputHistograms(
      fhicl::ParameterSet const &ps) {

    // -----------------------------------------------------------------------
    // Read configuration
    // -----------------------------------------------------------------------
    const std::string input_file_cfg = ps.get<std::string>("input_file");
    ENuBoundary   = ps.get<double>("ENuBoundary");
    baseline_tune = ps.get<std::string>("baseline_tune");
    target_tune   = ps.get<std::string>("target_tune");

    printf("[CCQETemplateReweightCalculator] ENuBoundary  = %1.2f GeV\n", ENuBoundary);
    printf("[CCQETemplateReweightCalculator] baseline_tune = %s\n", baseline_tune.c_str());
    printf("[CCQETemplateReweightCalculator] target_tune   = %s\n", target_tune.c_str());

    if (baseline_tune == target_tune) {
      std::cout << "[CCQETemplateReweightCalculator] WARNING: baseline_tune == target_tune ("
                << baseline_tune << "); all weights will be 1." << std::endl;
    }

    // Resolve relative paths under ${nusystematics_ROOT}/data/
    std::string input_file = systtools::expand_env_vars(input_file_cfg);
    printf("[CCQETemplateReweightCalculator] input_file    = %s\n", input_file.c_str());

    // -----------------------------------------------------------------------
    // Load baseline and target histograms for all four supported flavors
    // -----------------------------------------------------------------------
    const int er = 0;  // enu_range index for LowE
    const std::vector<int> supported_pdgs = {12, -12, 14, -14};
    int n_loaded = 0;

    for (const int pdg : supported_pdgs) {
      const std::string flavor = pdg_to_flavor_label(pdg);
      const std::string baseline_hist_name = flavor + "_" + baseline_tune;
      const std::string target_hist_name   = flavor + "_" + target_tune;

      printf("[CCQETemplateReweightCalculator] Loading histograms for %s (PDG %d):\n",
             flavor.c_str(), pdg);
      printf("  baseline: %s\n", baseline_hist_name.c_str());
      printf("  target  : %s\n", target_hist_name.c_str());

      // -- Baseline --
      TH3D* h_baseline = GetHistogram<TH3D>(input_file, baseline_hist_name);
      if (!h_baseline) {
        throw invalid_CCQE_Template_FILEPATH()
            << "[ERROR]: Failed to load baseline TH3D '" << baseline_hist_name
            << "' from '" << input_file << "' for flavor " << flavor
            << " (PDG " << pdg << "). "
            << "Run make_reweight_templates.py for this flavor with baseline_tune='"
            << baseline_tune << "' then rerun combine_flavor_templates.py.";
      }

      // -- Target --
      TH3D* h_target = GetHistogram<TH3D>(input_file, target_hist_name);
      if (!h_target) {
        throw invalid_CCQE_Template_FILEPATH()
            << "[ERROR]: Failed to load target TH3D '" << target_hist_name
            << "' from '" << input_file << "' for flavor " << flavor
            << " (PDG " << pdg << "). "
            << "Run make_reweight_templates.py for this flavor with target_tune='"
            << target_tune << "' then rerun combine_flavor_templates.py.";
      }

      // Transfer to maps
      map_PDG_ENuRange_to_BaselineXSec[pdg][er] = std::unique_ptr<TH3D>(h_baseline);
      map_PDG_ENuRange_to_TargetXSec  [pdg][er] = std::unique_ptr<TH3D>(h_target);

      // Cache axis bin-center limits for clamping
      const TAxis* xax = h_baseline->GetXaxis();
      const TAxis* yax = h_baseline->GetYaxis();
      const TAxis* zax = h_baseline->GetZaxis();

      pdg_x_FirstBinCenter[pdg][er] = xax->GetBinCenter(1);
      pdg_x_LastBinCenter [pdg][er] = xax->GetBinCenter(xax->GetNbins());
      pdg_y_FirstBinCenter[pdg][er] = yax->GetBinCenter(1);
      pdg_y_LastBinCenter [pdg][er] = yax->GetBinCenter(yax->GetNbins());
      pdg_z_FirstBinCenter[pdg][er] = zax->GetBinCenter(1);
      pdg_z_LastBinCenter [pdg][er] = zax->GetBinCenter(zax->GetNbins());

      printf("[CCQETemplateReweightCalculator]   %s x=[%.3f, %.3f] y=[%.3f, %.3f] z=[%.3f, %.3f]\n",
             flavor.c_str(),
             pdg_x_FirstBinCenter[pdg][er], pdg_x_LastBinCenter[pdg][er],
             pdg_y_FirstBinCenter[pdg][er], pdg_y_LastBinCenter[pdg][er],
             pdg_z_FirstBinCenter[pdg][er], pdg_z_LastBinCenter[pdg][er]);
      ++n_loaded;
    }

    printf("[CCQETemplateReweightCalculator] Loaded baseline+target histograms for %d flavor(s).\n",
           n_loaded);
  }

} // namespace nusyst

#endif
