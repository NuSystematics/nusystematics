/*******************************************************************************
 *   MECq0q3InterpWeighting_tool.hh
 ******************************************************************************/
#ifndef NUSYST_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH
#define NUSYST_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/responsecalculators/MECq0q3ResponseCalc.hh"
#include <unordered_map>
#include <memory>
#include <vector>
#include <limits>

namespace nusyst {

class MECq0q3InterpWeighting : public IGENIESystProvider_tool {
public:
  explicit MECq0q3InterpWeighting(const fhicl::ParameterSet& pset);

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t) override;

  bool SetupResponseCalculator(fhicl::ParameterSet const &) override;

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const& ev) override;

private:
  fhicl::ParameterSet tool_options;

  enum class Topo { np = 0, nn = 1, unknown = 2 };

  static Topo  ClassifyEvent(genie::EventRecord const&);
  static void  ComputeQ0Q3(genie::EventRecord const&, double& q0, double& q3,
                           double& Enu);

  // ---- configuration/state ----
  std::vector<double> fEgrid; ///< GeV


  

  // Clamp range
  double fWmin{0.0};
  double fWmax{5.0};


  // NEW: q0 apply window (gate). Weighting applies only if
  //      fQ0ApplyMin < q0 < fQ0ApplyMax. Defaults: [0, +inf)
  double fQ0ApplyMin{0.0};                                   // GeV
  double fQ0ApplyMax{std::numeric_limits<double>::infinity()}; // GeV

  // NEW: q3 apply window (gate). Weighting applies only if
  //      fQ3ApplyMin < q3 < fQ3ApplyMax. Defaults: [0, +inf)
  double fQ3ApplyMin{0.0};                                   // GeV
  double fQ3ApplyMax{std::numeric_limits<double>::infinity()}; // GeV

  // NEW: Additional Q0 selection range for fine-grained control
  // If set, reweighting is applied ONLY in this range; outside it, weight=1
  // This is independent of Q0ApplyMin/Q0ApplyMax above
  // Defaults: disabled (both 0.0 means no selection)
  double fQ0SelectMin{0.0};  // GeV
  double fQ0SelectMax{0.0};  // GeV

  // NEW: Q0 binning for multiple dials
  // If Q0Bins is provided, each bin gets its own dial
  // Format: [edge0, edge1, edge2, ...] creates bins [edge0,edge1), [edge1,edge2), ...
  // Empty vector means single dial mode (backward compatible)
  std::vector<double> fQ0Bins;  // GeV bin edges
  
  // Energy-guard window and snapping
  // Only energies within [fEnuMin, fEnuMax] are reweighted.
  // If |Enu - grid_point| <= fEnuSnapTol, use the exact map (no blending).
  double fEnuMin{0.4};      // GeV
  double fEnuMax{2.5};      // GeV
  double fEnuSnapTol{5e-3}; // GeV

  std::unordered_map<Topo,
      std::vector<std::unique_ptr<MECq0q3ResponseCalc>>> fCalcs;
  
  // Helper to determine which q0 bin (dial index) an event belongs to
  // Returns -1 if outside all bins
  int GetQ0BinIndex(double q0) const;
};

} // namespace nusyst

#endif // NUSYST_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH
