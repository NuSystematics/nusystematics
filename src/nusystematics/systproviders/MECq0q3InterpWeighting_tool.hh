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

  // Ramp guard: unity below Min; linear blend to full weight by Max
  // If Max == Min, behaves like a hard cutoff at Min
  double fQ0GuardMin{0.0};  // GeV
  double fQ0GuardMax{0.0};  // GeV

  // Energy-guard window and snapping
  // Only energies within [fEnuMin, fEnuMax] are reweighted.
  // If |Enu - grid_point| <= fEnuSnapTol, use the exact map (no blending).
  double fEnuMin{0.4};      // GeV
  double fEnuMax{2.5};      // GeV
  double fEnuSnapTol{5e-3}; // GeV

  std::unordered_map<Topo,
      std::vector<std::unique_ptr<MECq0q3ResponseCalc>>> fCalcs;
};

} // namespace nusyst

#endif // NUSYST_MEC_Q0Q3_INTERPWEIGHTING_TOOL_HH
