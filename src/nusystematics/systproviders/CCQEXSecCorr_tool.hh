#ifndef nusystematics_SYSTPROVIDERS_CCQEXSecCorr_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_CCQEXSecCorr_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/responsecalculators/MECq0q3ResponseCalc.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class CCQEXSecCorr : public nusyst::IGENIESystProvider_tool {

public:
  explicit CCQEXSecCorr(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~CCQEXSecCorr();

private:
  fhicl::ParameterSet tool_options;

  size_t ResponseParameterIdx;

  // Energy grid (GeV) — one calculator per energy point per flavor
  std::vector<double> fEgrid;

  // Calculators keyed by neutrino PDG code (12, -12, 14, -14).
  // Each vector is parallel to fEgrid.
  std::unordered_map<int,
      std::vector<std::unique_ptr<nusyst::MECq0q3ResponseCalc>>> fCalcs;

  // Weight clamp limits
  double fWmin{0.0};
  double fWmax{30.0};

  // Energy guard window and snapping
  double fEnuMin{0.20};
  double fEnuMax{3.00};
  double fEnuSnapTol{5e-3}; // GeV

  // Validation tree
  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W, weight;

  void InitValidTree();
};

#endif
