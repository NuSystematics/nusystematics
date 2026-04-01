#ifndef nusystematics_SYSTPROVIDERS_CCQETemplateReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_CCQETemplateReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/CCQETemplateReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class CCQETemplateReweight : public nusyst::IGENIESystProvider_tool {

  enum RWMode {
    NONE = 0,
    q3q0 = 1,
    PCTheta = 2,
    PSTheta = 3,
    PTheta = 4
  };

  std::unique_ptr<nusyst::CCQETemplateReweightCalculator> ccqeTemplateReweightCalculator;

public:
  explicit CCQETemplateReweight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~CCQETemplateReweight();

private:

  fhicl::ParameterSet tool_options;
  RWMode rwMode;

  void InitValidTree();
  
  int GetQ0BinIndex(double q0_value) const;

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;

  // q0 bins (GeV)
  // note that 10 GeV cutoff is ~dummy, since there is a Enu cutoff of 2 GeV
  static constexpr std::size_t kNumQ0Bins = 5;
  static constexpr std::array<double, kNumQ0Bins + 1> q0_bin_boundaries = {0.0, 0.05, 0.10, 0.20, 0.7, 10.0};
  std::array<size_t, kNumQ0Bins> ResponseParameterIdx_q0bin;
  std::array<std::string, kNumQ0Bins> q0bin_params_for_ccqe_sf_names;
};

#endif
