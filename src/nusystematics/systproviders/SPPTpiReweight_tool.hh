#ifndef nusystematics_SYSTPROVIDERS_SPPTpiReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_SPPTpiReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/SPPTpiReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class SPPTpiReweight : public nusyst::IGENIESystProvider_tool {

public:

  NEW_SYSTTOOLS_EXCEPT(invalid_engine_state);

  explicit SPPTpiReweight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~SPPTpiReweight();

private:

  fhicl::ParameterSet tool_options;

  size_t pidx_SPPTpiCVCorrection;
  size_t pidx_SPPTpiCorrectionRW;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
};

#endif
