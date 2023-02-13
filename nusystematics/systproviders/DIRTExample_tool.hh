#ifndef nusystematics_SYSTPROVIDERS_DIRTEXAMPLE_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_DIRTEXAMPLE_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/DIRTExample_Reweight.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class DIRTExample : public nusyst::IGENIESystProvider_tool {

public:
  explicit DIRTExample(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~DIRTExample();

private:

  fhicl::ParameterSet tool_options;

  size_t pidx_AltModelA;
  size_t pidx_AltModelB;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
};

#endif
