#ifndef nusystematics_SYSTPROVIDERS_CCQERPAReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_CCQERPAReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/CCQERPAReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class CCQERPAReweight : public nusyst::IGENIESystProvider_tool {

  enum RWMode {
    NONE = 0,
    q3q0 = 1,
    PCTheta = 2,
    PSTheta = 3,
    PTheta = 4
  };

  std::unique_ptr<nusyst::CCQERPAReweightCalculator> ccqeRPAReweightCalculator;

public:
  explicit CCQERPAReweight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~CCQERPAReweight();

private:

  fhicl::ParameterSet tool_options;
  RWMode rwMode;

  size_t ResponseParameterIdx;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
};

#endif
