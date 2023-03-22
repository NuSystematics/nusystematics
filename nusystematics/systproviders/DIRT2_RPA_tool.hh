#ifndef nusystematics_SYSTPROVIDERS_DIRT2_RPA_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_DIRT2_RPA_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/DIRT2_RPAEngine_Reweight.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class DIRT2_RPA : public nusyst::IGENIESystProvider_tool {

public:
  explicit DIRT2_RPA(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~DIRT2_RPA();

private:

  fhicl::ParameterSet tool_options;

  size_t pidx_RPA_dummy;
  
  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  // TH: change these!
  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W, Emiss, Emiss_preFSI;
  TVector3 pmiss, pmiss_preFSI;
};

#endif