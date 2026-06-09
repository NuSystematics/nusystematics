#ifndef nusystematics_SYSTPROVIDERS_FSIReweightMult_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_FSIReweightMult_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/FSIReweightCalculatorMult.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class FSIReweightMult : public nusyst::IGENIESystProvider_tool {

    std::unique_ptr<nusyst::FSIReweightCalculatorMult> fsiReweightCalculatorMult;

public:
  explicit FSIReweightMult(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~FSIReweightMult();

private:

  fhicl::ParameterSet tool_options;

  size_t ResponseParameterIdx;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  // option to save the 2D map from neutrino sample, which can be compared to the input template from hadron sample
  bool save_map;
  TFile* outfile_map;
  TH2D* h_KEini_Ebias;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
};

#endif
