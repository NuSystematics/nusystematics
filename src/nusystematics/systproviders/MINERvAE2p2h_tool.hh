#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/utility/GENIEUtils.hh"

// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/SppChannel.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class MINERvAE2p2h : public nusyst::IGENIESystProvider_tool {

public:
  explicit MINERvAE2p2h(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~MINERvAE2p2h();

private:
  fhicl::ParameterSet tool_options;
  std::pair<double, double> LimitWeights;

  bool ignore_parameter_dependence;

  size_t pidx_E2p2hResponse_nu, pidx_E2p2hResponse_nubar, pidx_E2p2hA_nu,
      pidx_E2p2hB_nu, pidx_E2p2hA_nubar, pidx_E2p2hB_nubar;

  double A_nu_CV, B_nu_CV, A_nubar_CV, B_nubar_CV;

  std::vector<double> A_nu_Variations, B_nu_Variations, A_nubar_Variations,
      B_nubar_Variations;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode;
  double Enu, Q2, weight;
};
