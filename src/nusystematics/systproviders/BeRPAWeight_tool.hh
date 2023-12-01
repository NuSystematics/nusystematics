#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class BeRPAWeight : public nusyst::IGENIESystProvider_tool {

public:
  explicit BeRPAWeight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~BeRPAWeight();

private:
  fhicl::ParameterSet tool_options;

  bool ApplyCV;

  bool ignore_parameter_dependence;

  size_t pidx_BeRPA_Response, pidx_BeRPA_A, pidx_BeRPA_B, pidx_BeRPA_D,
      pidx_BeRPA_E;
  double ACV, BCV, DCV, ECV;

  std::vector<double> AVariations, BVariations, DVariations, EVariations;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode;
  double Enu, Q2, weight;
};
