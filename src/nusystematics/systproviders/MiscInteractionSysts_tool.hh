#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class MiscInteractionSysts : public nusyst::IGENIESystProvider_tool {

public:
  explicit MiscInteractionSysts(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~MiscInteractionSysts();

private:
  fhicl::ParameterSet tool_options;

  size_t pidx_C12ToAr40_2p2hScaling_nu;
  size_t pidx_C12ToAr40_2p2hScaling_nubar;
  size_t pidx_nuenuebar_xsec_ratio;
  size_t pidx_nuenumu_xsec_ratio;
  size_t pidx_SPPLowQ2Suppression;

  std::vector<double>
  GetWeights_C12ToAr40_2p2hScaling(genie::EventRecord const &,
                                   std::vector<double> const &);
  std::vector<double>
  GetWeights_nuenuebar_xsec_ratio(genie::EventRecord const &,
                                  std::vector<double> const &);
  std::vector<double>
  GetWeights_nuenumu_xsec_ratio(genie::EventRecord const &,
                                std::vector<double> const &);
  std::vector<double>
  GetWeights_SPPLowQ2Suppression(genie::EventRecord const &,
                                 std::vector<double> const &);

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode;
  double Enu, Q2, W;
};
