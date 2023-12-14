#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class NOvAStyleNonResPionNorm : public nusyst::IGENIESystProvider_tool {

  struct channel_param {
    nusyst::NRPiChan_t channel;
    size_t paramidx;
  };
  std::vector<channel_param> ChannelParameterMapping;

public:
  explicit NOvAStyleNonResPionNorm(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~NOvAStyleNonResPionNorm();

private:
  fhicl::ParameterSet tool_options;

  // Dial has uniform response of OneSigmaResponse between WBegin and
  // WTransition, where it linearly reduces to a response of HighWResponse at
  // WEnd to W = infinity
  double WBegin, WEnd, WTransition, OneSigmaResponse, HighWResponse;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, NRPiChannel, NRPiChannel_param, NPi;
  double Enu, Q2, W, weight_m1, weight_p1;
};
