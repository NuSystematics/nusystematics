#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"

#include "nusystematics/responsecalculators/weightAMUvariations.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class AMUVariation : public nusyst::IGENIESystProvider_tool {

public:
  explicit AMUVariation(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~AMUVariation();

private:

  std::unique_ptr<nusyst::weightAMUvariations> AMUCalculator;

  size_t ResponseParameterIdx;
  fhicl::ParameterSet tool_options;

  void InitValidTree();



};

