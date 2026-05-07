#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"

#include "Physics/NuclearState/LocalFGM.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Algorithm/AlgConfigPool.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class DUNEDAS2026ExampleReweighter : public nusyst::IGENIESystProvider_tool {

public:
  explicit DUNEDAS2026ExampleReweighter(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString(){ return "DUNEDAS2026ExampleReweighter"; };

  ~DUNEDAS2026ExampleReweighter();

private:

  fhicl::ParameterSet tool_options;
  size_t pidx_DDAS;

};

