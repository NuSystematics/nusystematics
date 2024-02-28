#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "RwCalculators/GReWeightNuXSecCCRES.h"

#include "TFile.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <string>

//This file helps correctly transpose the appropriate .fcl file (format found in ResIso.fcl) into the usable .fcl file

class ResIso : public nusyst::IGENIESystProvider_tool {
public:
  explicit ResIso(fhicl::ParameterSet const &);

  //'First' configuration step: tool configuration
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  // Used to pass arbitrary FHiCL options from the tool configuration to the parameter headers.
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  // Can add as much or as little stateful information here for use when representing this instance as a string.
  std::string AsString() { return "ResIso"; }

  ~ResIso(){}

private:
  // arbitrary additional configuration from the tool configuration/parameter headers can be storeds here
  fhicl::ParameterSet tool_options;

  // The ParamHeaders id of the four free parameters provided by this systprovider
  std::vector<size_t> pidx_Params;

  // The configured variations to precalculate for the four parameters
  std::vector<double> CVs;
  std::vector<std::vector<double>> Variations;

  // Concrete example is more clear with some actual implementation, we will use some GENIE reweight engines.
  std::vector<std::vector<genie::rew::GReWeightNuXSecCCRES>> ReWeightEngines;

  struct DialInfo {
      std::string prettyname;
      int rescode;
      genie::rew::GSyst_t geniedial;
  };
    std::vector<DialInfo> dial_infos;

  // configurable verbosity as an example of some arbitrary systprovider configuration
  int verbosity_level;
};
