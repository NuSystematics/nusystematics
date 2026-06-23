#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "RwCalculators/GReWeightNuXSecCCQE.h"

#include "TFile.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <string>

class ZExpPCAWeighter : public nusyst::IGENIESystProvider_tool {
public:
  explicit ZExpPCAWeighter(YAML::Node const &);

  //'First' configuration step: tool configuration
  // - takes 'arbitrary' YAML configuration and
  //   produces SystMetaData object which can later be used to configure a
  //   specific set of parameter values to be calculated
  systtools::SystMetaData BuildSystMetaData(YAML::Node const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  // - Reads the preconstructed SystMetaData produced by BuildSystMetaData
  //   to configure an instance of this class to calculate weights
  // - Recieves a copy of the tool_options instance constructed by 
  //   BuildSystMetaData as an argument
  bool SetupResponseCalculator(YAML::Node const &);

  // Used to pass arbitrary YAML options from the tool configuration to the
  //   parameter headers.
  YAML::Node GetExtraToolOptions() { return tool_options; }

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  // Can add as much or as little stateful information here for use when
  // representing this instance as a string.
  std::string AsString() { return "ZExpPCAWeighter"; }

  ~ZExpPCAWeighter(){}

private:
  // arbitrary additional configuration from the tool configuration/parameter
  // headers can be storeds here
  YAML::Node tool_options;

  // The ParamHeaders id of the four free parameters provided by this
  // systprovider
  std::array<size_t, 4> pidx_Params;

  // The configured variations to precalculate for the four parameters
  std::array<double, 4> CVs;
  std::array<std::vector<double>, 4> Variations;

  // Concrete example is more clear with some actual implementation, we will use
  // some GENIE reweight engines.
  std::array<std::vector<std::unique_ptr<genie::rew::GReWeightNuXSecCCQE>>, 4> ReWeightEngines_new;
  
  
  // configurable verbosity as an example of some arbitrary systprovider
  // configuration
  int verbosity_level;
};