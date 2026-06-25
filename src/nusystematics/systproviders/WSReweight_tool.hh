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

class WSReweight : public nusyst::IGENIESystProvider_tool {

public:
  explicit WSReweight(YAML::Node const &);

  bool SetupResponseCalculator(YAML::Node const &);

  YAML::Node GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(YAML::Node const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~WSReweight();

private:

  YAML::Node tool_options;

  size_t pidx_nucleus_radius;
  size_t pidx_surface_thickness;

  bool estimate_emiss;

};

