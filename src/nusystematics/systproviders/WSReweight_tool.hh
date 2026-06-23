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

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  bool estimate_emiss;


  // TH: change these!
  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W, Emiss, Emiss_preFSI, KF_tree, radius, ref_prob_density, new_prob_density;
  TVector3 pmiss, pmiss_preFSI;
};

