// File: ValenciaExc2p2hReweighter_tool.hh
// Author: Zihao Lin zlin22@ur.rochester.edu

#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"
#include "nusystematics/utility/BDTReweighter_json.hh"
#include "Physics/NuclearState/LocalFGM.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "TFile.h"
#include "TTree.h"
#include <memory>
#include <string>

class ValenciaExc2p2hReweighter : public nusyst::IGENIESystProvider_tool {

public:
    explicit ValenciaExc2p2hReweighter(YAML::Node const &);

    bool SetupResponseCalculator(YAML::Node const &);

    // Extra options are saved in tool_options
    YAML::Node GetExtraToolOptions() { return tool_options; }
    systtools::SystMetaData BuildSystMetaData(YAML::Node const &, systtools::paramId_t);

    systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

    std::string AsString(){ return "ValenciaExc2p2hReweighter"; };

    ~ValenciaExc2p2hReweighter();

private:

    YAML::Node tool_options;
    size_t pidx_DialValencia;
    std::unique_ptr<BDTReweight::JSONReweighter> bdt_pp_reweighter;
    std::unique_ptr<BDTReweight::JSONReweighter> bdt_pn_reweighter;
    std::vector<double> BDTFeaturesWrapper(genie::EventRecord const &ev);
    double scale_factor_pp, scale_factor_pn;

};
