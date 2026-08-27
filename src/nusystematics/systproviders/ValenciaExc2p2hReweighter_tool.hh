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
    explicit ValenciaExc2p2hReweighter(fhicl::ParameterSet const &);
    
    bool SetupResponseCalculator(fhicl::ParameterSet const &);

    // Extra options are saved in tool_options
    fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }
    systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &, systtools::paramId_t);

    systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

    std::string AsString(){ return "ValenciaExc2p2hReweighter"; };

    ~ValenciaExc2p2hReweighter();

private:

    fhicl::ParameterSet tool_options;
    size_t pidx_DialA, pidx_DialB;
    std::unique_ptr<BDTReweight::JSONReweighter> bdt_reweighter;
    std::vector<double> BDTFeaturesWrapper(genie::EventRecord const &ev);

};
