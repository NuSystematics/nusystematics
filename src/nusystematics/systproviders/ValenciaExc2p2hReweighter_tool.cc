// File: ValenciaExc2p2hReweighter_tool.cc
// Author: Zihao Lin zlin22@ur.rochester.edu

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/systproviders/ValenciaExc2p2hReweighter_tool.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "Framework/GHEP/GHepParticle.h"
#include "TLorentzVector.h"
#include <cmath>


using namespace systtools;
using namespace nusyst;
using namespace fhicl;

ValenciaExc2p2hReweighter::ValenciaExc2p2hReweighter(ParameterSet const &params)
    :IGENIESystProvider_tool(params),
    pidx_DialValencia(systtools::kParamUnhandled<size_t>),
    scale_factor_pp(1.0), scale_factor_pn(1.0)
    {}

SystMetaData ValenciaExc2p2hReweighter::BuildSystMetaData(ParameterSet const &cfg, paramId_t firstId) {

    std::cout << "[ValenciaExc2p2hReweighter::BuildSystMetaData] Called" << std::endl;
    SystMetaData smd;

    // Name of the dials that are supported by this module
    // std::vector<std::string> AvailPNames = {"DialA", "DialB"};
    std::vector<std::string> AvailPNames = {"DialValencia"};

    // Loop over available names and check if they are specified in ToolConfig
    for(std::string const &pname: AvailPNames){
        systtools::SystParamHeader phdr;
        if (ParseFhiclToolConfigurationParameter(cfg, pname, phdr, firstId)) {
            printf("[ValenciaExc2p2hReweighter::BuildSystMetaData] %s is found from ToolConfig\n", pname.c_str());
            phdr.systParamId = firstId++;
            smd.push_back(phdr);
        }
    }
    if(smd.size()==0){
        std::cout << "[ValenciaExc2p2hReweighter::BuildSystMetaData] No dial is set" << std::endl;
    }

    // Set path to BDT json reweighters for pp, pn out-going states.
    std::string JSON_PP_REWEIGHTER_PATH = cfg.get<std::string>("JSON_PP_REWEIGHTER_PATH", "");
    tool_options.put("JSON_PP_REWEIGHTER_PATH", JSON_PP_REWEIGHTER_PATH);
    std::string JSON_PN_REWEIGHTER_PATH = cfg.get<std::string>("JSON_PN_REWEIGHTER_PATH", "");
    tool_options.put("JSON_PN_REWEIGHTER_PATH", JSON_PN_REWEIGHTER_PATH);
    // Set normalization scale factors for pp, pn out-going states.
    // BDT reweight is a shape-only reweight; need to control total weights manually.
    // The default values are determined using 3M CCMEC events weights such that sum(weights) = number of pp, pn events.
    scale_factor_pp = cfg.get<double>("SCALE_FACTOR_PP", 2.4642232350453983);
    tool_options.put("SCALE_FACTOR_PP", scale_factor_pp);
    scale_factor_pn = cfg.get<double>("SCALE_FACTOR_PN", 1.6862982652296663);
    tool_options.put("SCALE_FACTOR_PN", scale_factor_pn);
    
    return smd;
}

bool ValenciaExc2p2hReweighter::SetupResponseCalculator(fhicl::ParameterSet const &tool_options) {

    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] Called" << std::endl;
    systtools::SystMetaData const &md = GetSystMetaData();

    if(HasParam(md, "DialValencia")){
        pidx_DialValencia = GetParamIndex(md, "DialValencia");
    }
    
    
    // Ziggy: read json bdt reweighter
    std::string JSON_PP_REWEIGHTER_PATH = tool_options.get<std::string>("JSON_PP_REWEIGHTER_PATH", "");
    std::string JSON_PN_REWEIGHTER_PATH = tool_options.get<std::string>("JSON_PN_REWEIGHTER_PATH", "");
    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] JSON_PP_REWEIGHTER_PATH:"<< JSON_PP_REWEIGHTER_PATH<< std::endl;
    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] JSON_PN_REWEIGHTER_PATH:"<< JSON_PN_REWEIGHTER_PATH<< std::endl;
    
    if (JSON_PP_REWEIGHTER_PATH.empty()) {
        throw std::runtime_error("ValenciaExc2p2hReweighter: JSON_PP_REWEIGHTER_PATH is empty");
    }
    bdt_pp_reweighter = std::make_unique<BDTReweight::JSONReweighter>(JSON_PP_REWEIGHTER_PATH);    
    if (JSON_PN_REWEIGHTER_PATH.empty()) {
        throw std::runtime_error("ValenciaExc2p2hReweighter: JSON_PN_REWEIGHTER_PATH is empty");
    }
    bdt_pn_reweighter = std::make_unique<BDTReweight::JSONReweighter>(JSON_PN_REWEIGHTER_PATH);    

    scale_factor_pp = tool_options.get<double>("SCALE_FACTOR_PP", 2.4642232350453983);
    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] use SCALE_FACTOR_PP = "<< scale_factor_pp << std::endl;
    scale_factor_pn = tool_options.get<double>("SCALE_FACTOR_PN", 1.6862982652296663);
    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] use SCALE_FACTOR_PN = "<< scale_factor_pn << std::endl;
    
    return true;
}

event_unit_response_t ValenciaExc2p2hReweighter::GetEventResponse(genie::EventRecord const &ev) {
    
    // Process info
    genie::ProcessInfo const& procinfo = ev.Summary()->ProcInfo();
    bool IsCC = procinfo.IsWeakCC();
    bool IsMEC = procinfo.IsMEC();
    if(!IsCC || !IsMEC){
        return this->GetDefaultEventResponse();
    }

    std::vector<double> features = BDTFeaturesWrapper(ev);
    
    systtools::event_unit_response_t resp;
    systtools::SystMetaData const &md = GetSystMetaData();
    
    double weight = 1.0;    
    if (ev.Particle(5)->Pdg()==2000000202){ //pp state
        weight = (bdt_pp_reweighter->PredictWeight(features)) * scale_factor_pp;
    } else if (ev.Particle(5)->Pdg()==2000000201){ //pn state
        weight = (bdt_pn_reweighter->PredictWeight(features)) * scale_factor_pn;
    }
    if (pidx_DialValencia != systtools::kParamUnhandled<size_t>) {
        resp.push_back( {md[pidx_DialValencia].systParamId, {}} );
        for (double var : md[pidx_DialValencia].paramVariations) {
            // Implement variation: assume 0 <= var <= 1; 0 <= weight < infinity
            // var = 0 is CV, the default SuSAv2 -> varied_weight = 1.0
            // var = 1 is fully reweighted to Valencia exclusive -> varied_weight = weight
            // varied_weight = 1 + var * (weight - 1)
            resp.back().responses.push_back(1 + var * (weight - 1));
        } 
    }
    
    return resp;
}


std::vector<double> ValenciaExc2p2hReweighter::BDTFeaturesWrapper(genie::EventRecord const &ev) {

    genie::GHepParticle *muon = ev.Particle(4);
    // CCMEC vertex lepton: index 4
    // CCMEC vertex intermediate particle: index 5
    // CCMEC vertex out-going Nucleons, pre-FSI: indices 6, 7
    genie::GHepParticle *mother = ev.Particle(5);
    genie::GHepParticle *N1, *N2;
    if (mother->Pdg()==2000000201){ //n+p state
        if (ev.Particle(6)->Pdg()==2212){
            N1 = ev.Particle(6);
            N2 = ev.Particle(7);
        }else{
            N1 = ev.Particle(7);
            N2 = ev.Particle(6);    
        }
    }else if (mother->Pdg()==2000000202){ //p+p state
        if (ev.Particle(6)->KinE() > ev.Particle(7)->KinE()){
            N1 = ev.Particle(6);
            N2 = ev.Particle(7);
        }else{
            N1 = ev.Particle(7);
            N2 = ev.Particle(6);
        }
    }else{
        std::cout << "[ValenciaExc2p2hReweighter::BDTFeaturesWrapper] unknown intermediate mother particle: "<< mother->Pdg() << std::endl;
    }
    
    const double nucleon1_px = N1->Px();
    const double nucleon1_py = N1->Py();
    const double nucleon1_pz = N1->Pz();
    const double nucleon2_px = N2->Px();
    const double nucleon2_py = N2->Py();
    const double nucleon2_pz = N2->Pz();
    const double muon_px = muon->Px();
    const double muon_py = muon->Py();
    const double muon_pz = muon->Pz();
    
    double muon_py_new = - std::sqrt(muon_px*muon_px + muon_py*muon_py); // new muon py is negative of muon transverse momentum magnitude
    std::vector<double> vector_y = {muon_px / muon_py_new, muon_py / muon_py_new}; // unit vector y of reaction plane frame coordinates
    std::vector<double> vector_x = {-vector_y[1], vector_y[0]}; //vector x is rotating y clockwise by 90 deg
    double nucleon1_px_new = nucleon1_px * vector_x[0] + nucleon1_py * vector_x[1];
    double nucleon1_py_new = nucleon1_px * vector_y[0] + nucleon1_py * vector_y[1];
    double nucleon2_px_new = nucleon2_px * vector_x[0] + nucleon2_py * vector_x[1];
    double nucleon2_py_new = nucleon2_px * vector_y[0] + nucleon2_py * vector_y[1];
    
    std::vector<double> features = {
        nucleon1_px_new, nucleon1_py_new, nucleon1_pz, 
        nucleon2_px_new, nucleon2_py_new, nucleon2_pz, 
        muon_py_new, muon_pz
    };
    
    return features;
};


ValenciaExc2p2hReweighter::~ValenciaExc2p2hReweighter() {
} 
