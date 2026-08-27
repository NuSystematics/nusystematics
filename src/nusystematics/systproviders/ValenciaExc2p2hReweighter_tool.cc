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
    : IGENIESystProvider_tool(params),
      pidx_DialA(systtools::kParamUnhandled<size_t>),
      pidx_DialB(systtools::kParamUnhandled<size_t>){

}

SystMetaData ValenciaExc2p2hReweighter::BuildSystMetaData(ParameterSet const &cfg, paramId_t firstId) {

    std::cout << "[ValenciaExc2p2hReweighter::BuildSystMetaData] Called" << std::endl;
    SystMetaData smd;

    // Name of the dials that are supported by this module
    std::vector<std::string> AvailPNames = {"DialA", "DialB"};

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

    // You can extra parameters for the module;
    // Use get<T> function to retrieve the value from ToolConfig,
    // and then run "put" on "tool_options" (defined in the header of this class as a  member variable)
    // T can be string, bool, int, unsigned, float, double, std::string or even a new fhicl::ParameterSet
    std::string OPT_STRING = cfg.get<std::string>("OPT_STRING", ""); // second argument is the default when OPT_STRING does not exist
    tool_options.put("OPT_STRING", OPT_STRING);
    bool OPT_BOOL = cfg.get<bool>("OPT_BOOL", false);
    tool_options.put("OPT_BOOL", OPT_BOOL);
    fhicl::ParameterSet OPT_PSET = cfg.get<fhicl::ParameterSet>("OPT_PSET");
    tool_options.put("OPT_PSET", OPT_PSET);
    
    return smd;
}

bool ValenciaExc2p2hReweighter::SetupResponseCalculator(fhicl::ParameterSet const &tool_options) {

    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] Called" << std::endl;
    systtools::SystMetaData const &md = GetSystMetaData();
    
    if(HasParam(md, "DialA")){
        pidx_DialA = GetParamIndex(md, "DialA");
    }
    
    if(HasParam(md, "DialB")){
        pidx_DialB = GetParamIndex(md, "DialB");
    }
    
    // Ziggy: read json bdt reweighter
    std::string JSON_REWEIGHTER_PATH = tool_options.get<std::string>("JSON_REWEIGHTER_PATH", "");
    std::cout << "[ValenciaExc2p2hReweighter::SetupResponseCalculator] JSON_REWEIGHTER_PATH:"<< JSON_REWEIGHTER_PATH<< std::endl;
    if (JSON_REWEIGHTER_PATH.empty()) {
        throw std::runtime_error("ValenciaExc2p2hReweighter: JSON_REWEIGHTER_PATH is empty");
    }
    bdt_reweighter = std::make_unique<BDTReweight::JSONReweighter>(JSON_REWEIGHTER_PATH);    

    // Parameters in tool_options
    std::string OPT_STRING = tool_options.get<std::string>("OPT_STRING", "");
    bool OPT_BOOL = tool_options.get<bool>("OPT_BOOL", false);
    fhicl::ParameterSet OPT_PSET = tool_options.get<fhicl::ParameterSet>("OPT_PSET");
    std::string OPT_ROOTFileName = OPT_PSET.get<std::string>("ROOTFileName", "");
    std::string OPT_HistName = OPT_PSET.get<std::string>("HistName", "");
    if(OPT_ROOTFileName!="" && OPT_HistName!=""){
        printf("[ValenciaExc2p2hReweighter::SetupResponseCalculator] ROOTFileName: %s\n", OPT_ROOTFileName.c_str());
        printf("[ValenciaExc2p2hReweighter::SetupResponseCalculator] HistName: %s\n", OPT_HistName.c_str());
    }
    
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
    
    if (pidx_DialA != systtools::kParamUnhandled<size_t>) {
        resp.push_back( {md[pidx_DialA].systParamId, {}} );
        for (double var : md[pidx_DialA].paramVariations) {
            // resp.back().responses.push_back( GetReweight_DialA(Q2, var) );
            resp.back().responses.push_back( bdt_reweighter->PredictWeight(features) );
        } 
    }
    
    if (pidx_DialB != systtools::kParamUnhandled<size_t>) {
        resp.push_back( {md[pidx_DialB].systParamId, {}} );
        for (double var : md[pidx_DialB].paramVariations) {
            // resp.back().responses.push_back( GetReweight_DialB(Q2, var) );
            resp.back().responses.push_back( bdt_reweighter->PredictWeight(features) );
        }
    }
    
    return resp;
}


std::vector<double> ValenciaExc2p2hReweighter::BDTFeaturesWrapper(genie::EventRecord const &ev) {
    // const double nucleon1_px, const double nucleon1_py, const double nucleon1_pz,
    // const double nucleon2_px, const double nucleon2_py, const double nucleon2_pz,
    // const double muon_px, const double muon_py, const double muon_pz

    // CCMEC vertex muon: indices 4
    // CCMEC vertex out-going Nucleons: indices 6, 7
    genie::GHepParticle *muon = ev.Particle(4);
    genie::GHepParticle *N1 = ev.Particle(6);
    genie::GHepParticle *N2 = ev.Particle(7);
    
    double nucleon1_px = N1->Px();
    double nucleon1_py = N1->Py();
    double nucleon1_pz = N1->Pz();
    double nucleon2_px = N2->Px();
    double nucleon2_py = N2->Py();
    double nucleon2_pz = N2->Pz();
    double muon_px = muon->Px();
    double muon_py = muon->Py();
    double muon_pz = muon->Pz();
    std::cout << "[ValenciaExc2p2hReweighter::BDTFeaturesWrapper] muon, N1, N2 pdg:"<< muon->Pdg() << N1->Pdg() << N2->Pdg() << std::endl;

    double muon_py_new = - std::sqrt(muon_px*muon_px + muon_py*muon_py);
    std::vector<double> vector_y = {muon_px / muon_py_new, muon_py / muon_py_new};
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
