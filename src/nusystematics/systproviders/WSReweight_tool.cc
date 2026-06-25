#include "nusystematics/systproviders/WSReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/responsecalculators/WSReweight_calculator.hh"

#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace systtools;

WSReweight::WSReweight(YAML::Node const &params)
    : IGENIESystProvider_tool(params),
      pidx_nucleus_radius(systtools::kParamUnhandled<size_t>),
      pidx_surface_thickness(systtools::kParamUnhandled<size_t>) {}

SystMetaData WSReweight::BuildSystMetaData(YAML::Node const &cfg,
                                                     paramId_t firstId) {

  SystMetaData smd;

  for (std::string const &pname :
       {"nucleus_radius", "surface_thickness"}) {
    systtools::SystParamHeader phdr;
    if (ParseYAMLToolConfigurationParameter(cfg, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  if(cfg["estimate_emiss"]){
    tool_options["estimate_emiss"] = cfg["estimate_emiss"].as<bool>();
  }

  return smd;
}

bool WSReweight::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  systtools::SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "nucleus_radius")) {
    pidx_nucleus_radius = GetParamIndex(md, "nucleus_radius");
  }

  if (HasParam(md, "surface_thickness")) {
    pidx_surface_thickness = GetParamIndex(md, "surface_thickness");
  }

  estimate_emiss = tool_options["estimate_emiss"] ? tool_options["estimate_emiss"].as<bool>() : false;

  return true;
}

event_unit_response_t
WSReweight::GetEventResponse(genie::EventRecord const &ev) {

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  // TH: get Ermv directly from GENIE event record
  double Emiss_preFSI, Pmiss_preFSI;
  double Q = 0.01;
  int nucleon_PDG, target_PDG;
  double KF = 0;

    // GHepRecord::TargetNucleus() is designed to return nullptr for a free nucleon target (e.g., hydrogen)
  // If TargetNucleus() is available, use ev.TargetNucleus()->Pdg()
  // if not, force it to hydrogen
  target_PDG = ev.TargetNucleus() ? ev.TargetNucleus()->Pdg() : 1000010010;
  
  genie::GHepParticle *nucleon = ev.HitNucleon();
  if (nucleon == NULL){
    // TH: some events don't have an initial nucleon (e.g. coherent scattering)
    //     want to skip these events and not re-weight
    Emiss_preFSI = -999;
    nucleon_PDG = -999;
    target_PDG = -999;
  } else if (nucleon->Mass() > 1) //Getting rid of 2p2h for now
  {
    Emiss_preFSI = -999;
    nucleon_PDG = -999;
    target_PDG = -999;
  }
  else {
    if(estimate_emiss){
      Emiss_preFSI = nucleon->Mass() - nucleon->Energy();
    }
    else{
      Emiss_preFSI = nucleon->RemovalEnergy();
    }
    
    // Pmiss_preFSI = nucleon->P4()->Vect().Mag();
    Pmiss_preFSI = sqrt(nucleon->Px() * nucleon->Px() + nucleon->Py() * nucleon->Py() + nucleon->Pz() * nucleon->Pz());
    nucleon_PDG = nucleon->Pdg();
    KF = sqrt(pow(Emiss_preFSI + sqrt(pow(Pmiss_preFSI, 2) + pow(nucleon->Mass(), 2)) - Q, 2) - pow(nucleon->Mass(), 2));
  }

  bool isProton = (nucleon_PDG == 2212);

  // std::cout << *nucleon << std::endl;
  // std::cout << "Mass: " << ev.Summary()->InitState().Tgt().HitNucMass() << std::endl;

  // now make the output
  systtools::event_unit_response_t resp;
  systtools::SystMetaData const &md = GetSystMetaData();

  if (pidx_nucleus_radius != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_nucleus_radius].systParamId, {}} );
    if (target_PDG == 1000180400){
      for (double var : md[pidx_nucleus_radius].paramVariations) {
        resp.back().responses.push_back( GetWeightFomKF(KF, var, kAr40SkinDepth, isProton) );
      } 
    }
    else{
      for (unsigned int i = 0; i < md[pidx_nucleus_radius].paramVariations.size(); i++) {
        resp.back().responses.push_back(1);
      }
    }
  }

  if (pidx_surface_thickness != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_surface_thickness].systParamId, {}} );
    if (target_PDG == 1000180400){
      for (double var : md[pidx_surface_thickness].paramVariations) {
        resp.back().responses.push_back( GetWeightFomKF(KF, kAr40Radius, var, isProton) );
      } 
    }
    else{
      for (unsigned int i = 0; i < md[pidx_surface_thickness].paramVariations.size(); i++) {
        resp.back().responses.push_back(1);
      }
    }
  }

  return resp;
}

std::string WSReweight::AsString() { return ""; }

WSReweight::~WSReweight() {
} 
