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
      pidx_surface_thickness(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

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

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  std::string OPT_STRING = cfg["OPT_STRING"] ? cfg["OPT_STRING"].as<std::string>() : ""; // second argument is the default when OPT_STRING does not exist
  tool_options["OPT_STRING"] = OPT_STRING;

  bool OPT_BOOL = cfg["OPT_BOOL"] ? cfg["OPT_BOOL"].as<bool>() : false;
  tool_options["OPT_BOOL"] = OPT_BOOL;

  fill_valid_tree = cfg["fill_valid_tree"] ? cfg["fill_valid_tree"].as<bool>() : false;
  tool_options["fill_valid_tree"] = fill_valid_tree;

  return smd;
}

bool WSReweight::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  std::cout << "[WSReweight::SetupResponseCalculator] OPT_STRING = " << (tool_options["OPT_STRING"] ? tool_options["OPT_STRING"].as<std::string>() : "") << std::endl;
  std::cout << "[WSReweight::SetupResponseCalculator] OPT_BOOL = " << (tool_options["OPT_BOOL"] ? tool_options["OPT_BOOL"].as<bool>() : false) << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();



  if (HasParam(md, "nucleus_radius")) {
    pidx_nucleus_radius = GetParamIndex(md, "nucleus_radius");
  }

  if (HasParam(md, "surface_thickness")) {
    pidx_surface_thickness = GetParamIndex(md, "surface_thickness");
  }

  fill_valid_tree = tool_options["fill_valid_tree"] ? tool_options["fill_valid_tree"].as<bool>() : false;
  if (fill_valid_tree) {
    InitValidTree();
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

  if (fill_valid_tree) {

    pdgfslep = ev.FinalStatePrimaryLepton()->Pdg();
    momfslep = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();

    Pdgnu = ISLep->Pdg();
    NEUTMode = 0;
    NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);

    QELikeTarget_t qel_targ = GetQELikeTarget(ev);
    QELTarget = e2i(qel_targ);

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    Emiss = Emiss_preFSI;
    if(nucleon){
      pmiss = nucleon->P4()->Vect();
    }
    else{
      pmiss = TLorentzVector(0,0,0,0).Vect();
    }
    KF_tree = KF;
    radius = GetRadiusFromKF(KF, isProton);
    radius = std::max(radius, 0.0); // No negative radii!
    ref_prob_density = genie::utils::nuclear::Density(radius, 40);
    new_prob_density = genie::utils::nuclear::DensityWoodsSaxon(radius, 4, kAr40SkinDepth);

    valid_tree->Fill();
  }

  return resp;
}

std::string WSReweight::AsString() { return ""; }

void WSReweight::InitValidTree() {

  valid_file = new TFile("WSReweightWeights_validTree.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("QELTarget", &QELTarget);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("pmiss", &pmiss);
  valid_tree->Branch("Emiss", &Emiss);
  valid_tree->Branch("KF", &KF_tree);
  valid_tree->Branch("radius", &radius);
  valid_tree->Branch("ref_prob_density", &ref_prob_density);
  valid_tree->Branch("new_prob_density", &new_prob_density);
}

WSReweight::~WSReweight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
} 
