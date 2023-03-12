#include "nusystematics/systproviders/DIRT2_Emiss_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#ifndef NO_ART
#include "art/Utilities/ToolMacros.h"
#endif

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

DIRT2_Emiss::DIRT2_Emiss(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_Emiss_CorrTail(systtools::kParamUnhandled<size_t>),
      pidx_Emiss_Linear(systtools::kParamUnhandled<size_t>),
      pidx_Emiss_ShiftPeak(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

#ifndef NO_ART
DEFINE_ART_CLASS_TOOL(DIRT2_Emiss)
#endif

SystMetaData DIRT2_Emiss::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  SystMetaData smd;

  for (std::string const &pname :
       {"Emiss_CorrTail", "Emiss_Linear", "Emiss_ShiftPeak"}) {
    systtools::SystParamHeader phdr;
    if (ParseFHiCLSimpleToolConfigurationParameter(cfg, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  std::string OPT_STRING = cfg.get<std::string>("OPT_STRING", ""); // second argument is the default when OPT_STRING does not exist
  tool_options.put("OPT_STRING", OPT_STRING);

  bool OPT_BOOL = cfg.get<bool>("OPT_BOOL", false);
  tool_options.put("OPT_BOOL", OPT_BOOL);

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool DIRT2_Emiss::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[DIRT2_Emiss::SetupResponseCalculator] OPT_STRING = " << tool_options.get<std::string>("OPT_STRING") << std::endl;
  std::cout << "[DIRT2_Emiss::SetupResponseCalculator] OPT_BOOL = " << tool_options.get<bool>("OPT_BOOL") << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "Emiss_CorrTail")) {
    pidx_Emiss_CorrTail =
        GetParamIndex(md, "Emiss_CorrTail");
  }

  if (HasParam(md, "Emiss_Linear")) {
    pidx_Emiss_Linear =
        GetParamIndex(md, "Emiss_Linear");
  }

  if (HasParam(md, "Emiss_ShiftPeak")) {
    pidx_Emiss_ShiftPeak =
        GetParamIndex(md, "Emiss_ShiftPeak");
  }

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
DIRT2_Emiss::GetEventResponse(genie::EventRecord const &ev) {

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  // TH: get Ermv directly from GENIE event record
  double Emiss_preFSI;
  genie::GHepParticle *nucleon = ev.HitNucleon();
  if (nucleon == NULL){
    Emiss_preFSI = -999;
  }
  else {
    Emiss_preFSI = nucleon->RemovalEnergy();
  }  
  
  // now make the output
  systtools::event_unit_response_t resp;
  systtools::SystMetaData const &md = GetSystMetaData();

  if (pidx_Emiss_CorrTail != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_Emiss_CorrTail].systParamId, {}} );
    for (double var : md[pidx_Emiss_CorrTail].paramVariations) {
      resp.back().responses.push_back( GetEmissCorrTailRW( Emiss_preFSI, var) );
    }
  }

  if (pidx_Emiss_Linear != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_Emiss_Linear].systParamId, {}} );
    for (double var : md[pidx_Emiss_Linear].paramVariations) {
      resp.back().responses.push_back( GetEmissLinearRW( Emiss_preFSI, var) );
    }
  }

  if (pidx_Emiss_ShiftPeak != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_Emiss_ShiftPeak].systParamId, {}} );
    for (double var : md[pidx_Emiss_ShiftPeak].paramVariations) {
      resp.back().responses.push_back( GetEmissShiftPeakRW( Emiss_preFSI, var) );
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

    valid_tree->Fill();
  }

  return resp;
}

std::string DIRT2_Emiss::AsString() { return ""; }

void DIRT2_Emiss::InitValidTree() {

  valid_file = new TFile("DIRT2_EmissWeights_validTree.root", "RECREATE");
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
}

DIRT2_Emiss::~DIRT2_Emiss() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
} 
