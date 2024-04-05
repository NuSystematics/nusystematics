#include "nusystematics/systproviders/SPPTpiReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

SPPTpiReweight::SPPTpiReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_SPPTpiCVCorrection(systtools::kParamUnhandled<size_t>),
      pidx_SPPTpiCorrectionRW(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData SPPTpiReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[SPPTpiReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;
  for (std::string const &pname :
       {"SPPTpiCVCorrection", "SPPTpiCorrectionRW"}){
    systtools::SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(cfg, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool SPPTpiReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[SPPTpiReweight::SetupResponseCalculator] called" << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "SPPTpiCVCorrection")) {
    pidx_SPPTpiCVCorrection = 
        GetParamIndex(md, "SPPTpiCVCorrection");

    // this must be "isCorrection"
    if( ! GetParam(md, pidx_SPPTpiCVCorrection).isCorrection ){

      throw invalid_engine_state()
          << "SPPTpiCVCorrection from SPPTpiReweight module must be a correction dial, but it is not. Check your config/paramheader file";
    }

  }
  if (HasParam(md, "SPPTpiCorrectionRW")) {
    pidx_SPPTpiCorrectionRW =
        GetParamIndex(md, "SPPTpiCorrectionRW");
  }

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
SPPTpiReweight::GetEventResponse(genie::EventRecord const &ev) {

  // when the event is not applicable for this type of reweighting,
  // use GetDefaultEventResponse() to return an auto-1.-filled vector

  bool IsCOH = ev.Summary()->ProcInfo().IsCoherentProduction();
  if(IsCOH){
    return this->GetDefaultEventResponse();
  }

  // loop over particles
  int ip=-1;
  genie::GHepParticle * p = 0;
  TIter event_iter(&ev);

  // find highest momentum final-state particle
  double max_mom_fs_pip = -999.;
  genie::GHepParticle *ptl_hm_fs_pip = 0;
  int genie_n_photons = 0;
  int genie_n_mesons = 0;
  int nPip=0;

  while ( (p = dynamic_cast<genie::GHepParticle *>(event_iter.Next())) ) {

    ip++;

    // Skip particles not rescattered by the actual hadron transport code
    int  pdgc       = p->Pdg();
    bool is_pion    = genie::pdg::IsPion   (pdgc);
    genie::GHepStatus_t ist  = p->Status();
    if( ist!=genie::kIStStableFinalState ) continue;

    // photon
    if( pdgc==22 ){
      TLorentzVector* ptl_P4 = p->P4();
      if(ptl_P4->E() > 0.010){
        genie_n_photons++;
      }
    }

    // count all mesons
    if (abs(pdgc) == 211 || //pi+-
             pdgc == 111 ||  // pi0
             abs(pdgc) == 321 || // K-
             abs(pdgc) == 323 || // K*+-
             pdgc == 130 || // KL0
             pdgc == 310 || // KS0
             pdgc == 311 || // K0
             pdgc == 313 || // K*0
             abs(pdgc) == 221 || // eta
             abs(pdgc) == 331 // eta' (958)
             ) {
      genie_n_mesons++;
    }

    // count pip and save thet momentum
    if( pdgc==genie::kPdgPiP ){

      nPip++;

      TLorentzVector* ptl_P4 = p->P4();
      double mom_mag = ptl_P4->Vect().Mag();

      if( mom_mag > max_mom_fs_pip ){
        ptl_hm_fs_pip = p;
      }

    }
    else{
      continue;
    }

  }//p

  // SPP
  if( nPip != 1 || genie_n_mesons!= 1 ){
    return this->GetDefaultEventResponse();
  }

  // has photon with E>10 MeV; following MINERvA CC1pip signal definition
  if(genie_n_photons!=0){
    return this->GetDefaultEventResponse();
  }


  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double this_Tpi_GeV = ptl_hm_fs_pip->KinE();
  double this_Q2_GeV2 = -emTransfer.Mag2();
  //std::cout << "Leading pi+ T = " << this_Tpi_GeV << " GeV" << std::endl;

  // now make the output
  systtools::event_unit_response_t resp;

  systtools::SystMetaData const &md = GetSystMetaData();

  int TargetA = ev.Summary()->InitState().Tgt().A();
  bool IsH = TargetA==1;

  if (pidx_SPPTpiCVCorrection != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_SPPTpiCVCorrection].systParamId, {}} );
    double this_reweight = GetSPPTpiCVCorrection(this_Q2_GeV2, this_Tpi_GeV);
    if(IsH) resp.back().responses.push_back( 1. );
    else resp.back().responses.push_back( this_reweight );
  }

  if (pidx_SPPTpiCorrectionRW != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_SPPTpiCorrectionRW].systParamId, {}} );
    for (double var : md[pidx_SPPTpiCorrectionRW].paramVariations) {
      double this_reweight = GetSPPTpiCorrectionRW(this_Q2_GeV2, this_Tpi_GeV, var);
      if(IsH) resp.back().responses.push_back( 1. );
      else resp.back().responses.push_back( this_reweight );
    }
  }


  if (fill_valid_tree) {

    pdgfslep = ev.FinalStatePrimaryLepton()->Pdg();
    momfslep = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();

    Pdgnu = ISLep->Pdg();
    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

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

std::string SPPTpiReweight::AsString() { return ""; }

void SPPTpiReweight::InitValidTree() {

  valid_file = new TFile("MINERvAq3q0WeightValid.root", "RECREATE");
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
}

SPPTpiReweight::~SPPTpiReweight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
