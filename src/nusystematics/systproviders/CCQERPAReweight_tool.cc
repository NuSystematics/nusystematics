#include "nusystematics/systproviders/CCQERPAReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

CCQERPAReweight::CCQERPAReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      ccqeRPAReweightCalculator(nullptr),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData CCQERPAReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[CCQERPAReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "CCQERPAReweight",
                                                 phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet templateManifest =
      cfg.get<fhicl::ParameterSet>("CCQERPAReweight_input_manifest");
  tool_options.put("CCQERPAReweight_input_manifest", templateManifest);

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool CCQERPAReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[CCQERPAReweight::SetupResponseCalculator] called" << std::endl;

  fhicl::ParameterSet templateManifest =
      tool_options.get<fhicl::ParameterSet>("CCQERPAReweight_input_manifest");

  std::string rwmode_str = templateManifest.get<std::string>("RWMode");
  std::string kin_Y_str(""), kin_Z_str("");
  if(rwmode_str=="q3q0"){
    rwMode = q3q0;
    kin_Y_str = "q3";
    kin_Z_str = "q0";
  }
  else if(rwmode_str=="PCTheta"){
    rwMode = PCTheta;
    kin_Y_str = "P";
    kin_Z_str = "CTheta";
  }
  else if(rwmode_str=="PSTheta"){
    rwMode = PSTheta;
    kin_Y_str = "P";
    kin_Z_str = "STheta";
  }
  else if(rwmode_str=="PTheta"){
    rwMode = PTheta;
    kin_Y_str = "P";
    kin_Z_str = "Theta";
  }
  else{
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: RWMode is wrong: " << rwmode_str;
  }

  std::cout << "[CCQERPAReweight::SetupResponseCalculator] RWMode: " << rwmode_str << std::endl;
  std::cout << "[CCQERPAReweight::SetupResponseCalculator] Template binnings are" << std::endl;
  std::cout << "[CCQERPAReweight::SetupResponseCalculator] x: E_nu" << std::endl;
  std::cout << "[CCQERPAReweight::SetupResponseCalculator] y: " << kin_Y_str << std::endl;
  std::cout << "[CCQERPAReweight::SetupResponseCalculator] z: " << kin_Z_str << std::endl;

  ccqeRPAReweightCalculator = std::make_unique<CCQERPAReweightCalculator>( templateManifest );

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "CCQERPAReweight");

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
CCQERPAReweight::GetEventResponse(genie::EventRecord const &ev) {

  // when the event is not applicable for this type of reweighting,
  // use GetDefaultEventResponse() to return an auto-1.-filled vector

  if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return this->GetDefaultEventResponse();
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  TLorentzVector FSLepP4 = *FSLep->P4(); // l
  TLorentzVector ISLepP4 = *ISLep->P4(); // nu
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double AngleLeps = FSLepP4.Vect().Angle( ISLepP4.Vect() );
  double CAngleLeps = TMath::Cos(AngleLeps);
  double SAngleLeps = TMath::Sin(AngleLeps);
  if(AngleLeps>=M_PI/2.) SAngleLeps *= -1.;

  std::array<double, 2> bin_kin;
  if(rwMode==q3q0){
    bin_kin = {emTransfer.Vect().Mag(), emTransfer.E()};
  }
  else if(rwMode==PCTheta){
    bin_kin = {FSLepP4.Vect().Mag(), CAngleLeps};
  }
  else if(rwMode==PSTheta){
    bin_kin = {FSLepP4.Vect().Mag(), SAngleLeps};
  }
  else if(rwMode==PTheta){
    bin_kin = {FSLepP4.Vect().Mag(), AngleLeps};
  }
  else{
    throw invalid_ToolConfigurationFHiCL() 
        << "[ERROR]: RWMode is wrong: " << rwMode;
  }

  // now make the output
  systtools::event_unit_response_t resp;

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  resp.push_back( {hdr.systParamId, {}} );
  for (double var : hdr.paramVariations) {
    double this_reweight = ccqeRPAReweightCalculator->GetRPAReweight( 
      ISLepP4.E(),
      bin_kin,
      var
    );
    resp.back().responses.push_back( this_reweight );
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

std::string CCQERPAReweight::AsString() { return ""; }

void CCQERPAReweight::InitValidTree() {

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

CCQERPAReweight::~CCQERPAReweight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
