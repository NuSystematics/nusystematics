#include "nusystematics/systproviders/CCQETemplateReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;


CCQETemplateReweight::CCQETemplateReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      ccqeTemplateReweightCalculator(nullptr),
      valid_file(nullptr), valid_tree(nullptr) {

  // knob names
  for (size_t i = 0; i < kNumQ0Bins; ++i) {
      q0bin_params_for_ccqe_sf_names[i] = "q0bin" + std::to_string(i + 1);
  }

  for (size_t i = 0; i < kNumQ0Bins; ++i) {
    ResponseParameterIdx_q0bin[i] = systtools::kParamUnhandled<size_t>;
  }
}

SystMetaData CCQETemplateReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[CCQETemplateReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  for (size_t i = 0; i < kNumQ0Bins; ++i) {
    systtools::SystParamHeader phdr;
    std::cout << "[CCQETemplateReweight::BuildSystMetaData] Attempting to parse parameter: " 
              << q0bin_params_for_ccqe_sf_names[i] << std::endl;
    if (ParseFhiclToolConfigurationParameter(cfg, q0bin_params_for_ccqe_sf_names[i],
                                                   phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
      std::cout << "[CCQETemplateReweight::BuildSystMetaData] Successfully parsed parameter: " 
                << q0bin_params_for_ccqe_sf_names[i] << " with ID: " << phdr.systParamId << std::endl;
    } else {
      std::cout << "[CCQETemplateReweight::BuildSystMetaData] Failed to parse parameter: " 
                << q0bin_params_for_ccqe_sf_names[i] << std::endl;
    }
  }

  fhicl::ParameterSet templateManifest =
      cfg.get<fhicl::ParameterSet>("CCQETemplateReweight_input_manifest");
  tool_options.put("CCQETemplateReweight_input_manifest", templateManifest);

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  std::cout << "[CCQETemplateReweight::BuildSystMetaData] Total parameters parsed: " << smd.size() << std::endl;

  return smd;
}

bool CCQETemplateReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] called" << std::endl;

  fhicl::ParameterSet templateManifest =
      tool_options.get<fhicl::ParameterSet>("CCQETemplateReweight_input_manifest");

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

  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] RWMode: " << rwmode_str << std::endl;
  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] Template binnings are" << std::endl;
  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] x: E_nu" << std::endl;
  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] y: " << kin_Y_str << std::endl;
  std::cout << "[CCQETemplateReweight::SetupResponseCalculator] z: " << kin_Z_str << std::endl;

  ccqeTemplateReweightCalculator = std::make_unique<CCQETemplateReweightCalculator>( templateManifest );

  for (size_t i = 0; i < kNumQ0Bins; ++i) {
    ResponseParameterIdx_q0bin[i] = GetParamIndex(GetSystMetaData(), q0bin_params_for_ccqe_sf_names[i]);
    std::cout << "[CCQETemplateReweight::SetupResponseCalculator] " << q0bin_params_for_ccqe_sf_names[i] 
              << " parameter index: " << ResponseParameterIdx_q0bin[i] << std::endl;
  }

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
CCQETemplateReweight::GetEventResponse(genie::EventRecord const &ev) {

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
    //bin_kin = {emTransfer.Vect().Mag(), emTransfer.E()};
    // TEST; cutoff for non-physical reweights
    double this_q3 = emTransfer.Vect().Mag();
    double this_q0 = emTransfer.E()>q0_bin_boundaries[kNumQ0Bins] ? q0_bin_boundaries[kNumQ0Bins] : emTransfer.E();
    bin_kin = {this_q3, this_q0};
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

  // Determine which q0 bin this event belongs to
  double q0_value = emTransfer.E();
  int q0_bin_index = GetQ0BinIndex(q0_value);
  // std::cout << "[CCQETemplateReweight::GetEventResponse] q0=" << q0_value 
  //          << " -> bin " << q0_bin_index << " -> resp index " << ResponseParameterIdx_q0bin[q0_bin_index] << std::endl;

  systtools::event_unit_response_t resp;

  // put in default (1) for all of the weights
  for (size_t i = 0; i < kNumQ0Bins; ++i) {
    SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx_q0bin[i]];
    resp.push_back( {hdr.systParamId, {}} );
    for (double var : hdr.paramVariations) resp.back().responses.push_back(1.);
  }

  // reweighting for the corresponding q0 bin
  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx_q0bin[q0_bin_index]];
  unsigned i_var = 0;
  for (double var : hdr.paramVariations) {
    double this_reweight = ccqeTemplateReweightCalculator->GetTemplateReweight( 
      ISLepP4.E(),
      bin_kin,
      var
    );
    resp[q0_bin_index].responses[i_var] = this_reweight;
    i_var ++;
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

std::string CCQETemplateReweight::AsString() { return ""; }

int CCQETemplateReweight::GetQ0BinIndex(double q0_value) const {

  double q0_cutoff = (q0_value > q0_bin_boundaries[kNumQ0Bins]) ? q0_bin_boundaries[kNumQ0Bins] : q0_value;
  
  for (size_t i = 0; i < kNumQ0Bins; ++i) {
    if (q0_cutoff >= q0_bin_boundaries[i] && q0_cutoff < q0_bin_boundaries[i + 1]) {
      return static_cast<int>(i);
    }
  }
  
  if (q0_cutoff >= q0_bin_boundaries[kNumQ0Bins]) {
    return kNumQ0Bins;
  }
  
  // Fallback 
  // TODO: shouldn't happen -- but should we add a handler?
  std::cout << "[CCQETemplateReweight::GetQ0BinIndex] WARNING: q0=" << q0_value 
            << " (cutoff=" << q0_cutoff << ") not in any bin, using bin 0" << std::endl;
  return 0;
}

void CCQETemplateReweight::InitValidTree() {

  valid_file = new TFile("CCQETemplateReweight.root", "RECREATE");
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

CCQETemplateReweight::~CCQETemplateReweight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
