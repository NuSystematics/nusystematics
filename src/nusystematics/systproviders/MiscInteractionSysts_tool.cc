#include "nusystematics/systproviders/MiscInteractionSysts_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "nusystematics/responsecalculators/C12ToAr40_2p2hScaling.hh"
#include "nusystematics/responsecalculators/SPPLowQ2Suppression.hh"
#include "nusystematics/responsecalculators/nuenuebar_xsec_ratio.hh"
#include "nusystematics/responsecalculators/nuenumu_xsec_ratio.hh"

using namespace nusyst;

MiscInteractionSysts::MiscInteractionSysts(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_C12ToAr40_2p2hScaling_nu(systtools::kParamUnhandled<size_t>),
      pidx_C12ToAr40_2p2hScaling_nubar(systtools::kParamUnhandled<size_t>),
      pidx_nuenuebar_xsec_ratio(systtools::kParamUnhandled<size_t>),
      pidx_nuenumu_xsec_ratio(systtools::kParamUnhandled<size_t>),
      pidx_SPPLowQ2Suppression(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

systtools::SystMetaData
MiscInteractionSysts::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                        systtools::paramId_t firstId) {

  systtools::SystMetaData smd;

  for (std::string const &pname :
       {"C12ToAr40_2p2hScaling_nu", "C12ToAr40_2p2hScaling_nubar",
        "nuenuebar_xsec_ratio", "nuenumu_xsec_ratio", "SPPLowQ2Suppression"}) {
    systtools::SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(ps, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  fill_valid_tree = ps.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool MiscInteractionSysts::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  systtools::SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "C12ToAr40_2p2hScaling_nu")) {
    pidx_C12ToAr40_2p2hScaling_nu =
        GetParamIndex(md, "C12ToAr40_2p2hScaling_nu");
  }
  if (HasParam(md, "C12ToAr40_2p2hScaling_nubar")) {
    pidx_C12ToAr40_2p2hScaling_nubar =
        GetParamIndex(md, "C12ToAr40_2p2hScaling_nubar");
  }
  if (HasParam(md, "nuenuebar_xsec_ratio")) {
    pidx_nuenuebar_xsec_ratio = GetParamIndex(md, "nuenuebar_xsec_ratio");
  }
  if (HasParam(md, "nuenumu_xsec_ratio")) {
    pidx_nuenumu_xsec_ratio = GetParamIndex(md, "nuenumu_xsec_ratio");
  }
  if (HasParam(md, "SPPLowQ2Suppression")) {
    pidx_SPPLowQ2Suppression = GetParamIndex(md, "SPPLowQ2Suppression");
  }
  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);

  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

std::vector<double> MiscInteractionSysts::GetWeights_C12ToAr40_2p2hScaling(
    genie::EventRecord const &ev, std::vector<double> const &vals) {

  std::vector<double> resp;

  QELikeTarget_t mec_topology = GetQELikeTarget(ev);

  if ((mec_topology == nusyst::QELikeTarget_t::kQE) ||
      (mec_topology == nusyst::QELikeTarget_t::kInvalidTopology)) {
    return std::vector<double>(vals.size(), 1.);
  }

  for (double v : vals) {
    resp.push_back(GetC_Ar2p2hScalingWeight(v));
  }

  return resp;
}

std::vector<double> MiscInteractionSysts::GetWeights_nuenuebar_xsec_ratio(
    genie::EventRecord const &ev, std::vector<double> const &vals) {

  std::vector<double> resp;

  genie::GHepParticle *ISLep = ev.Probe();
  if (abs(ISLep->Pdg()) != 12) {
    return std::vector<double>(vals.size(), 1.);
  }

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return std::vector<double>(vals.size(), 1.);
  }

  int pdgnu = ISLep->Pdg();
  double enu = ISLep->P4()->E();

  for (double v : vals) {
    resp.push_back(GetNueNueBarXSecRatioWeight(pdgnu, true, enu, v));
  }

  return resp;
}
std::vector<double> MiscInteractionSysts::GetWeights_nuenumu_xsec_ratio(
    genie::EventRecord const &ev, std::vector<double> const &vals) {

  std::vector<double> resp;

  genie::GHepParticle *ISLep = ev.Probe();
  if (abs(ISLep->Pdg()) != 12) {
    return std::vector<double>(vals.size(), 1.);
  }

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return std::vector<double>(vals.size(), 1.);
  }

  int pdgnu = ISLep->Pdg();

  TLorentzVector ISLepP4 = *ISLep->P4();

  double enu = ISLepP4.E();

  TLorentzVector FSLepP4 = *ev.FinalStatePrimaryLepton()->P4();

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double q0_GeV = emTransfer.E();
  double q3_GeV = emTransfer.Vect().Mag();

  for (double v : vals) {
    resp.push_back(GetNueNumuRatioWeight(pdgnu, true, enu, q0_GeV, q3_GeV, v));
  }

  return resp;
}
std::vector<double> MiscInteractionSysts::GetWeights_SPPLowQ2Suppression(
    genie::EventRecord const &ev, std::vector<double> const &vals) {

  std::vector<double> resp;

  if (SPPChannelFromGHep(ev) == genie::kSppNull) {
    return std::vector<double>(vals.size(), 1.);
  }

  TLorentzVector ISLepP4 = *ev.Probe()->P4();

  TLorentzVector FSLepP4 = *ev.FinalStatePrimaryLepton()->P4();

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
  double Q2_GeV = -emTransfer.Mag2();

  for (double v : vals) {
    resp.push_back(GetMINERvASPPLowQ2SuppressionWeight(e2i(GetSimbMode(ev)),
                                                       true, Q2_GeV, v));
  }

  return resp;
}

systtools::event_unit_response_t
MiscInteractionSysts::GetEventResponse(genie::EventRecord const &ev) {

  systtools::event_unit_response_t resp;

  systtools::SystMetaData const &md = GetSystMetaData();

  if (pidx_C12ToAr40_2p2hScaling_nu != systtools::kParamUnhandled<size_t>) {
    std::vector<double> wght = (ev.Probe()->Pdg() > 0) ? GetWeights_C12ToAr40_2p2hScaling(ev, md[pidx_C12ToAr40_2p2hScaling_nu].paramVariations)
                                                       : std::vector<double>(md[pidx_C12ToAr40_2p2hScaling_nu].paramVariations.size(), 1.);
    if (wght.size()) {
      resp.push_back(
          {md[pidx_C12ToAr40_2p2hScaling_nu].systParamId, std::move(wght)});
    }
  }
  if (pidx_C12ToAr40_2p2hScaling_nubar != systtools::kParamUnhandled<size_t>) {
    std::vector<double> wght = (ev.Probe()->Pdg() < 0) ? GetWeights_C12ToAr40_2p2hScaling(ev, md[pidx_C12ToAr40_2p2hScaling_nubar].paramVariations)
                                                       : std::vector<double>(md[pidx_C12ToAr40_2p2hScaling_nubar].paramVariations.size(), 1.);
    if (wght.size()) {
      resp.push_back(
          {md[pidx_C12ToAr40_2p2hScaling_nubar].systParamId, std::move(wght)});
    }
  }
  if (pidx_nuenuebar_xsec_ratio != systtools::kParamUnhandled<size_t>) {
    std::vector<double> wght = GetWeights_nuenuebar_xsec_ratio(
        ev, md[pidx_nuenuebar_xsec_ratio].paramVariations);
    if (wght.size()) {
      resp.push_back(
          {md[pidx_nuenuebar_xsec_ratio].systParamId, std::move(wght)});
    }
  }
  if (pidx_nuenumu_xsec_ratio != systtools::kParamUnhandled<size_t>) {
    std::vector<double> wght = GetWeights_nuenumu_xsec_ratio(
        ev, md[pidx_nuenumu_xsec_ratio].paramVariations);
    if (wght.size()) {
      resp.push_back(
          {md[pidx_nuenumu_xsec_ratio].systParamId, std::move(wght)});
    }
  }
  if (pidx_SPPLowQ2Suppression != systtools::kParamUnhandled<size_t>) {
    std::vector<double> wght = GetWeights_SPPLowQ2Suppression(
        ev, md[pidx_SPPLowQ2Suppression].paramVariations);
    if (wght.size()) {
      resp.push_back(
          {md[pidx_SPPLowQ2Suppression].systParamId, std::move(wght)});
    }
  }

  // Remove any empty responses
  for (size_t r_it = 0; r_it < resp.size();) {
    if (!resp[r_it].responses.size()) {
      resp.erase(resp.begin() + r_it);
    } else {
      ++r_it;
    }
  }

  if (fill_valid_tree) {
    genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
    genie::GHepParticle *ISLep = ev.Probe();
    TLorentzVector FSLepP4 = *FSLep->P4();
    TLorentzVector ISLepP4 = *ISLep->P4();
    TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

    int Pdgnu = ISLep->Pdg();

    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();

    W = ev.Summary()->Kine().W(true);

    valid_tree->Fill();
  }

  return resp;
}
std::string MiscInteractionSysts::AsString() { return "MiscInteractionSysts"; }

void MiscInteractionSysts::InitValidTree() {
  valid_file = new TFile("MiscInteractionSysts_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
}

MiscInteractionSysts::~MiscInteractionSysts() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
