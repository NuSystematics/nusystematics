#include "nusystematics/systproviders/MKSinglePiTemplate_tool.hh"

#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/exceptions.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "Framework/Messenger/Messenger.h"

using namespace systtools;
using namespace nusyst;
using namespace systtools;

// #define DEBUG_MKSINGLEPI

MKSinglePiTemplate::MKSinglePiTemplate(YAML::Node const &params)
    : IGENIESystProvider_tool(params),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

namespace {
struct channel_id {
  std::string name;
  genie::SppChannel_t channel;
};
} // namespace

SystMetaData MKSinglePiTemplate::BuildSystMetaData(YAML::Node const &cfg,
                                                   paramId_t firstId) {

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseYAMLToolConfigurationParameter(cfg, "MKSPP_ReWeight", phdr,
                                                 firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  if (!cfg["MKSPP_Template_input_manifest"] ||
      !cfg["MKSPP_Template_input_manifest"].IsMap()) {
    throw invalid_ToolConfigurationYAML()
        << "[ERROR]: When configuring calculated variations for "
           "MKSPP_Template, expected to find a YAML table keyed by "
           "MKSPP_Template_input_manifest describing the location of the "
           "histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh "
           "for the layout.";
  }
  YAML::Node templateManifest =
      cfg["MKSPP_Template_input_manifest"];
  tool_options["MKSPP_Template_input_manifest"] = templateManifest;

  size_t NNuChannels = 0;
  for (std::string const &ch : {"NumuPPiPlus", "NumuPPi0", "NumuNPiPlus"}) {

    if (!templateManifest[ch]) {
      continue;
    }

    NNuChannels++;
  }

  size_t NAntiNuChannels = 0;
  for (std::string const &ch :
       {"NumuBNPiMinus", "NumuBNPi0", "NumuBPPiMinus"}) {

    if (!templateManifest[ch]) {
      continue;
    }

    NAntiNuChannels++;
  }

  size_t NChannels = NNuChannels + NAntiNuChannels;
  if (!NChannels) {
    throw invalid_ToolConfigurationYAML()
        << "[ERROR]: When configuring a MKSPP_Template reweighting instance, "
           "failed to find any configured channels. Input templates must be "
           "described by in a YAML table keyed MKSPP_Template_input_manifest with "
           "the layout follows that consumed by "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh";
  }

  SuppressNeutrinoBkgSPP = NNuChannels;
  SuppressAntiNeutrinoBkgSPP = NAntiNuChannels;

  tool_options["SuppressNeutrinoBkgSPP"] = SuppressNeutrinoBkgSPP;
  tool_options["SuppressAntiNeutrinoBkgSPP"] = SuppressAntiNeutrinoBkgSPP;

  fill_valid_tree = cfg["fill_valid_tree"] ? cfg["fill_valid_tree"].as<bool>() : false;
  tool_options["fill_valid_tree"] = fill_valid_tree;

  use_Q2W_templates = cfg["use_Q2W_templates"] ? cfg["use_Q2W_templates"].as<bool>() : true;
  tool_options["use_Q2W_templates"] = use_Q2W_templates;

  Q2_or_q0_is_x = cfg["Q2_or_q0_is_x"] ? cfg["Q2_or_q0_is_x"].as<bool>() : false;
  tool_options["Q2_or_q0_is_x"] = Q2_or_q0_is_x;

  return smd;
}

bool MKSinglePiTemplate::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  if (!HasParam(GetSystMetaData(), "MKSPP_ReWeight")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("MKSPP_ReWeight");
  }

  if (!tool_options["MKSPP_Template_input_manifest"]) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: MKSPP_ReWeight parameter exists in the SystMetaData, but "
           "no MKSPP_Template_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "MKSinglePiTemplate::BuildSystMetaData, but wasn't, this is a bug, "
           "please report to the maintiner.";
  }

  YAML::Node const &templateManifest =
      tool_options["MKSPP_Template_input_manifest"];

  ResponseParameterIdx = GetParamIndex(GetSystMetaData(), "MKSPP_ReWeight");

  for (channel_id const &ch :
       std::vector<channel_id>{{"NumuPPiPlus", genie::kSpp_vp_cc_10100},
                               {"NumuPPi0", genie::kSpp_vn_cc_10010},
                               {"NumuNPiPlus", genie::kSpp_vn_cc_01100},

                               {"NumuBNPiMinus", genie::kSpp_vbn_cc_01001},
                               {"NumuBNPi0", genie::kSpp_vbp_cc_01010},
                               {"NumuBPPiMinus", genie::kSpp_vbp_cc_10001}}) {

    if (!templateManifest[ch.name]) {
      continue;
    }

    TemplateHelper th;
    th.Template = std::make_unique<MKSinglePiTemplate_ReWeight>(
        templateManifest[ch.name]);
    th.ZeroIsValid = th.Template->IsValidVariation(0);

    ChannelParameterMapping.emplace(ch.channel, std::move(th));
  }

  SuppressNeutrinoBkgSPP = tool_options["SuppressNeutrinoBkgSPP"] ? tool_options["SuppressNeutrinoBkgSPP"].as<bool>() : false;
  SuppressAntiNeutrinoBkgSPP =
      tool_options["SuppressAntiNeutrinoBkgSPP"] ? tool_options["SuppressAntiNeutrinoBkgSPP"].as<bool>() : false;

  fill_valid_tree = tool_options["fill_valid_tree"] ? tool_options["fill_valid_tree"].as<bool>() : false;
  use_Q2W_templates = tool_options["use_Q2W_templates"] ? tool_options["use_Q2W_templates"].as<bool>() : true;
  Q2_or_q0_is_x = tool_options["Q2_or_q0_is_x"] ? tool_options["Q2_or_q0_is_x"].as<bool>() : true;

  if (fill_valid_tree) {
    InitValidTree();
  }

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");
  genie::Messenger::Instance()->SetPriorityLevel("GHepUtils",
                                                 log4cpp::Priority::FATAL);

  return true;
}

event_unit_response_t
MKSinglePiTemplate::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return resp;
  }

  bool is_res = ev.Summary()->ProcInfo().IsResonant();

  if (!(is_res || ev.Summary()->ProcInfo().IsDeepInelastic())) {
    return resp;
  }

  if (ev.Summary()->Kine().W(true) > 1.7) {
    return resp;
  }

  bool is_nu = (ev.Probe()->Pdg() > 0);

  // Only suppress non-resonant background channels when we have templates to
  // reweight to.
  if (!is_res && ((is_nu && !SuppressNeutrinoBkgSPP) ||
                  (!is_nu && !SuppressAntiNeutrinoBkgSPP))) {
    return resp;
  }

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  genie::SppChannel_t chan = genie::kSppNull;

  genie::Target const &tgt = ev.Summary()->InitState().Tgt();
  if (!tgt.HitNucIsSet()) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to get hit nucleon kinematics as it was not "
           "included in this GHep event. This is a fatal error.";
  }

  TLorentzVector NucP4 = tgt.HitNucP4();

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();

  FSLepP4.Boost(-NucP4.BoostVector());
  ISLepP4.Boost(-NucP4.BoostVector());
  NucP4.Boost(-NucP4.BoostVector());

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  if (ev.Summary()->ProcInfo().IsResonant()) {

    chan = SPPChannelFromGHep(ev);

    if ((chan == genie::kSppNull) ||
        (ChannelParameterMapping.find(chan) == ChannelParameterMapping.end())) {

#ifdef DEBUG_MKSINGLEPI
      int neut_code = abs(genie::utils::ghep::NeutReactionCode(&ev));
      if ((neut_code == 21)    // Multi-pi
          || (neut_code == 22) // Eta-production
          || (neut_code == 23) // Kaon-production
          || (neut_code == 17) // gamma-production
      ) {
        return resp;
      }
      std::cout << "[INFO]: RES event failed to find SPP Channel  (NEUT: "
                << genie::utils::ghep::NeutReactionCode(&ev) << ") "
                << std::endl
                << DumpGENIEEv(ev) << std::endl;
#endif
      return resp;
    }

    std::array<double, 2> kinematics;
    kinematics[0] = use_Q2W_templates ? -emTransfer.Mag2() : emTransfer.E();
    kinematics[1] = use_Q2W_templates ? ev.Summary()->Kine().W(true)
                                      : emTransfer.Vect().Mag();

    if (Q2_or_q0_is_x) {
      std::swap(kinematics[0], kinematics[1]);
    }

    resp.push_back({hdr.systParamId, {}});
    for (double val : hdr.paramVariations) {

      if ((val == 0) && !ChannelParameterMapping[chan].ZeroIsValid) {
        resp.back().responses.push_back(1);
      } else {
        resp.back().responses.push_back(
            ChannelParameterMapping[chan].Template->GetVariation(
                val, ISLepP4.E(), kinematics));
      }
    }
  } else { // Non-resonant background has to die off as MK is turned on, as the
           // MK prediction includes the coupled background channels
    resp.push_back({hdr.systParamId, {}});
    for (double val : hdr.paramVariations) {
      val = std::min(fabs(val), 1.0);
      resp.back().responses.push_back(1 - val);
    }
  }

  if (fill_valid_tree) {

    q0_nuc_rest_frame = emTransfer.E();
    q3_nuc_rest_frame = emTransfer.Vect().Mag();
    Enu_nuc_rest_frame = ISLepP4.E();

    ISLepP4 = *ISLep->P4();
    FSLepP4 = *FSLep->P4();
    emTransfer = (ISLepP4 - FSLepP4);

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
    IsDIS = ev.Summary()->ProcInfo().IsDeepInelastic();

    SppChannel = chan;

    momhmfspi = 0;
    size_t NParts = ev.GetEntries();
    for (size_t it = 0; it < NParts; ++it) {
      genie::GHepParticle *part = ev.Particle(it);
      if (part->Status() != genie::kIStStableFinalState) {
        continue;
      }
      if ((part->Pdg() != 211) && (part->Pdg() != 111)) {
        continue;
      }
      if (momhmfspi < part->P4()->Vect().Mag()) {
        pdghmfspi = part->Pdg();
        momhmfspi = part->P4()->Vect().Mag();
        cthetahmfspi = part->P4()->Vect().CosTheta();
      }
    }

    weight = resp.back().responses.back();

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    std::array<double, 2> kinematics;
    kinematics[0] = use_Q2W_templates ? -emTransfer.Mag2() : emTransfer.E();
    kinematics[1] = use_Q2W_templates ? ev.Summary()->Kine().W(true)
                                      : emTransfer.Vect().Mag();

    valid_tree->Fill();
  }

  return resp;
}

std::string MKSinglePiTemplate::AsString() { return ""; }

void MKSinglePiTemplate::InitValidTree() {
  valid_file = new TFile("MKSinglePiTemplate_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("SppChannel", &SppChannel);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Pdg_HMFSPi", &pdghmfspi);
  valid_tree->Branch("P_HMFSPi", &momhmfspi);
  valid_tree->Branch("CosTheta_HMFSPi", &cthetahmfspi);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("Enu_nuc_rest_frame", &Enu_nuc_rest_frame);
  valid_tree->Branch("q0_nuc_rest_frame", &q0_nuc_rest_frame);
  valid_tree->Branch("q3_nuc_rest_frame", &q3_nuc_rest_frame);
  valid_tree->Branch("weight", &weight);
  valid_tree->Branch("IsDIS", &IsDIS);
}

MKSinglePiTemplate::~MKSinglePiTemplate() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
