#include "nusystematics/systproviders/FSILikeEAvailSmearing_tool.hh"

#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/exceptions.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "Framework/Messenger/Messenger.h"

using namespace systtools;
using namespace nusyst;

// #define DEBUG_MKSINGLEPI

FSILikeEAvailSmearing::FSILikeEAvailSmearing(YAML::Node const &params)
    : IGENIESystProvider_tool(params),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>) {}

SystMetaData FSILikeEAvailSmearing::BuildSystMetaData(YAML::Node const &cfg,
                                                      paramId_t firstId) {

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseYAMLToolConfigurationParameter(cfg, "FSILikeEAvailSmearing",
                                          phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  YAML::Node templateManifest =
      cfg["FSILikeEAvailSmearing_input_manifest"];

  if (!cfg["FSILikeEAvailSmearing_input_manifest"] ||
      !cfg["FSILikeEAvailSmearing_input_manifest"].IsMap()) {
    throw invalid_ToolConfigurationYAML()
        << "[ERROR]: When configuring calculated variations for "
           "FSILikeEAvailSmearing, expected to find a YAML table keyed by "
           "FSILikeEAvailSmearing_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }

  size_t NChannels = 0;
  for (std::string const &ch : {"CCQE", "CCRes", "CCDIS", "CCMEC", "CCQE_bar",
                                "CCRes_bar", "CCDIS_bar", "CCMEC_bar", "NC"}) {

    if (!templateManifest[ch]) {
      continue;
    }

    NChannels++;
  }

  if (!NChannels) {
    throw invalid_ToolConfigurationYAML()
        << "[ERROR]: When configuring a FSILikeEAvailSmearing reweighting "
           "instance, failed to find any configured channels. Input templates "
           "must be described by in a table keyed "
           "FSILikeEAvailSmearing_input_manifest with the layout follows that "
           "consumed by "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh";
  }

  tool_options["FSILikeEAvailSmearing_input_manifest"] = templateManifest;

  LimitWeights = cfg["LimitWeights"] ?
    cfg["LimitWeights"].as<std::pair<double, double>>() :
    std::pair<double, double>{0, std::numeric_limits<double>::max()};
  tool_options["LimitWeights"] = YAML::Node();
  tool_options["LimitWeights"].push_back(LimitWeights.first);
  tool_options["LimitWeights"].push_back(LimitWeights.second);

  return smd;
}

namespace {
struct channel_id {
  std::string name;
  FSILikeEAvailSmearing::chan channel;
};
} // namespace

bool FSILikeEAvailSmearing::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  if (!HasParam(GetSystMetaData(), "FSILikeEAvailSmearing")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("FSILikeEAvailSmearing");
  }

  if (!tool_options["FSILikeEAvailSmearing_input_manifest"]) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: FSILikeEAvailSmearing parameter exists in the "
           "SystMetaData, but "
           "no FSILikeEAvailSmearing_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "FSILikeEAvailSmearing::BuildSystMetaData, but wasn't, this is a "
           "bug, "
           "please report to the maintiner.";
  }

  YAML::Node const &templateManifest =
      tool_options["FSILikeEAvailSmearing_input_manifest"];

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "FSILikeEAvailSmearing");

  for (channel_id const &ch :
       std::vector<channel_id>{{"CCQE", chan::kCCQE},
                               {"CCRes", chan::kCCRes},
                               {"CCDIS", chan::kCCDIS},
                               {"CCMEC", chan::kCCMEC},
                               {"CCQE_bar", chan::kCCQE_bar},
                               {"CCRes_bar", chan::kCCRes_bar},
                               {"CCDIS_bar", chan::kCCDIS_bar},
                               {"CCMEC_bar", chan::kCCMEC_bar},
                               {"NC", chan::kNC}}) {

    if (!templateManifest[ch.name]) {
      continue;
    }

    TemplateHelper th;
    th.Template = std::make_unique<FSILikeEAvailSmearing_ReWeight>();
    th.Template->LoadInputHistograms(
        templateManifest[ch.name]);
    th.ZeroIsValid = th.Template->IsValidVariation(0);

    ChannelParameterMapping.emplace(ch.channel, std::move(th));
  }

  LimitWeights = tool_options["LimitWeights"] ?
    std::pair<double, double>{tool_options["LimitWeights"][0].as<double>(), tool_options["LimitWeights"][1].as<double>()} :
    std::pair<double, double>{0, std::numeric_limits<double>::max()};

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");
  genie::Messenger::Instance()->SetPriorityLevel("GHepUtils",
                                                 log4cpp::Priority::FATAL);

  return true;
}

FSILikeEAvailSmearing::chan GetChan(simb_mode_copy mode, bool IsCC, bool IsNu) {

  if (!IsCC) {
    return FSILikeEAvailSmearing::chan::kNC;
  }

  switch (mode) {
  case simb_mode_copy::kQE: {
    return IsNu ? FSILikeEAvailSmearing::chan::kCCQE
                : FSILikeEAvailSmearing::chan::kCCQE_bar;
  }
  case simb_mode_copy::kRes: {
    return IsNu ? FSILikeEAvailSmearing::chan::kCCRes
                : FSILikeEAvailSmearing::chan::kCCRes_bar;
  }
  case simb_mode_copy::kDIS: {
    return IsNu ? FSILikeEAvailSmearing::chan::kCCDIS
                : FSILikeEAvailSmearing::chan::kCCDIS_bar;
  }
  case simb_mode_copy::kMEC: {
    return IsNu ? FSILikeEAvailSmearing::chan::kCCMEC
                : FSILikeEAvailSmearing::chan::kCCMEC_bar;
  }
  default: { return FSILikeEAvailSmearing::chan::kBadChan; }
  }
}

event_unit_response_t
FSILikeEAvailSmearing::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;

  // Ignore Coherent
  simb_mode_copy mode = GetSimbMode(ev);
  if (mode == simb_mode_copy::kCoh) {
    return resp;
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  chan evch =
      GetChan(mode, ev.Summary()->ProcInfo().IsWeakCC(), ISLep->Pdg() > 0);

  if (ChannelParameterMapping.find(evch) == ChannelParameterMapping.end()) {
    return resp;
  }

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();

  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  std::array<double, 3> kinematics;
  kinematics[0] = emTransfer.Vect().Mag();
  kinematics[1] = emTransfer[3];
  kinematics[2] = GetErecoil_MINERvA_LowRecoil(ev) / kinematics[1];

  resp.push_back({hdr.systParamId, {}});
  for (double val : hdr.paramVariations) {

    if ((val == 0) && !ChannelParameterMapping[evch].ZeroIsValid) {
      resp.back().responses.push_back(1);
    } else {
      double wght =
          ChannelParameterMapping[evch].Template->GetVariation(val, kinematics);

      wght = (wght < LimitWeights.first) ? LimitWeights.first : wght;
      wght = (wght > LimitWeights.second) ? LimitWeights.second : wght;

      resp.back().responses.push_back(wght);
    }
  }
  return resp;
}

std::string FSILikeEAvailSmearing::AsString() { return ""; }

FSILikeEAvailSmearing::~FSILikeEAvailSmearing() {}
