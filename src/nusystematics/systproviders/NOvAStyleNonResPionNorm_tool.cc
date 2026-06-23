#include "nusystematics/systproviders/NOvAStyleNonResPionNorm_tool.hh"

#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

using namespace nusyst;

NOvAStyleNonResPionNorm::NOvAStyleNonResPionNorm(
    YAML::Node const &params)
    : IGENIESystProvider_tool(params), valid_file(nullptr),
      valid_tree(nullptr) {}

systtools::SystMetaData
NOvAStyleNonResPionNorm::BuildSystMetaData(YAML::Node const &yamlnd,
                                           systtools::paramId_t firstId) {

  systtools::SystMetaData smd;

  for (bool IsNeutrino : {true, false}) {
    for (bool IsCC : {true, false}) {
      for (int TargetNucleon : {1, 2, 3}) { // neutron, proton, either
        for (size_t NPi : {1, 2, 3}) {
          NRPiChan_t chan =
              BuildNRPiChannel(IsNeutrino, IsCC, TargetNucleon, NPi);

          systtools::SystParamHeader ch_param;
            if (ParseYAMLToolConfigurationParameter(
              yamlnd, GetNRPiChannelName(chan), ch_param, firstId)) {
            ch_param.systParamId = firstId++;
            smd.push_back(ch_param);
          }
        }
      }
    }
  }

  WBegin = yamlnd["WBegin"] ? yamlnd["WBegin"].as<double>() : 0;
  WEnd = yamlnd["WEnd"] ? yamlnd["WEnd"].as<double>() : 5;
  WTransition = yamlnd["WTransition"] ? yamlnd["WTransition"].as<double>() : 3;
  OneSigmaResponse = yamlnd["OneSigmaResponse"] ? yamlnd["OneSigmaResponse"].as<double>() : 0.5;
  HighWResponse = yamlnd["HighWResponse"] ? yamlnd["HighWResponse"].as<double>() : 0.05;
  fill_valid_tree = yamlnd["fill_valid_tree"] ? yamlnd["fill_valid_tree"].as<bool>() : false;

  tool_options["WBegin"] = WBegin;
  tool_options["WEnd"] = WEnd;
  tool_options["WTransition"] = WTransition;
  tool_options["OneSigmaResponse"] = OneSigmaResponse;
  tool_options["HighWResponse"] = HighWResponse;
  tool_options["fill_valid_tree"] = fill_valid_tree;

  return smd;
}

bool NOvAStyleNonResPionNorm::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  systtools::SystMetaData const &md = GetSystMetaData();

  for (bool IsNeutrino : {true, false}) {
    for (bool IsCC : {true, false}) {
      for (int TargetNucleon : {1, 2, 3}) { // neutron, proton, either
        for (size_t NPi : {1, 2, 3}) {
          NRPiChan_t chan =
              BuildNRPiChannel(IsNeutrino, IsCC, TargetNucleon, NPi);

          std::string ch_name = GetNRPiChannelName(chan);
          if (!HasParam(md, ch_name)) {
            continue;
          }
          ChannelParameterMapping.push_back({chan, GetParamIndex(md, ch_name)});
        }
      }
    }
  }

  WBegin = tool_options["WBegin"] ? tool_options["WBegin"].as<double>() : 0;
  WEnd = tool_options["WEnd"] ? tool_options["WEnd"].as<double>() : 5;
  WTransition = tool_options["WTransition"] ? tool_options["WTransition"].as<double>() : 3;
  OneSigmaResponse = tool_options["OneSigmaResponse"] ? tool_options["OneSigmaResponse"].as<double>() : 0.5;
  HighWResponse = tool_options["HighWResponse"] ? tool_options["HighWResponse"].as<double>() : 0.05;
  fill_valid_tree = tool_options["fill_valid_tree"] ? tool_options["fill_valid_tree"].as<bool>() : false;

  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

systtools::event_unit_response_t
NOvAStyleNonResPionNorm::GetEventResponse(genie::EventRecord const &ev) {

  systtools::event_unit_response_t resp = this->GetDefaultEventResponse();

  if (!ev.Summary()->ProcInfo().IsDeepInelastic()) {
    return resp;
  }
  double WTrue = ev.Summary()->Kine().W(true);
  if (WTrue < WBegin) {
    return resp;
  }

  NRPiChan_t chan = GetNRPiChannel(ev);

  if (!chan) {
    return resp;
  }

  double OneSigResp = OneSigmaResponse;

  if (WTrue > WEnd) {
    OneSigResp = HighWResponse;
  } else if (WTrue > WTransition) {
    OneSigResp -= (OneSigmaResponse - HighWResponse) *
                  ((WTrue - WTransition) / (WEnd - WTransition));
  }

  NRPiChan_t param_channel;
  systtools::SystParamHeader const *hdr = nullptr;
  size_t smdInx(-1);
  for (channel_param const &chpar : ChannelParameterMapping) {
    if (ChannelsAreEquivalent(chpar.channel, chan, 3)) {
      hdr = &GetSystMetaData()[chpar.paramidx];
      param_channel = chpar.channel;
      smdInx = chpar.paramidx;
      break;
    }
  }

  if (!hdr) {
    return resp;
  }

  resp[smdInx].responses.clear();
  for (double vals : hdr->paramVariations) {
    resp[smdInx].responses.push_back( std::max(0., 1 + vals * OneSigResp) );
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
    W = WTrue;

    NRPiChannel = chan;
    NRPiChannel_param = param_channel;
    NPi = GetNRPiChanNPi(chan);

    if (hdr->paramVariations.size() == 7) { // standard +/- 3 sigma
      weight_m1 = 1 + hdr->paramVariations[2]*OneSigResp;
      weight_p1 = 1 + hdr->paramVariations[4]*OneSigResp;
    }
    valid_tree->Fill();
  }

  return resp;
}
std::string NOvAStyleNonResPionNorm::AsString() {
  return "NOvAStyleNonResPionNorm";
}

void NOvAStyleNonResPionNorm::InitValidTree() {
  valid_file = new TFile("NOvAStyleNonResPionNorm_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("NRPiChannel", &NRPiChannel);
  valid_tree->Branch("NRPiChannel_param", &NRPiChannel_param);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("NPi", &NPi);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("weight_m1", &weight_m1);
  valid_tree->Branch("weight_p1", &weight_p1);
}

NOvAStyleNonResPionNorm::~NOvAStyleNonResPionNorm() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
