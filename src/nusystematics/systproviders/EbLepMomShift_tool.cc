#include "nusystematics/systproviders/EbLepMomShift_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"

using namespace nusyst;
using namespace systtools;

EbLepMomShift::EbLepMomShift(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData EbLepMomShift::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                              paramId_t firstId) {

  SystMetaData smd;

  SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(ps, "EbFSLepMomShift", phdr,
                                                 firstId)) {
    phdr.systParamId = firstId++;
    phdr.isWeightSystematicVariation = false;
    smd.push_back(phdr);
  }

  if (!ps.has_key("EbLepMomShift_Template_input_manifest") ||
      !ps.is_key_to_table("EbLepMomShift_Template_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for EbLepMomShift, "
           "expected to find a FHiCL table keyed by "
           "EbLepMomShift_Template_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }
  fhicl::ParameterSet templateManifest =
      ps.get<fhicl::ParameterSet>("EbLepMomShift_Template_input_manifest");
  tool_options.put("EbLepMomShift_Template_input_manifest", templateManifest);

  tool_options.put("fill_valid_tree", ps.get<bool>("fill_valid_tree", false));

  return smd;
}

bool EbLepMomShift::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  SystMetaData const &md = GetSystMetaData();

  if (!HasParam(md, "EbFSLepMomShift")) {
    throw incorrectly_configured()
        << "[ERROR]: Expected to find parameter named "
        << std::quoted("EbFSLepMomShift");
  }

  if (!tool_options.has_key("EbLepMomShift_Template_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: EbFSLepMomShift parameter exists in the SystMetaData, but "
           "no EbLepMomShift_Template_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "EbLepMomShift::BuildSystMetaData, but wasn't, this is a bug, "
           "please report to the maintiner.";
  }

  fhicl::ParameterSet const &templateManifest =
      tool_options.get<fhicl::ParameterSet>(
          "EbLepMomShift_Template_input_manifest");

  ResponseParameterIdx = GetParamIndex(md, "EbFSLepMomShift");

  EbTemplate.LoadInputHistograms(templateManifest);

  fill_valid_tree = tool_options.get("fill_valid_tree", false);

  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
EbLepMomShift::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC() || ev.Summary()->ExclTag().IsCharmEvent()) {
    return resp;
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();

  Enu = ISLepP4.E();
  FSLep_ctheta = FSLepP4.Vect().CosTheta();

  int bin = EbTemplate.GetBin({{Enu, FSLep_ctheta}});
  if (bin != kBinOutsideRange) {
    resp.push_back({md[ResponseParameterIdx].systParamId, {}});
    for (double v : md[ResponseParameterIdx].paramVariations) {
      if ((v == 0) && !EbTemplate.IsValidVariation(0)) {
        resp.back().responses.push_back(0);
      } else {
        resp.back().responses.push_back(
            EbTemplate.GetVariation(v, {{Enu, FSLep_ctheta}}));
      }
    }
  }

  if (fill_valid_tree) {

    int Pdgnu = ISLep->Pdg();

    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

    FSLep_pmu = FSLepP4.Vect().Mag();
    shift = resp.front().responses[3];

    BinOutsideRange = (bin == kBinOutsideRange);

    valid_tree->Fill();
  }

  return resp;
}
std::string EbLepMomShift::AsString() { return "EbLepMomShift"; }

void EbLepMomShift::InitValidTree() {
  valid_file = new TFile("EbLepMomShift_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("FSLep_pmu", &FSLep_pmu);
  valid_tree->Branch("FSLep_ctheta", &FSLep_ctheta);
  valid_tree->Branch("BinOutsideRange", &BinOutsideRange);
  valid_tree->Branch("shift", &shift);
}

EbLepMomShift::~EbLepMomShift() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
