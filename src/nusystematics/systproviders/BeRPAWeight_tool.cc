#include "nusystematics/systproviders/BeRPAWeight_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "systematicstools/utility/ResponselessParamUtility.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/responsecalculators/BeRPA.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"

using namespace nusyst;
using namespace systtools;

// #define BERPAWEIGHT_DEBUG

BeRPAWeight::BeRPAWeight(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_BeRPA_Response(kParamUnhandled<size_t>),
      pidx_BeRPA_A(kParamUnhandled<size_t>),
      pidx_BeRPA_B(kParamUnhandled<size_t>),
      pidx_BeRPA_D(kParamUnhandled<size_t>),
      pidx_BeRPA_E(kParamUnhandled<size_t>), valid_file(nullptr),
      valid_tree(nullptr) {}

SystMetaData BeRPAWeight::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                            paramId_t firstId) {

  SystMetaData smd;

  ignore_parameter_dependence =
      ps.get<bool>("ignore_parameter_dependence", false);
  ApplyCV = ps.get<bool>("ApplyCV", false);

  SystParamHeader responseParam;
  std::vector<std::string> dependentParamNames;
  if (!ignore_parameter_dependence) {
    responseParam.prettyName = "BeRPA_Response";
    responseParam.systParamId = firstId++;
  }

  for (std::string const &pname :
       {"BeRPA_A", "BeRPA_B", "BeRPA_D", "BeRPA_E"}) {
    SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(ps, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      if (!ignore_parameter_dependence) {
        dependentParamNames.push_back(phdr.prettyName);
        phdr.isResponselessParam = true;
        phdr.responseParamId = responseParam.systParamId;
        if (phdr.isSplineable) {
          throw invalid_ToolConfigurationFHiCL()
              << "[ERROR]: Attempted to build spline from "
                 "parameter "
              << phdr.prettyName
              << ", which enters into an intrinsically "
                 "multi-parameter response calculation. Either run in "
                 "random "
                 "throw mode or set \"ignore_parameter_dependence\" "
                 "in the "
                 "BeRPAWeight configuration.";
        }
      }
      smd.push_back(phdr);
    }
  }

  if (!ignore_parameter_dependence) {
    smd.push_back(responseParam);
    FinalizeAndValidateDependentParameters(smd, smd.back().prettyName,
                                           dependentParamNames);
  }

  tool_options.put("fill_valid_tree", ps.get<bool>("fill_valid_tree", false));
  tool_options.put("ignore_parameter_dependence", ignore_parameter_dependence);
  tool_options.put("ApplyCV", ApplyCV);

  return smd;
}

bool BeRPAWeight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  ignore_parameter_dependence =
      tool_options.get<bool>("ignore_parameter_dependence", false);

  SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "BeRPA_Response")) {
    pidx_BeRPA_Response = GetParamIndex(md, "BeRPA_Response");
  }

  ACV = 0;
  if (HasParam(md, "BeRPA_A")) {
    pidx_BeRPA_A = GetParamIndex(md, "BeRPA_A");
    ACV = md[pidx_BeRPA_A].centralParamValue;
    if (ACV == kDefaultDouble) {
      ACV = 0;
    }
    std::copy(md[pidx_BeRPA_A].paramVariations.begin(),
              md[pidx_BeRPA_A].paramVariations.end(),
              std::back_inserter(AVariations));
  }

  BCV = 0;
  if (HasParam(md, "BeRPA_B")) {
    pidx_BeRPA_B = GetParamIndex(md, "BeRPA_B");
    BCV = md[pidx_BeRPA_B].centralParamValue;
    if (BCV == kDefaultDouble) {
      BCV = 0;
    }
    std::copy(md[pidx_BeRPA_B].paramVariations.begin(),
              md[pidx_BeRPA_B].paramVariations.end(),
              std::back_inserter(BVariations));
  }

  DCV = 0;
  if (HasParam(md, "BeRPA_D")) {
    pidx_BeRPA_D = GetParamIndex(md, "BeRPA_D");
    DCV = md[pidx_BeRPA_D].centralParamValue;
    if (DCV == kDefaultDouble) {
      DCV = 0;
    }
    std::copy(md[pidx_BeRPA_D].paramVariations.begin(),
              md[pidx_BeRPA_D].paramVariations.end(),
              std::back_inserter(DVariations));
  }

  ECV = 0;
  if (HasParam(md, "BeRPA_E")) {
    pidx_BeRPA_E = GetParamIndex(md, "BeRPA_E");
    ECV = md[pidx_BeRPA_E].centralParamValue;
    if (ECV == kDefaultDouble) {
      ECV = 0;
    }
    std::copy(md[pidx_BeRPA_E].paramVariations.begin(),
              md[pidx_BeRPA_E].paramVariations.end(),
              std::back_inserter(EVariations));
  }

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  ApplyCV = tool_options.get<bool>("ApplyCV", false);

  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
BeRPAWeight::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC() ||
      ev.Summary()->ExclTag().IsCharmEvent()) {
    return resp;
  }

#ifdef BERPAWEIGHT_DEBUG
  if (ev.Summary()->Kine().W(true) > 1) {
    std::cout << "[INFO]: QE event with high W (NEUT: "
              << genie::utils::ghep::NeutReactionCode(&ev) << ") " << std::endl
              << DumpGENIEEv(ev) << std::endl;
    return resp;
  }
#endif

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
  Q2 = -emTransfer.Mag2();

  // Only want the CV response to be used in one of the dials, after the first
  // dial is found, all other dial responses should be /= CVResponse.
  double CVResponse =
      GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, ACV, BCV, DCV, ECV);

#ifdef BERPAWEIGHT_DEBUG
  std::cout << "[CV Response @ " << ACV << ", " << BCV << ", " << DCV << ", "
            << ECV << "] = " << CVResponse << std::endl;
#endif

  if (!ignore_parameter_dependence) {
    resp.push_back({md[pidx_BeRPA_Response].systParamId, {}});

    for (size_t univ = 0; univ < md[pidx_BeRPA_Response].paramVariations.size();
         ++univ) {

      double Aval = AVariations[univ];
      double Bval = BVariations[univ];
      double Dval = DVariations[univ];
      double Eval = EVariations[univ];

      double weight = GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, Aval,
                                     Bval, Dval, Eval);
      if (!ApplyCV) {
        weight /= CVResponse;
      }

      resp.back().responses.push_back(weight);
    }
  } else {

    bool UsedADial = false;
    if (pidx_BeRPA_A != kParamUnhandled<size_t>) {
      resp.push_back({md[pidx_BeRPA_A].systParamId, {}});
      for (double av : AVariations) {
        double weight = GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, av,
                                       BCV, DCV, ECV);
#ifdef BERPAWEIGHT_DEBUG
        std::cout << "[ weight @ " << av << ", " << BCV << ", " << DCV << ", "
                  << ECV << "] = " << weight << std::endl;
#endif
        if (!ApplyCV) {
          weight /= CVResponse;
        }
        resp.back().responses.push_back(weight);
      }
      UsedADial = true;
    }
    if (pidx_BeRPA_B != kParamUnhandled<size_t>) {
      resp.push_back({md[pidx_BeRPA_B].systParamId, {}});
      for (double bv : BVariations) {
        double weight = GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, ACV,
                                       bv, DCV, ECV);
#ifdef BERPAWEIGHT_DEBUG
        std::cout << "[ weight @ " << ACV << ", " << bv << ", " << DCV << ", "
                  << ECV << "] = " << weight << std::endl;
#endif
        if (!ApplyCV || UsedADial) {
          weight /= CVResponse;
        }
        resp.back().responses.push_back(weight);
      }
      UsedADial = true;
    }
    if (pidx_BeRPA_D != kParamUnhandled<size_t>) {
      resp.push_back({md[pidx_BeRPA_D].systParamId, {}});
      for (double dv : DVariations) {
        double weight = GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, ACV,
                                       BCV, dv, ECV);
#ifdef BERPAWEIGHT_DEBUG
        std::cout << "[ weight @ " << ACV << ", " << BCV << ", " << dv << ", "
                  << ECV << "] = " << weight << std::endl;
#endif
        if (!ApplyCV || UsedADial) {
          weight /= CVResponse;
        }
        resp.back().responses.push_back(weight);
      }
      UsedADial = true;
    }
    if (pidx_BeRPA_E != kParamUnhandled<size_t>) {
      resp.push_back({md[pidx_BeRPA_E].systParamId, {}});
      for (double eval : EVariations) {
        double weight = GetBeRPAWeight(e2i(simb_mode_copy::kQE), true, Q2, ACV,
                                       BCV, DCV, eval);
#ifdef BERPAWEIGHT_DEBUG
        std::cout << "[ weight @ " << ACV << ", " << BCV << ", " << DCV << ", "
                  << eval << "] = " << weight << std::endl;
#endif
        if (!ApplyCV || UsedADial) {
          weight /= CVResponse;
        }
        resp.back().responses.push_back(weight);
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

    Enu = ISLepP4.E();
    weight = 1;

    for (auto const &eur : resp) {
      weight *= eur.responses[3];
    }

    valid_tree->Fill();
  }

  return resp;
}
std::string BeRPAWeight::AsString() { return "BeRPAWeight"; }

void BeRPAWeight::InitValidTree() {
  valid_file = new TFile("BeRPAWeight_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("weight", &weight);
}

BeRPAWeight::~BeRPAWeight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
