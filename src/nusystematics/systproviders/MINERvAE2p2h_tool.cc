#include "nusystematics/systproviders/MINERvAE2p2h_tool.hh"

#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"
#include "systematicstools/utility/ResponselessParamUtility.hh"

#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/responsecalculators/MINERvA2p2hEnergyDependenceScale.hh"

using namespace nusyst;
using namespace systtools;

// #define MINERVAE2p2h_DEBUG

MINERvAE2p2h::MINERvAE2p2h(YAML::Node const &params)
    : IGENIESystProvider_tool(params),
      pidx_E2p2hResponse_nu(kParamUnhandled<size_t>),
      pidx_E2p2hResponse_nubar(kParamUnhandled<size_t>),
      pidx_E2p2hA_nu(kParamUnhandled<size_t>),
      pidx_E2p2hB_nu(kParamUnhandled<size_t>),
      pidx_E2p2hA_nubar(kParamUnhandled<size_t>),
      pidx_E2p2hB_nubar(kParamUnhandled<size_t>), valid_file(nullptr),
      valid_tree(nullptr) {}

SystMetaData MINERvAE2p2h::BuildSystMetaData(YAML::Node const &yamlnd,
                                             paramId_t firstId) {

  SystMetaData smd;

  ignore_parameter_dependence =
      yamlnd["ignore_parameter_dependence"] ? yamlnd["ignore_parameter_dependence"].as<bool>() : false;

  for (std::string const &nu : {"nu", "nubar"}) {
    SystParamHeader responseParam;
    std::vector<std::string> dependentParamNames;
    if (!ignore_parameter_dependence) {
      responseParam.prettyName = std::string("E2p2hResponse_") + nu;
      responseParam.systParamId = firstId++;
    }

    for (std::string const &p : {"E2p2h_A", "E2p2h_B"}) {
      std::string pname = p + "_" + nu;
      SystParamHeader phdr;
      if (ParseYAMLToolConfigurationParameter(yamlnd, pname, phdr,
                      firstId)) {
        phdr.systParamId = firstId++;
        if (!ignore_parameter_dependence) {
          dependentParamNames.push_back(phdr.prettyName);
          phdr.isResponselessParam = true;
          phdr.responseParamId = responseParam.systParamId;
          if (phdr.isSplineable) {
            throw invalid_ToolConfigurationYAML()
                << "[ERROR]: Attempted to build spline from "
                   "parameter "
                << phdr.prettyName
                << ", which enters into an intrinsically "
                   "multi-parameter response calculation. Either run in "
                   "random "
                   "throw mode or set \"ignore_parameter_dependence\" "
                   "in the "
                   "MINERvAE2p2h configuration.";
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
  }

  tool_options["fill_valid_tree"] = yamlnd["fill_valid_tree"] ? yamlnd["fill_valid_tree"].as<bool>() : false;
  tool_options["ignore_parameter_dependence"] = ignore_parameter_dependence;

  LimitWeights = yamlnd["LimitWeights"] ? yamlnd["LimitWeights"].as<std::pair<double, double>>() : std::pair<double, double>{0, std::numeric_limits<double>::max()};

  tool_options["LimitWeights"] = YAML::Node();
  tool_options["LimitWeights"].push_back(LimitWeights.first);
  tool_options["LimitWeights"].push_back(LimitWeights.second);

  return smd;
}

bool MINERvAE2p2h::SetupResponseCalculator(
    YAML::Node const &tool_options) {

  SystMetaData const &md = GetSystMetaData();

  ignore_parameter_dependence =
      tool_options["ignore_parameter_dependence"] ? tool_options["ignore_parameter_dependence"].as<bool>() : false;

  if (HasParam(md, "E2p2hResponse_nu")) {
    pidx_E2p2hResponse_nu = GetParamIndex(md, "E2p2hResponse_nu");
  }

  A_nu_CV = 0;
  if (HasParam(md, "E2p2h_A_nu")) {
    pidx_E2p2hA_nu = GetParamIndex(md, "E2p2h_A_nu");
    A_nu_CV = md[pidx_E2p2hA_nu].centralParamValue;
    if (A_nu_CV == kDefaultDouble) {
      A_nu_CV = 0;
    }
    std::copy(md[pidx_E2p2hA_nu].paramVariations.begin(),
              md[pidx_E2p2hA_nu].paramVariations.end(),
              std::back_inserter(A_nu_Variations));
  }
  B_nu_CV = 0;
  if (HasParam(md, "E2p2h_B_nu")) {
    pidx_E2p2hB_nu = GetParamIndex(md, "E2p2h_B_nu");
    B_nu_CV = md[pidx_E2p2hB_nu].centralParamValue;
    if (B_nu_CV == kDefaultDouble) {
      B_nu_CV = 0;
    }
    std::copy(md[pidx_E2p2hB_nu].paramVariations.begin(),
              md[pidx_E2p2hB_nu].paramVariations.end(),
              std::back_inserter(B_nu_Variations));
  }

  if (HasParam(md, "E2p2hResponse_nubar")) {
    pidx_E2p2hResponse_nubar = GetParamIndex(md, "E2p2hResponse_nubar");
  }
  A_nubar_CV = 0;
  if (HasParam(md, "E2p2h_A_nubar")) {
    pidx_E2p2hA_nubar = GetParamIndex(md, "E2p2h_A_nubar");
    A_nubar_CV = md[pidx_E2p2hA_nubar].centralParamValue;
    if (A_nubar_CV == kDefaultDouble) {
      A_nubar_CV = 0;
    }
    std::copy(md[pidx_E2p2hA_nubar].paramVariations.begin(),
              md[pidx_E2p2hA_nubar].paramVariations.end(),
              std::back_inserter(A_nubar_Variations));
  }
  B_nubar_CV = 0;
  if (HasParam(md, "E2p2h_B_nubar")) {
    pidx_E2p2hB_nubar = GetParamIndex(md, "E2p2h_B_nubar");
    B_nubar_CV = md[pidx_E2p2hB_nubar].centralParamValue;
    if (B_nubar_CV == kDefaultDouble) {
      B_nubar_CV = 0;
    }
    std::copy(md[pidx_E2p2hB_nubar].paramVariations.begin(),
              md[pidx_E2p2hB_nubar].paramVariations.end(),
              std::back_inserter(B_nubar_Variations));
  }

  fill_valid_tree = tool_options["fill_valid_tree"] ? tool_options["fill_valid_tree"].as<bool>() : false;

  LimitWeights = tool_options["LimitWeights"] ? 
    std::pair<double, double>{tool_options["LimitWeights"][0].as<double>(), tool_options["LimitWeights"][1].as<double>()} :
    std::pair<double, double>{0, std::numeric_limits<double>::max()};

  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
MINERvAE2p2h::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  if (!ev.Summary()->ProcInfo().IsMEC() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return this->GetDefaultEventResponse();
  }

  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector ISLepP4 = *ISLep->P4();

  size_t pidx_Response, pidx_A, pidx_B;
  std::vector<double> *A_var, *B_var;
  double ACV, BCV;
  Enu = ISLepP4.E();

  for (int const &nu_pdgsign : {+1, -1}) {

    pidx_Response = nu_pdgsign>0 ? pidx_E2p2hResponse_nu : pidx_E2p2hResponse_nubar;
    pidx_A        = nu_pdgsign>0 ? pidx_E2p2hA_nu        : pidx_E2p2hA_nubar;
    pidx_B        = nu_pdgsign>0 ? pidx_E2p2hB_nu        : pidx_E2p2hB_nubar;
    A_var         = nu_pdgsign>0 ? &A_nu_Variations      : &A_nubar_Variations;
    B_var         = nu_pdgsign>0 ? &B_nu_Variations      : &B_nubar_Variations;
    ACV           = nu_pdgsign>0 ? A_nu_CV               : A_nubar_CV;
    BCV           = nu_pdgsign>0 ? B_nu_CV               : B_nubar_CV;

    bool nuMatched = (ISLep->Pdg() * nu_pdgsign > 0);

    if (!ignore_parameter_dependence) {

      resp.push_back({md[pidx_Response].systParamId, {}});

      for (size_t univ = 0; univ < md[pidx_Response].paramVariations.size();
           ++univ) {

        if(!nuMatched){
          resp.back().responses.push_back(1.);
          continue;
        }

        double Aval = A_var->at(univ);
        double Bval = B_var->at(univ);

        double weight = Get_MINERvA2p2h2EnergyDependencyScaling(
            e2i(simb_mode_copy::kMEC), true, Enu, Aval, Bval);

        weight = (weight < LimitWeights.first) ? LimitWeights.first : weight;
        weight = (weight > LimitWeights.second) ? LimitWeights.second : weight;

        resp.back().responses.push_back(weight);
      }

    } else {

      // Only want the CV response to be used in one of the dials, after the first
      // dial is found, all other dial responses should be /= CVResponse.
      double CVResponse = Get_MINERvA2p2h2EnergyDependencyScaling(
          e2i(simb_mode_copy::kMEC), true, Enu, ACV, BCV);

      CVResponse =
          (CVResponse < LimitWeights.first) ? LimitWeights.first : CVResponse;
      CVResponse =
          (CVResponse > LimitWeights.second) ? LimitWeights.second : CVResponse;

#ifdef MINERVAE2p2h_DEBUG
      std::cout << "[CV Response @ " << Enu << ", " << (ISLep->Pdg()) << ", "
                << ACV << ", " << BCV << "] = " << CVResponse << std::endl;
#endif

      bool UsedADial = false;
      if (pidx_A != kParamUnhandled<size_t>) {
        resp.push_back({md[pidx_A].systParamId, {}});
        for (double av : (*A_var)) {

          if(!nuMatched){
            resp.back().responses.push_back(1.);
            continue;
          }

          double weight = Get_MINERvA2p2h2EnergyDependencyScaling(
              e2i(simb_mode_copy::kMEC), true, Enu, av, BCV);
#ifdef MINERVAE2p2h_DEBUG
          std::cout << "[weight @ " << Enu << ", " << av << ", " << BCV
                    << "] = " << weight << std::endl;
#endif
          weight = (weight < LimitWeights.first) ? LimitWeights.first : weight;
          weight = (weight > LimitWeights.second) ? LimitWeights.second : weight;

          resp.back().responses.push_back(weight);
        }
        UsedADial = true;
      }
      if (pidx_B != kParamUnhandled<size_t>) {
        resp.push_back({md[pidx_B].systParamId, {}});
        for (double bv : (*B_var)) {

          if(!nuMatched){
            resp.back().responses.push_back(1.);
            continue;
          }

          double weight = Get_MINERvA2p2h2EnergyDependencyScaling(
              e2i(simb_mode_copy::kMEC), true, Enu, ACV, bv);
#ifdef MINERVAE2p2h_DEBUG
          std::cout << "[weight @ " << Enu << ", " << ACV << ", " << bv
                    << "] = " << weight << std::endl;
#endif

          weight = (weight < LimitWeights.first) ? LimitWeights.first : weight;
          weight = (weight > LimitWeights.second) ? LimitWeights.second : weight;

          if (UsedADial) {
            weight /= CVResponse;
          }
          resp.back().responses.push_back(weight);
        }
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

    weight = 1;

    for (auto const &eur : resp) {
      weight *= eur.responses[2];
    }

    valid_tree->Fill();
  }

  return resp;
}
std::string MINERvAE2p2h::AsString() { return "MINERvAE2p2h"; }

void MINERvAE2p2h::InitValidTree() {
  valid_file = new TFile("MINERvAE2p2h_valid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("weight", &weight);
}

MINERvAE2p2h::~MINERvAE2p2h() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
