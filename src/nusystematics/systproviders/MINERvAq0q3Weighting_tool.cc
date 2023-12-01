#include "nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh"

#include "nusystematics/responsecalculators/MINERvA2p2hq0q3.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

MINERvAq0q3Weighting::MINERvAq0q3Weighting(ParameterSet const &params)
    : IGENIESystProvider_tool(params), RPATemplateReweighter(nullptr),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData MINERvAq0q3Weighting::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  SystMetaData smd;
  if (cfg.get<bool>("use_MINERvA_RPA_tunes", false)) { // RPA
    systtools::SystParamHeader param;
    param.systParamId = firstId++;

    param.prettyName = "MINERvATune_RPA";
    param.centralParamValue = 0;
    param.paramVariations = std::vector<double>{-1, 0, 1};

    if (!cfg.has_key("MINERvATune_RPA_input_manifest") ||
        !cfg.is_key_to_table("MINERvATune_RPA_input_manifest")) {
      throw invalid_ToolConfigurationFHiCL()
          << "[ERROR]: When configuring calculated variations for "
             "MINERvATune_RPA, expected to find a FHiCL table keyed by "
             "MINERvATune_RPA_input_manifest describing the location of the "
             "histogram inputs. See "
             "nusystematics/responsecalculators/"
             "TemplateResponseCalculatorBase.hh "
             "for the layout.";
    }

    fhicl::ParameterSet ps =
        cfg.get<fhicl::ParameterSet>("MINERvATune_RPA_input_manifest");
    tool_options.put("MINERvATune_RPA_input_manifest", ps);

    smd.push_back(param);
  }

  if (cfg.get<bool>("use_MINERvA_2p2h_tunes", false)) { // 2p2h

    parameter_per_2p2h_universe =
        cfg.get<bool>("parameter_per_2p2h_universe", false);

    tool_options.put("parameter_per_2p2h_universe",
                     parameter_per_2p2h_universe);

    if (parameter_per_2p2h_universe) {
      systtools::SystParamHeader param_CV, param_NN, param_np, param_QE;
      if (ParseFhiclToolConfigurationParameter(
              cfg, "Mnv2p2hGaussEnhancement_CV", param_CV, firstId)) {
        param_CV.systParamId = firstId++;
        smd.push_back(param_CV);
      }
      if (ParseFhiclToolConfigurationParameter(
              cfg, "Mnv2p2hGaussEnhancement_NN", param_NN, firstId)) {
        param_NN.systParamId = firstId++;
        smd.push_back(param_NN);
      }
      if (ParseFhiclToolConfigurationParameter(
              cfg, "Mnv2p2hGaussEnhancement_np", param_np, firstId)) {
        param_np.systParamId = firstId++;
        smd.push_back(param_np);
      }
      if (ParseFhiclToolConfigurationParameter(
              cfg, "Mnv2p2hGaussEnhancement_QE", param_QE, firstId)) {
        param_QE.systParamId = firstId++;
        smd.push_back(param_QE);
      }

    } else {
      systtools::SystParamHeader param;
      if (ParseFhiclToolConfigurationParameter(
              cfg, "Mnv2p2hGaussEnhancement", param, firstId)) {
        param.systParamId = firstId++;
      } else {
        param.systParamId = firstId++;
        param.centralParamValue = 1;
        param.paramVariations = {1, 2, 3, 4};
        param.prettyName = "Mnv2p2hGaussEnhancement";
      }
      smd.push_back(param);
    }
  }

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  MEC_LimitWeights = cfg.get<std::pair<double, double>>(
      "Mnv2p2hGaussEnhancement_LimitWeights",
      {0, std::numeric_limits<double>::max()});
  tool_options.put(
      "Mnv2p2hGaussEnhancement_LimitWeights",
      std::vector<double>{MEC_LimitWeights.first, MEC_LimitWeights.second});

  return smd;
}

bool MINERvAq0q3Weighting::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  if (HasParam(GetSystMetaData(), "MINERvATune_RPA")) {
    ConfiguredParameters[param_t::kMINERvARPA] =
        GetParamIndex(GetSystMetaData(), "MINERvATune_RPA");

    if (!tool_options.has_key("MINERvATune_RPA_input_manifest")) {
      throw systtools::invalid_ToolOptions()
          << "[ERROR]: MINERvATune_RPA parameter exists in the SystMetaData, "
             "but no MINERvATune_RPA_input_manifest key can be found on the "
             "tool_options table. This reweighting requires input histograms "
             "that must be specified. This should have been caught by  "
             "MINERvAq0q3Weighting::BuildSystMetaData, but wasn't, this is a "
             "bug, please report to the maintiner.";
    }

    RPATemplateReweighter = std::make_unique<MINERvARPAq0q3_ReWeight>(
        tool_options.get<fhicl::ParameterSet>(
            "MINERvATune_RPA_input_manifest"));
  }

  if (HasParam(GetSystMetaData(), "Mnv2p2hGaussEnhancement")) {
    ConfiguredParameters[param_t::kMINERvA2p2h] =
        GetParamIndex(GetSystMetaData(), "Mnv2p2hGaussEnhancement");

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h]];

    if (hdr.isCorrection) {
      vals_2p2hTotal.push_back(hdr.centralParamValue);
    } else {
      vals_2p2hTotal = hdr.paramVariations;
    }
  }
  if (HasParam(GetSystMetaData(), "Mnv2p2hGaussEnhancement_CV")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_CV] =
        GetParamIndex(GetSystMetaData(), "Mnv2p2hGaussEnhancement_CV");

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_CV]];

    if (hdr.isCorrection) {
      vals_2p2hCV.push_back(hdr.centralParamValue);
    } else {
      vals_2p2hCV = hdr.paramVariations;
    }
  }
  if (HasParam(GetSystMetaData(), "Mnv2p2hGaussEnhancement_NN")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_NN] =
        GetParamIndex(GetSystMetaData(), "Mnv2p2hGaussEnhancement_NN");

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_NN]];

    if (hdr.isCorrection) {
      vals_2p2hNN.push_back(hdr.centralParamValue);
    } else {
      vals_2p2hNN = hdr.paramVariations;
    }
  }
  if (HasParam(GetSystMetaData(), "Mnv2p2hGaussEnhancement_np")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_np] =
        GetParamIndex(GetSystMetaData(), "Mnv2p2hGaussEnhancement_np");

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_np]];

    if (hdr.isCorrection) {
      vals_2p2hnp.push_back(hdr.centralParamValue);
    } else {
      vals_2p2hnp = hdr.paramVariations;
    }
  }
  if (HasParam(GetSystMetaData(), "Mnv2p2hGaussEnhancement_QE")) {
    ConfiguredParameters[param_t::kMINERvA2p2h_QE] =
        GetParamIndex(GetSystMetaData(), "Mnv2p2hGaussEnhancement_QE");

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_QE]];

    if (hdr.isCorrection) {
      vals_2p2hQE.push_back(hdr.centralParamValue);
    } else {
      vals_2p2hQE = hdr.paramVariations;
    }
  }

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  MEC_LimitWeights = tool_options.get<std::pair<double, double>>(
      "Mnv2p2hGaussEnhancement_LimitWeights",
      {0, std::numeric_limits<double>::max()});

  return true;
}

double MINERvAq0q3Weighting::GetMINERvARPATuneWeight(double val, double q0,
                                                     double q3) {
  MINERvARPAq0q3_ReWeight::RPATweak_t tval;
  if (val == 0) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kCV;
  } else if (val == 1) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kPlus1;
  } else if (val == -1) {
    tval = MINERvARPAq0q3_ReWeight::RPATweak_t::kMinus1;
  } else {
    throw invalid_parameter_value()
        << "[ERROR]: When applying MINERvA RPA tune expected to find "
           "parameter values of [ 0 == CV, 1 == Plus1, -1 == Minus1 ], "
           "but found "
        << val;
  }

  return RPATemplateReweighter->GetWeight(q0, q3, tval);
}

double MINERvAq0q3Weighting::GetMINERvA2p2hTuneEnhancement(
    int val, double q0, double q3, QELikeTarget_t QELTarget) {
  std::array<double, 6> const *GaussParams;
  if (val == 1) {
    if (!((QELTarget == QELikeTarget_t::kNN) ||
          (QELTarget == QELikeTarget_t::knp))) {
      return 1;
    }

    GaussParams = &Gauss2DParams_CV;
  } else if (val == 2) {

    if (QELTarget != QELikeTarget_t::kNN) {
      return 1;
    }

    GaussParams = &Gauss2DParams_NNOnly;
  } else if (val == 3) {

    if (QELTarget != QELikeTarget_t::knp) {
      return 1;
    }

    GaussParams = &Gauss2DParams_npOnly;
  } else if (val == 4) {

    if (QELTarget != QELikeTarget_t::kQE) {
      return 1;
    }

    GaussParams = &Gauss2DParams_1p1hOnly;
  } else {
    throw invalid_parameter_value()
        << "[ERROR]: When applying MINERvA 2p2h tune expected to find "
           "parameter values of [ 1 == CV, 2 == NN, 3 == np, 4 == 1p1h ], "
           "but found "
        << val;
  }
  return 1 + Gaussian2D(q0, q3, *GaussParams);
}

event_unit_response_t
MINERvAq0q3Weighting::GetEventResponse(genie::EventRecord const &ev) {

  // make default response for configured parameter
  event_unit_response_t resp;

  if (!ev.Summary()->ProcInfo().IsWeakCC()) {
    return this->GetDefaultEventResponse();
  }

  if (!(ev.Summary()->ProcInfo().IsQuasiElastic() ||
        ev.Summary()->ProcInfo().IsMEC()) ||
      ev.Summary()->ExclTag().IsCharmEvent()) {
    return this->GetDefaultEventResponse();
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
  std::array<double, 2> q0q3{{emTransfer.E(), emTransfer.Vect().Mag()}};

  if (ConfiguredParameters.find(param_t::kMINERvARPA) !=
      ConfiguredParameters.end()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvARPA]];

    resp.push_back({hdr.systParamId, {}});
    if (hdr.isCorrection) {
      resp.back().responses.push_back(
          GetMINERvARPATuneWeight(hdr.centralParamValue, q0q3[0], q0q3[1]));
    } else {
      for (double var : hdr.paramVariations) {
        resp.back().responses.push_back(
            GetMINERvARPATuneWeight(var, q0q3[0], q0q3[1]));
      }
    }
  }

  QELikeTarget_t qel_targ = GetQELikeTarget(ev);

  // Only ever applies to 2p2h/qe events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h) !=
       ConfiguredParameters.end()) &&
      (qel_targ != QELikeTarget_t::kInvalidTopology)) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h]];

    resp.push_back({hdr.systParamId, {}});
    for (double var : vals_2p2hTotal) {
      double wght =
          GetMINERvA2p2hTuneEnhancement(var, q0q3[0], q0q3[1], qel_targ);
      wght = (wght < MEC_LimitWeights.first) ? MEC_LimitWeights.first : wght;
      wght = (wght > MEC_LimitWeights.second) ? MEC_LimitWeights.second : wght;
      resp.back().responses.push_back(wght);
    }
  }

  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_CV) !=
       ConfiguredParameters.end()) &&
      ev.Summary()->ProcInfo().IsMEC()) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_CV]];

    resp.push_back({hdr.systParamId, {}});
    for (double v : vals_2p2hCV) {
      double cv_weight =
          1 + v * GetMINERvA2p2hTuneEnhancement(1, q0q3[0], q0q3[1], qel_targ);

      cv_weight = (cv_weight < MEC_LimitWeights.first) ? MEC_LimitWeights.first
                                                       : cv_weight;
      cv_weight = (cv_weight > MEC_LimitWeights.second)
                      ? MEC_LimitWeights.second
                      : cv_weight;

      resp.back().responses.push_back(cv_weight);
    }
  }
  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_NN) !=
       ConfiguredParameters.end()) &&
      qel_targ == QELikeTarget_t::kNN) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_NN]];

    resp.push_back({hdr.systParamId, {}});
    for (double v : vals_2p2hNN) {
      double tune_ench =
          v * GetMINERvA2p2hTuneEnhancement(2, q0q3[0], q0q3[1], qel_targ);

      tune_ench = (tune_ench < MEC_LimitWeights.first) ? MEC_LimitWeights.first
                                                       : tune_ench;
      tune_ench = (tune_ench > MEC_LimitWeights.second)
                      ? MEC_LimitWeights.second
                      : tune_ench;

      resp.back().responses.push_back(tune_ench);
    }
  }
  // Only ever applies to 2p2h events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_np) !=
       ConfiguredParameters.end()) &&
      qel_targ == QELikeTarget_t::knp) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_np]];

    resp.push_back({hdr.systParamId, {}});
    for (double v : vals_2p2hnp) {
      double tune_ench =
          v * GetMINERvA2p2hTuneEnhancement(3, q0q3[0], q0q3[1], qel_targ);

      tune_ench = (tune_ench < MEC_LimitWeights.first) ? MEC_LimitWeights.first
                                                       : tune_ench;
      tune_ench = (tune_ench > MEC_LimitWeights.second)
                      ? MEC_LimitWeights.second
                      : tune_ench;

      resp.back().responses.push_back(tune_ench);
    }
  }
  // Only ever applies to qe events
  if ((ConfiguredParameters.find(param_t::kMINERvA2p2h_QE) !=
       ConfiguredParameters.end()) &&
      qel_targ == QELikeTarget_t::kQE) {

    SystParamHeader const &hdr =
        GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h_QE]];

    resp.push_back({hdr.systParamId, {}});
    for (double v : vals_2p2hQE) {
      double tune_ench =
          v * GetMINERvA2p2hTuneEnhancement(4, q0q3[0], q0q3[1], qel_targ);

      tune_ench = (tune_ench < MEC_LimitWeights.first) ? MEC_LimitWeights.first
                                                       : tune_ench;
      tune_ench = (tune_ench > MEC_LimitWeights.second)
                      ? MEC_LimitWeights.second
                      : tune_ench;

      resp.back().responses.push_back(tune_ench);
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

    QELTarget = e2i(qel_targ);

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    RPA_weights.clear();
    MEC_weights.clear();

    if (ConfiguredParameters.find(param_t::kMINERvARPA) !=
        ConfiguredParameters.end()) {
      paramId_t RPA_param =
          GetSystMetaData()[ConfiguredParameters[param_t::kMINERvARPA]]
              .systParamId;
      RPA_weights = GetParamElementFromContainer(resp, RPA_param).responses;
    }
    if ((ConfiguredParameters.find(param_t::kMINERvA2p2h) !=
         ConfiguredParameters.end()) &&
        (ev.Summary()->ProcInfo().IsQuasiElastic() ||
         ev.Summary()->ProcInfo().IsMEC())) {
      paramId_t MEC_param =
          GetSystMetaData()[ConfiguredParameters[param_t::kMINERvA2p2h]]
              .systParamId;
      MEC_weights = GetParamElementFromContainer(resp, MEC_param).responses;
    }
    for (param_t tune_2p2h_universe :
         {param_t::kMINERvA2p2h_CV, param_t::kMINERvA2p2h_NN,
          param_t::kMINERvA2p2h_np, param_t::kMINERvA2p2h_QE}) {
      if ((ConfiguredParameters.find(tune_2p2h_universe) !=
           ConfiguredParameters.end()) &&
          (ev.Summary()->ProcInfo().IsQuasiElastic() ||
           ev.Summary()->ProcInfo().IsMEC())) {
        paramId_t MEC_param =
            GetSystMetaData()[ConfiguredParameters[tune_2p2h_universe]]
                .systParamId;
        size_t idx = GetParamContainerIndex(resp, MEC_param);
        if (idx != kParamUnhandled<size_t>) {
          MEC_weights.push_back(resp[idx].responses.front());
        }
      }
    }

    nRPA_weights = RPA_weights.size();
    nMEC_weights = MEC_weights.size();
    valid_tree->Fill();
  }

  return resp;
}

std::string MINERvAq0q3Weighting::AsString() { return ""; }

void MINERvAq0q3Weighting::InitValidTree() {

  valid_file = new TFile("MINERvAq0q3WeightValid.root", "RECREATE");
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
  valid_tree->Branch("nRPA_weights", &nRPA_weights);
  valid_tree->Branch("nMEC_weights", &nMEC_weights);
  valid_tree->Branch("RPA_weights", &RPA_weights);
  valid_tree->Branch("MEC_weights", &MEC_weights);
}

MINERvAq0q3Weighting::~MINERvAq0q3Weighting() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
