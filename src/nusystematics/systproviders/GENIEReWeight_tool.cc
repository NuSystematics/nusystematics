#include "nusystematics/systproviders/GENIEReWeight_tool.hh"

#include "nusystematics/systproviders/GENIEReWeightEngineConfig.hh"
#include "nusystematics/systproviders/GENIEReWeightParamConfig.hh"

#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
#include "nusystematics/utility/GENIEUtils.hh"

// GENIE
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

#include "TH1.h"

#include <sstream>
#include <fstream>

using namespace fhicl;
using namespace systtools;
using namespace nusyst;

GENIEReWeight::GENIEReWeight(ParameterSet const &params)
    : IGENIESystProvider_tool(params), fHaveReconfiguredOneOfTheHERG(false),
      valid_file(nullptr), valid_tree(nullptr) {}

std::string GENIEReWeight::AsString() {
  CheckHaveMetaData();
  return "";
}

SystMetaData GENIEReWeight::BuildSystMetaData(ParameterSet const &params,
                                              paramId_t firstParamId) {

  tool_options = fhicl::ParameterSet();

  bool ignore_parameter_dependence =
      params.get<bool>("ignore_parameter_dependence", false);
  tool_options.put("ignore_parameter_dependence", ignore_parameter_dependence);

  bool UseFullHERG = params.get<bool>("UseFullHERG", false);
  tool_options.put("UseFullHERG", UseFullHERG);

  std::string genie_tune_name = params.get<std::string>("genie_tune_name",
                                                   "${GENIE_XSEC_TUNE}");
  tool_options.put("genie_tune_name",genie_tune_name);

  std::string evgen_list_name = params.get<std::string>("EventGeneratorList", "");
  tool_options.put("evgen_list_name",evgen_list_name);

  fill_valid_tree = params.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  SystMetaData QEmd =
      ConfigureQEParameterHeaders(params, firstParamId, tool_options);
  firstParamId += QEmd.size();

  SystMetaData MECmd =
      ConfigureMECParameterHeaders(params, firstParamId, tool_options);
  firstParamId += MECmd.size();

  SystMetaData NCELmd =
      ConfigureNCELParameterHeaders(params, firstParamId, tool_options);
  firstParamId += NCELmd.size();

  SystMetaData RESmd =
      ConfigureRESParameterHeaders(params, firstParamId, tool_options);
  firstParamId += RESmd.size();

  SystMetaData COHmd =
      ConfigureCOHParameterHeaders(params, firstParamId, tool_options);
  firstParamId += COHmd.size();

  SystMetaData DISmd =
      ConfigureDISParameterHeaders(params, firstParamId, tool_options);
  firstParamId += DISmd.size();

  SystMetaData FSImd =
      ConfigureFSIParameterHeaders(params, firstParamId, tool_options);
  firstParamId += FSImd.size();

  SystMetaData Othermd = ConfigureOtherParameterHeaders(params, firstParamId);
  firstParamId += Othermd.size();

  // Don't extend inline to make firstParamId incrementing more clear.
  ExtendSystMetaData(QEmd, DISmd);
  ExtendSystMetaData(QEmd, MECmd);
  ExtendSystMetaData(QEmd, NCELmd);
  ExtendSystMetaData(QEmd, RESmd);
  ExtendSystMetaData(QEmd, COHmd);
  ExtendSystMetaData(QEmd, FSImd);
  ExtendSystMetaData(QEmd, Othermd);

  return QEmd;
}

void GENIEReWeight::extend_ResponseToGENIEParameters(
    std::vector<GENIEResponseParameter> &&other) {
  for (auto &&o : other) {
    for (auto const &configured : ResponseToGENIEParameters) {
      if (configured.pidx == o.pidx) {
        std::cout << "[ERROR]: Attempted to merge GENIE GSyst response map, "
                     "but found duplicate response parameter index: "
                  << o.pidx << ", which corresponds to parameter: "
                  << std::quoted(GetSystMetaData()[o.pidx].prettyName)
                  << std::endl;
        throw;
      }
    }
    ResponseToGENIEParameters.push_back(std::move(o));
  }
}

bool GENIEReWeight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[INFO]: Setting up GENIE ReWeight instances..." << std::endl;

  //==== taken from https://github.com/LArSoft/larsim/blob/5b8ddeef9a556aa8a6bfa61915a16af25fa33a74/larsim/EventWeight/App/EventWeight_module.cc#L63-L78
  std::string genie_tune_name = tool_options.get<std::string>("genie_tune_name");
  std::string evgen_list_name = tool_options.get<std::string>("EventGeneratorList", "");
  std::cout << "[INFO]: genie_tune_name = " << genie_tune_name << std::endl;
  std::cout << "[INFO]: evgen_list_name = " << evgen_list_name << std::endl;

  // Constructor automatically calls grunopt->Init();
  genie::RunOpt* grunopt = genie::RunOpt::Instance();
  // SetEventGeneratorList wasn't introduced until R-3
  std::string expEvtGenListName = nusyst::ExpandEnvVar(evgen_list_name);
  if ( expEvtGenListName != "" ) {
    grunopt->SetEventGeneratorList(expEvtGenListName);
  }

  std::string expTuneName = nusyst::ExpandEnvVar(genie_tune_name);
  if ( expTuneName != genie_tune_name ) {
    std::cout << "TuneName started as '" << genie_tune_name << "' "
                             << " converted to " << expTuneName << std::endl;;
  }

  // If the XSecSplineList returns a non-empty string as the current tune name,
  // then genie::RunOpt::BuildTune() has already been called.
  std::string current_tune = genie::XSecSplineList::Instance()->CurrentTune();
  if ( current_tune.empty() ) {
    // We need to build the GENIE tune config
      std::cout << "Configuring GENIE tune \""
        << expTuneName << '\"' << std::endl;;
      grunopt->SetTuneName( expTuneName );
      grunopt->BuildTune();
      std::cout << *(grunopt->Tune()) << std::endl;
    }
    else {
      // It has already been built, so just check consistency
      if ( expTuneName != current_tune) {
        std::cout << "[TuneNameMismatch] Requested GENIE tune \""
          << expTuneName << "\" does not match previously built tune \""
          << current_tune << '\"' << std::endl;
      }
    }

  // Quiet mode
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_laconic.xml");

  extend_ResponseToGENIEParameters(
      ConfigureQEWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureMECWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureNCELWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureRESWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureCOHWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureDISWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureFSIWeightEngine(GetSystMetaData(), tool_options));

  extend_ResponseToGENIEParameters(
      ConfigureOtherWeightEngine(GetSystMetaData(), tool_options));

  std::cout << "[INFO]: Done!" << std::endl;

  fill_valid_tree = tool_options.get("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");
  return true;
}

// #define GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG

systtools::event_unit_response_t
GENIEReWeight::GetEventResponse(genie::EventRecord const &gev) {

  systtools::event_unit_response_t event_responses;
  size_t NResps = ResponseToGENIEParameters.size();

  for (size_t resp_idx = 0; resp_idx < NResps; ++resp_idx) {
    event_responses.push_back(GetEventGENIEParameterResponse(gev, resp_idx));
  }
  if (fill_valid_tree) {
    TLorentzVector FSLepP4 = gev.Summary()->Kine().FSLeptonP4();
    TLorentzVector ISLepP4 =
        *gev.Summary()->InitState().GetProbeP4(genie::kRfLab);
    TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

    Pdgnu = gev.Summary()->InitState().ProbePdg();
    NEUTMode = 0;
    if (gev.Summary()->ProcInfo().IsMEC() &&
        gev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&gev);
    }

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = gev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();
    valid_tree->Fill();
  }

  return event_responses;
}

double GENIEReWeight::GetEventWeightResponse(
    genie::EventRecord const &gev,
    systtools::param_value_list_t const &set_params) {

  double weight = 1;

  fHaveReconfiguredOneOfTheHERG = true;

  for (auto &GENIEResponse : ResponseToGENIEParameters) {
    for (auto const &dep : GENIEResponse.dependents) {
      SystParamHeader const &hdr = GetSystMetaData()[dep.pidx];
      double pval = hdr.centralParamValue;
      if (ContainterHasParam(set_params, dep.pidx)) {
        pval = GetParamElementFromContainer(set_params, dep.pidx).val;
      }

      GENIEResponse.Herg.front()->Systematics().Set(dep.gdial, pval);
    }
    GENIEResponse.Herg.front()->Reconfigure();
    bool is_set_dir = TH1::AddDirectoryStatus();
    if (!is_set_dir) {
      TH1::AddDirectory(true);
    }
    weight *= GENIEResponse.Herg.front()->CalcWeight(gev);
    if (!is_set_dir) {
      TH1::AddDirectory(false);
    }
  }

  return weight;
}

systtools::event_unit_response_t
GENIEReWeight::GetEventResponse(genie::EventRecord const &gev,
                                systtools::paramId_t pid) {

  size_t NResps = ResponseToGENIEParameters.size();
  for (size_t resp_idx = 0; resp_idx < NResps; ++resp_idx) {
    if (GetSystMetaData()[ResponseToGENIEParameters[resp_idx].pidx]
            .systParamId == pid) {
      return systtools::event_unit_response_t{
          {GetEventGENIEParameterResponse(gev, resp_idx)}};
    }
  }

  return systtools::event_unit_response_t();
}

namespace {
std::streambuf *default_cout;
std::ofstream redirect_stream("/dev/null");
} // namespace

systtools::ParamResponses
GENIEReWeight::GetEventGENIEParameterResponse(genie::EventRecord const &gev,
                                              size_t idx) {

  GENIEResponseParameter &GENIEResponse = ResponseToGENIEParameters[idx];
  systtools::SystParamHeader const &hdr = GetSystMetaData()[GENIEResponse.pidx];

  size_t NVars = hdr.isCorrection ? 1 : hdr.paramVariations.size();
#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
  std::cout << "[INFO]: Resp dial: " << hdr.prettyName << " with " << NVars
            << " variations of " << GENIEResponse.dependents.size()
            << " dependent parameters." << std::endl;
#endif

  // Have one GENIEReWeight per response rather than per variation, must
  // reconfigure.
  bool IsReducedHERG = (NVars > GENIEResponse.Herg.size());

  if (fHaveReconfiguredOneOfTheHERG && !IsReducedHERG) {
    throw invalid_engine_state()
        << "[ERROR]: GENIEReWeight_tool is configured to instantiate and "
           "Reconfigure one GENIE genie::GReWeight per variation. Because this "
           "instance has been used to get an externally specified variation by "
           "GetEventWeightResponse, these engines will no longer be correctly "
           "configured. It is advised to only use GetEventWeightResponse on a "
           "separate GENIEReWeight_tool instance than the one used for "
           "calculating the pre-configured event responses via "
           "GetEventResponse.";
  }

  ParamResponses presp{hdr.systParamId, {}};

  for (size_t var_it = 0; var_it < NVars; ++var_it) {

    if (IsReducedHERG) { // Need a reconfigure for each variation
      for (auto const &dep : GENIEResponse.dependents) {
        SystParamHeader const &hdr = GetSystMetaData()[dep.pidx];
        GENIEResponse.Herg.front()->Systematics().Set(
            dep.gdial, hdr.isCorrection ? hdr.centralParamValue
                                        : hdr.paramVariations[var_it]);
#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
        std::cout << "\t\t Var = "
                  << (hdr.isCorrection ? hdr.centralParamValue
                                       : hdr.paramVariations[var_it])
                  << " GDial: " << genie::rew::GSyst::AsString(dep.gdial)
                  << " at "
                  << GENIEResponse.Herg.front()
                         ->Systematics()
                         .Info(dep.gdial)
                         ->CurValue
                  << std::endl;

#endif
      }
      GENIEResponse.Herg.front()->Reconfigure();
      bool is_set_dir = TH1::AddDirectoryStatus();
      if (!is_set_dir) {
        TH1::AddDirectory(true);
      }
      presp.responses.push_back(GENIEResponse.Herg.front()->CalcWeight(gev));
      if (!is_set_dir) {
        TH1::AddDirectory(false);
      }
    } else { // Is full HERG, no reconfigure needed
#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
      for (GENIEResponseParameter::DependentParameter const &dep :
           GENIEResponse.dependents) {
        SystParamHeader const &hdr = GetSystMetaData()[dep.pidx];
        std::cout << "\t\t Var = "
                  << (hdr.isCorrection ? hdr.centralParamValue
                                       : hdr.paramVariations[var_it])
                  << " GDial: " << genie::rew::GSyst::AsString(dep.gdial)
                  << " at "
                  << GENIEResponse.Herg[var_it]
                         ->Systematics()
                         .Info(dep.gdial)
                         ->CurValue
                  << std::endl;
      }
#endif

      // As some GENIE dials are very slow, it is worth checking if we have
      // already calculated this variation before
      size_t pindx = std::numeric_limits<size_t>::max();
      for (size_t v_it = 0; v_it < var_it; ++v_it) {
        if (fabs(hdr.paramVariations[v_it] - hdr.paramVariations[var_it]) <
            1E-5) {
          pindx = v_it;
          break;
        }
      }
      if (pindx !=
          std::numeric_limits<size_t>::max()) { // Have already calculated this
        // value, just use that
        presp.responses.push_back(presp.responses[pindx]);
      } else { // must calculate
        bool is_set_dir = TH1::AddDirectoryStatus();
        if (!is_set_dir) {
          TH1::AddDirectory(true);
        }

        ::default_cout = std::cout.rdbuf();
        std::cout.rdbuf(::redirect_stream.rdbuf());
        presp.responses.push_back(GENIEResponse.Herg[var_it]->CalcWeight(gev));
        std::cout.rdbuf(::default_cout);

        if (!is_set_dir) {
          TH1::AddDirectory(false);
        }
      }
    }
#ifdef GENIEREWEIGHT_GETEVENTRESPONSE_DEBUG
    std::cout << "\t -> " << presp.responses.back() << std::endl;
#endif
  }

  return presp;
}

void GENIEReWeight::InitValidTree() {

  valid_file = new TFile("GENIEReWeightValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdgnu", &Pdgnu);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("weights", &weights);
}

GENIEReWeight::~GENIEReWeight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
