#include "nusystematics/systproviders/CCQEXSecCorr_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"

#include "TLorentzVector.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

NEW_SYSTTOOLS_EXCEPT(invalid_QECorr_DataBaseDir_FILEPATH);

// ---------------------------------------------------------------------------
// Flavor string used in the ROOT file naming convention
static std::string FlavorTag(int pdg) {
  switch (pdg) {
    case  12: return "nue";
    case -12: return "nuebar";
    case  14: return "numu";
    case -14: return "numubar";
    default:  return "";
  }
}

// ---------------------------------------------------------------------------
CCQEXSecCorr::CCQEXSecCorr(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

// ---------------------------------------------------------------------------
SystMetaData CCQEXSecCorr::BuildSystMetaData(ParameterSet const &cfg,
                                              paramId_t firstId) {

  std::cout << "[CCQEXSecCorr::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "CCQEXSecCorr", phdr,
                                           firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet manifest =
      cfg.get<fhicl::ParameterSet>("CCQEXSecCorr_input_manifest");
  tool_options.put("CCQEXSecCorr_input_manifest", manifest);

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

// ---------------------------------------------------------------------------
bool CCQEXSecCorr::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[CCQEXSecCorr::SetupResponseCalculator] called" << std::endl;

  // --- Read manifest ---
  fhicl::ParameterSet manifest =
      tool_options.get<fhicl::ParameterSet>("CCQEXSecCorr_input_manifest");

  std::string dataBaseDir = manifest.get<std::string>("DataBaseDir");
  std::string filePattern =
      manifest.get<std::string>("file_pattern",
          "reweight_data_v3_04_02a_to_v3_06_02_sbn2_CCQE_%FLAVOR%_%ENERGY%_flat.root");
  std::string histName =
      manifest.get<std::string>("histogram_name", "h_weights_map");

  dataBaseDir = systtools::expand_env_vars(dataBaseDir);

  // Energy grid
  fEgrid = manifest.get<std::vector<double>>("EnergyGrid");
  if (fEgrid.empty())
    throw std::runtime_error("[CCQEXSecCorr] EnergyGrid must be non-empty");

  // Weight clamp
  if (manifest.has_key("WeightLimits")) {
    auto wl = manifest.get<std::vector<double>>("WeightLimits");
    if (wl.size() >= 2) {
      fWmin = std::min(wl.front(), wl.back());
      fWmax = std::max(wl.front(), wl.back());
    }
  }

  // Energy window & snap tolerance
  fEnuMin     = manifest.get<double>("EnuMin", 0.20);
  fEnuMax     = manifest.get<double>("EnuMax", 3.00);
  fEnuSnapTol = manifest.get<double>("EnuSnapTol", 5e-3);

  // Histogram interpretation options
  const bool useNearestBin = manifest.get<bool>("UseNearestBin", true);
  const bool edgeClamp     = manifest.get<bool>("EdgeClamp", true);

  std::cout << "  DataBaseDir    : " << dataBaseDir << "\n"
            << "  EnergyGrid size: " << fEgrid.size() << "\n"
            << "  WeightLimits   : [" << fWmin << ", " << fWmax << "]\n"
            << "  Enu window     : [" << fEnuMin << ", " << fEnuMax << "] GeV\n"
            << "  Enu snap tol   : " << fEnuSnapTol << " GeV\n"
            << "  UseNearestBin  : " << (useNearestBin ? "true" : "false") << "\n"
            << "  EdgeClamp      : " << (edgeClamp ? "true" : "false") << "\n";

  // --- Load histograms for each flavor x energy ---
  const std::vector<int> pdg_codes = {12, -12, 14, -14};

  fCalcs.clear();
  for (int pdg : pdg_codes) {
    std::string flav = FlavorTag(pdg);
    auto &vec = fCalcs[pdg];
    vec.reserve(fEgrid.size());

    for (double E : fEgrid) {
      // Build filename: substitute %FLAVOR% and %ENERGY%
      std::string fname = filePattern;

      // Replace %FLAVOR%
      auto pos_f = fname.find("%FLAVOR%");
      if (pos_f != std::string::npos)
        fname.replace(pos_f, 8, flav);

      // Replace %ENERGY% with energy formatted as "X.XX"
      auto pos_e = fname.find("%ENERGY%");
      if (pos_e != std::string::npos)
        fname.replace(pos_e, 8, Form("%0.2f", E));

      std::string fullpath = dataBaseDir + "/" + fname;

      TFile fin(fullpath.c_str(), "READ");
      if (!fin.IsOpen())
        throw std::runtime_error(
            "[CCQEXSecCorr] Cannot open file: " + fullpath);

      TH2D *h = dynamic_cast<TH2D *>(fin.Get(histName.c_str()));
      if (!h)
        throw std::runtime_error(
            "[CCQEXSecCorr] Missing histogram '" + histName +
            "' in file " + fullpath);

      std::cout << "  Loaded " << fullpath
                << "  X:[" << h->GetXaxis()->GetXmin() << ","
                << h->GetXaxis()->GetXmax() << "]"
                << "  Y:[" << h->GetYaxis()->GetXmin() << ","
                << h->GetYaxis()->GetXmax() << "]\n";

      auto calc =
          std::make_unique<MECq0q3ResponseCalc>(h, fWmin, fWmax);
      calc->SetUseNearestBin(useNearestBin);
      calc->SetEdgeClamp(edgeClamp);
      calc->SetOutOfRangeWeight(1.0); // weight=1 outside histogram domain
      vec.emplace_back(std::move(calc));

      fin.Close();
    }
  }

  // --- Parameter index ---
  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "CCQEXSecCorr");

  // --- Validation tree ---
  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  std::cout << "[CCQEXSecCorr::SetupResponseCalculator] done" << std::endl;
  return true;
}

// ---------------------------------------------------------------------------
event_unit_response_t
CCQEXSecCorr::GetEventResponse(genie::EventRecord const &ev) {

  // Only apply to CC-QE events
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

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = ISLepP4 - FSLepP4;

  double ev_Enu = ISLepP4.E();
  double ev_q0  = emTransfer.E();
  double ev_q3  = emTransfer.Vect().Mag();
  int    nu_pdg = ISLep->Pdg();

  // Check that we have calculators for this flavor
  auto it = fCalcs.find(nu_pdg);
  if (it == fCalcs.end()) {
    return this->GetDefaultEventResponse();
  }

  // Energy guard: weight=1 outside configured range
  if (ev_Enu < fEnuMin - 1e-6 || ev_Enu > fEnuMax + 1e-6) {
    return this->GetDefaultEventResponse();
  }

  // --- Energy interpolation ---
  const auto &calcs = it->second;

  auto it_hi = std::lower_bound(fEgrid.begin(), fEgrid.end(), ev_Enu);
  size_t ih = (it_hi == fEgrid.end()) ? fEgrid.size() - 1
                                      : std::distance(fEgrid.begin(), it_hi);
  size_t il = (ih == 0) ? 0 : ih - 1;

  const double Elo = fEgrid[il];
  const double Ehi = fEgrid[ih];

  // Snap to exact grid point if within tolerance
  double t = 0.0;
  if (std::fabs(ev_Enu - Elo) <= fEnuSnapTol) {
    ih = il;
    t  = 0.0;
  } else if (std::fabs(ev_Enu - Ehi) <= fEnuSnapTol) {
    il = ih;
    t  = 0.0;
  } else {
    t = (ih == il || Ehi <= Elo)
            ? 0.0
            : std::clamp((ev_Enu - Elo) / (Ehi - Elo), 0.0, 1.0);
  }

  const double w_lo = calcs[il]->GetCentralWeight(ev_q0, ev_q3);
  const double w_hi = calcs[ih]->GetCentralWeight(ev_q0, ev_q3);
  double w_blend = (1.0 - t) * w_lo + t * w_hi;

  // Clamp
  w_blend = std::clamp(w_blend, fWmin, fWmax);

  // --- Build response ---
  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  event_unit_response_t resp;
  resp.push_back({hdr.systParamId, {}});
  for (double var : hdr.paramVariations) {
    (void)var; // correction dial — one weight regardless of variation value
    resp.back().responses.push_back(w_blend);
  }
  // If no paramVariations (pure correction), still return the weight
  if (resp.back().responses.empty())
    resp.back().responses.push_back(w_blend);

  // --- Validation tree ---
  if (fill_valid_tree) {
    pdgfslep    = FSLep->Pdg();
    momfslep    = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();
    Pdgnu       = nu_pdg;
    Enu         = ev_Enu;
    Q2          = -emTransfer.Mag2();
    W           = ev.Summary()->Kine().W(true);
    q0          = ev_q0;
    q3          = ev_q3;
    weight      = w_blend;
    NEUTMode    = genie::utils::ghep::NeutReactionCode(&ev);
    valid_tree->Fill();
  }

  return resp;
}

// ---------------------------------------------------------------------------
std::string CCQEXSecCorr::AsString() { return "CCQEXSecCorr"; }

// ---------------------------------------------------------------------------
void CCQEXSecCorr::InitValidTree() {
  valid_file = new TFile("CCQEXSecCorrValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
  valid_tree->Branch("weight", &weight);
}

// ---------------------------------------------------------------------------
CCQEXSecCorr::~CCQEXSecCorr() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
