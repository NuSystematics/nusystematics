/*******************************************************************************
 * MECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "MECq0q3InterpWeighting_tool.hh"

#include <fhiclcpp/ParameterSet.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TString.h>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

// SystematicsTools helper
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

// GENIE
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"

using namespace nusyst;
using namespace systtools;

MECq0q3InterpWeighting::MECq0q3InterpWeighting(
    const fhicl::ParameterSet& p)
  : IGENIESystProvider_tool(p) {}

// ---------------------------------------------------------------------------
// Build metadata (standard NuSyst header parsing)
SystMetaData
MECq0q3InterpWeighting::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                                  systtools::paramId_t firstId) {

  std::cout << "[MECq0q3InterpWeighting::BuildSystMetaData] Called\n";

  SystMetaData smd;
  
  // Check if Q0 binning is enabled
  auto man = ps.get<fhicl::ParameterSet>("MECResponse_input_manifest");
  std::vector<double> q0Bins;
  if (man.has_key("Q0Bins")) {
    q0Bins = man.get<std::vector<double>>("Q0Bins");
  }
  
  // If Q0 binning is enabled, create multiple dials (one per bin)
  if (!q0Bins.empty() && q0Bins.size() >= 2) {
    std::cout << "  Q0-binned mode: creating " << (q0Bins.size() - 1) << " dials\n";
    
    for (size_t i = 0; i < q0Bins.size() - 1; ++i) {
      std::string dialName = "MECResponse_q0bin" + std::to_string(i);
      SystParamHeader phdr;
      if (ParseFhiclToolConfigurationParameter(ps, dialName, phdr, firstId)) {
        phdr.systParamId = firstId++;
        smd.push_back(phdr);
        std::cout << "    Created dial: " << dialName 
                  << " for q0 [" << q0Bins[i] << ", " << q0Bins[i+1] << ") GeV\n";
      }
    }
  } else {
    // Single dial mode (backward compatible)
    std::cout << "  Single-dial mode (backward compatible)\n";
    SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(ps, "MECResponse", phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  // stash manifest for SetupResponseCalculator
  tool_options.put("MECResponse_input_manifest", man);

  return smd;
}

// ---------------------------------------------------------------------------
// Read manifest and build calculators
bool
MECq0q3InterpWeighting::SetupResponseCalculator(fhicl::ParameterSet const &tool_opts)
{
  std::cout << "[MECq0q3InterpWeighting] SetupResponseCalculator begin\n";

  const auto manifest =
      tool_opts.get<fhicl::ParameterSet>("MECResponse_input_manifest");

  // Energy grid (required)
  if (!manifest.has_key("EnergyGrid"))
    throw std::runtime_error("Missing EnergyGrid");
  fEgrid = manifest.get<std::vector<double>>("EnergyGrid");
  if (fEgrid.empty())
    throw std::runtime_error("EnergyGrid must be non-empty");

  // Weight clamp (optional)
  if (manifest.has_key("WeightLimits")) {
    auto wl = manifest.get<std::vector<double>>("WeightLimits");
    if (wl.size() >= 2) {
      fWmin = std::min(wl.front(), wl.back());
      fWmax = std::max(wl.front(), wl.back());
    }
  }
  if (!(std::isfinite(fWmin) && std::isfinite(fWmax) && fWmin >= 0.0 && fWmax >= fWmin))
    throw std::runtime_error("Invalid WeightLimits");


  // NEW: q0 apply window (optional, defaults to [0, +inf))
  fQ0ApplyMin = manifest.get<double>("Q0ApplyMin", 0.0);
  if (manifest.has_key("Q0ApplyMax")) {
    fQ0ApplyMax = manifest.get<double>("Q0ApplyMax");
    if (!(std::isfinite(fQ0ApplyMax)))
      throw std::runtime_error("Q0ApplyMax must be finite if provided");
  } else {
    fQ0ApplyMax = std::numeric_limits<double>::infinity();
  }
  if (!(std::isfinite(fQ0ApplyMin) && fQ0ApplyMin >= 0.0))
    throw std::runtime_error("Q0ApplyMin must be finite and >= 0");
  if (!(fQ0ApplyMax > fQ0ApplyMin))
    throw std::runtime_error("Q0ApplyMax must be > Q0ApplyMin");

  // NEW: q3 apply window (optional, defaults to [0, +inf))
  fQ3ApplyMin = manifest.get<double>("Q3ApplyMin", 0.0);
  if (manifest.has_key("Q3ApplyMax")) {
    fQ3ApplyMax = manifest.get<double>("Q3ApplyMax");
    if (!(std::isfinite(fQ3ApplyMax)))
      throw std::runtime_error("Q3ApplyMax must be finite if provided");
  } else {
    fQ3ApplyMax = std::numeric_limits<double>::infinity();
  }
  if (!(std::isfinite(fQ3ApplyMin) && fQ3ApplyMin >= 0.0))
    throw std::runtime_error("Q3ApplyMin must be finite and >= 0");
  if (!(fQ3ApplyMax > fQ3ApplyMin))
    throw std::runtime_error("Q3ApplyMax must be > Q3ApplyMin");

  // NEW: Q0 binning for multiple dials (optional)
  if (manifest.has_key("Q0Bins")) {
    fQ0Bins = manifest.get<std::vector<double>>("Q0Bins");
    if (fQ0Bins.size() < 2)
      throw std::runtime_error("Q0Bins must have at least 2 edges to define bins");
    
    // Validate bins are in ascending order
    for (size_t i = 1; i < fQ0Bins.size(); ++i) {
      if (fQ0Bins[i] <= fQ0Bins[i-1])
        throw std::runtime_error("Q0Bins must be in strictly ascending order");
    }
    
    std::cout << "  Q0 binning enabled: " << (fQ0Bins.size() - 1) << " bins\n";
    std::cout << "  Bin edges:";
    for (double edge : fQ0Bins) std::cout << " " << edge;
    std::cout << " GeV\n";
    
    // Q0Select range incompatible with Q0Bins (they serve different purposes)
    if (manifest.has_key("Q0SelectMin") || manifest.has_key("Q0SelectMax")) {
      std::cout << "  WARNING: Q0SelectMin/Max ignored when Q0Bins is set\n";
    }
  } else {
    // Single dial mode - legacy Q0 selection range still works
    fQ0Bins.clear();
  }

  // NEW: Q0 selection range (optional, defaults to disabled)
  fQ0SelectMin = manifest.get<double>("Q0SelectMin", 0.0);
  fQ0SelectMax = manifest.get<double>("Q0SelectMax", 0.0);
  if (!(std::isfinite(fQ0SelectMin) && fQ0SelectMin >= 0.0))
    throw std::runtime_error("Q0SelectMin must be finite and >= 0");
  if (!(std::isfinite(fQ0SelectMax) && fQ0SelectMax >= 0.0))
    throw std::runtime_error("Q0SelectMax must be finite and >= 0");
  // If both are zero, selection is disabled; otherwise Max must be > Min
  if (fQ0SelectMax > 0.0 && fQ0SelectMax <= fQ0SelectMin)
    throw std::runtime_error("Q0SelectMax must be > Q0SelectMin when enabled");
  if (fQ0SelectMin > 0.0 || fQ0SelectMax > 0.0) {
    std::cout << "  Q0 selection range enabled: [" << fQ0SelectMin << ", " << fQ0SelectMax << "] GeV\n";
  }

  // Energy window & snap tolerance (defaults implement your request)
  fEnuMin     = manifest.get<double>("EnuMin",     0.4);
  fEnuMax     = manifest.get<double>("EnuMax",     2.5);
  fEnuSnapTol = manifest.get<double>("EnuSnapTol", 5e-3); // 5 MeV
  if (!(std::isfinite(fEnuMin) && std::isfinite(fEnuMax) && fEnuMax >= fEnuMin))
    throw std::runtime_error("EnuMin/EnuMax must be finite and EnuMax>=EnuMin");
  if (!(std::isfinite(fEnuSnapTol) && fEnuSnapTol >= 0.0))
    throw std::runtime_error("EnuSnapTol must be finite and >= 0");

  const bool mapIsQ3xQ0 = manifest.get<bool>("MapIsQ3xQ0", false);
  const bool useNearestBin = manifest.get<bool>("UseNearestBin", true);  // turn ON new behavior
  const bool edgeClamp     = manifest.get<bool>("EdgeClamp",     true);  // clamp OOR to edge bin
  
  // Read Model parameter early to determine default out-of-range behavior
  std::string model = manifest.get<std::string>("Model", "");
  
  // Out-of-range weight: default depends on model
  // Valencia: 0.0 (suppress beyond q3~1.2 GeV, q0~??? GeV)
  // Martini:  0.0 (suppress beyond q0~0.995 GeV)
  // Custom:   configurable via OutOfRangeWeight parameter
  double outOfRangeWeight = 1.0;  // backward compatible default
  if (manifest.has_key("OutOfRangeWeight")) {
    outOfRangeWeight = manifest.get<double>("OutOfRangeWeight");
  } else if (!model.empty()) {
    // Auto-set based on model to match native generator behavior
    if (model == "valencia" || model == "martini") {
      outOfRangeWeight = 0.0;  // Suppress events outside model's phase space
    }
  }
  
  std::cout << "  MapIsQ3xQ0     : " << (mapIsQ3xQ0 ? "true" : "false") << "\n";
  std::cout << "  UseNearestBin  : " << (useNearestBin ? "true" : "false") << "\n";
  std::cout << "  EdgeClamp      : " << (edgeClamp ? "true" : "false") << "\n";
  std::cout << "  OutOfRangeWeight: " << outOfRangeWeight << "\n";



  std::cout << "  EnergyGrid size: " << fEgrid.size() << "\n"
            << "  WeightLimits   : [" << fWmin << ", " << fWmax << "]\n"
            << "  Enu window     : [" << fEnuMin << ", " << fEnuMax << "] GeV\n"
            << "  Enu snap tol   : " << fEnuSnapTol << " GeV\n";





  // Histogram sourcing:
  //  (A) WeightFile with conventional names: h_weights_map_{np|nn}_<E>GeV
  //  (B) Arrays np_files/nn_files with explicit histogram names: HistNameNP/HistNameNN
  //  (C) Flavor-selected external files under NuisFlatWeightMapBaseDir
  //  (D) Legacy auto-generated file paths based on Model and DataBaseDir
  const bool haveSingle = manifest.has_key("WeightFile");

  std::string dataBaseDir =
      manifest.get<std::string>("DataBaseDir", "");
  dataBaseDir = systtools::expand_env_vars(dataBaseDir);

  const std::string nuisFlatBaseDir =
      manifest.get<std::string>("NuisFlatWeightMapBaseDir", "");

  std::vector<std::string> np_files, nn_files;

  auto set_legacy_model_info = [&](std::string& modelDir,
                                   std::string& filePrefix) {
    if (model == "valencia") {
      modelDir = "ValenciaMECq0q3";
      filePrefix = "reweight_data_SuSAv2_to_valencia";
    } else if (model == "martini") {
      modelDir = "martini_2p2h_weights";
      filePrefix = "reweight_data_SuSAv2_to_martini";
    } else {
      throw std::runtime_error("Unknown Model: '" + model + "'. Expected 'valencia' or 'martini'");
    }
  };

  auto external_model_tune = [&]() -> std::string {
    if (model == "valencia")
      return "g18_10a";
    if (model == "martini")
      return "g24_12a";
    throw std::runtime_error("Unknown Model: '" + model + "'. Expected 'valencia' or 'martini'");
  };

  auto flavor_tag = [](const std::string& flavor) -> std::string {
    if (flavor == "numu")
      return "NUMU";
    if (flavor == "nue")
      return "NUE";
    if (flavor == "numubar")
      return "NUMUBAR";
    if (flavor == "nuebar")
      return "NUEBAR";
    throw std::runtime_error("Unknown Flavor: '" + flavor + "'. Expected 'numu', 'nue', 'numubar', or 'nuebar'");
  };

  auto energy_label = [](double E) -> std::string {
    return std::string(Form("%0.2f", E));
  };

  auto is_numu_local_energy = [](double E) -> bool {
    constexpr double tol = 1e-6;
    return E >= (0.40 - tol) && E <= (2.50 + tol);
  };

  auto require_existing_map = [&](const std::string& fname,
                                  const std::string& flavor,
                                  const std::string& topo,
                                  double E) {
    std::ifstream fin(fname.c_str());
    if (!fin.good()) {
      throw std::runtime_error("Missing MEC q0/q3 weight map for Model '" +
                               model + "', Flavor '" + flavor +
                               "', topology '" + topo + "', energy " +
                               energy_label(E) + " GeV: " + fname);
    }
  };

  if (!haveSingle) {
    if (manifest.has_key("np_files") || manifest.has_key("nn_files")) {
      if (!(manifest.has_key("np_files") && manifest.has_key("nn_files"))) {
        throw std::runtime_error("Both np_files and nn_files are required when either is specified");
      }
      np_files = manifest.get<std::vector<std::string>>("np_files");
      nn_files = manifest.get<std::vector<std::string>>("nn_files");
    } else if (!nuisFlatBaseDir.empty()) {
      if (model.empty()) {
        throw std::runtime_error("NuisFlatWeightMapBaseDir requires Model to be 'valencia' or 'martini'");
      }

      const std::string flavor =
          manifest.get<std::string>("Flavor", "numu");
      const std::string flavorTag = flavor_tag(flavor);
      const std::string tune = external_model_tune();

      std::string legacyModelDir, legacyFilePrefix;
      set_legacy_model_info(legacyModelDir, legacyFilePrefix);

      std::cout << "[MECq0q3InterpWeighting] Auto-selecting model: " << model << "\n";
      std::cout << "  Flavor: " << flavor << "\n";
      std::cout << "  External weight map base: " << nuisFlatBaseDir << "\n";
      if (flavor == "numu") {
        std::cout << "  numu middle-grid folder: "
                  << nuisFlatBaseDir << "/" << model << "/" << flavor << "\n";
      }

      auto external_path = [&](const std::string& topo, double E) {
        return nuisFlatBaseDir + "/" + model + "/" + flavor +
               "/weight_map_ar23_to_" + tune + "_" + flavorTag + "_" +
               topo + "_" + energy_label(E) + "GeV.root";
      };

      auto numu_middle_grid_path = [&](const std::string& topo, double E) {
        return nuisFlatBaseDir + "/" + model + "/" + flavor + "/" + legacyFilePrefix +
               "_" + topo + "_" + energy_label(E) + "GeV.root";
      };

      for (double E : fEgrid) {
        const bool useNumuMiddleGrid = (flavor == "numu" && is_numu_local_energy(E));
        const std::string np_file = useNumuMiddleGrid ? numu_middle_grid_path("np", E)
                                                      : external_path("np", E);
        const std::string nn_file = useNumuMiddleGrid ? numu_middle_grid_path("nn", E)
                                                      : external_path("nn", E);
        require_existing_map(np_file, flavor, "np", E);
        require_existing_map(nn_file, flavor, "nn", E);
        np_files.push_back(np_file);
        nn_files.push_back(nn_file);
        std::cout << "  Generated np file: " << np_file << "\n";
        std::cout << "  Generated nn file: " << nn_file << "\n";
      }
    } else if (!model.empty() && !dataBaseDir.empty()) {
      std::string modelDir, filePrefix;
      set_legacy_model_info(modelDir, filePrefix);

      std::cout << "[MECq0q3InterpWeighting] Auto-selecting model: " << model << "\n";
      std::cout << "  Model directory: " << modelDir << "\n";

      for (double E : fEgrid) {
        const std::string np_file = dataBaseDir + "/" + modelDir + "/" +
            filePrefix + "_np_" + energy_label(E) + "GeV.root";
        const std::string nn_file = dataBaseDir + "/" + modelDir + "/" +
            filePrefix + "_nn_" + energy_label(E) + "GeV.root";
        np_files.push_back(np_file);
        nn_files.push_back(nn_file);
        std::cout << "  Generated np file: " << np_file << "\n";
        std::cout << "  Generated nn file: " << nn_file << "\n";
      }
    }
  }

  const bool haveArrays = !np_files.empty() && !nn_files.empty();
  if (!haveSingle && !haveArrays)
    throw std::runtime_error("Need either WeightFile, (np_files & nn_files), (NuisFlatWeightMapBaseDir & Model), or (Model & DataBaseDir)");

  fCalcs.clear();

  if (haveSingle) {
    const std::string fname = manifest.get<std::string>("WeightFile");
    TFile fin(fname.c_str(), "READ");
    if (!fin.IsOpen())
      throw std::runtime_error("Cannot open WeightFile: " + fname);

    for (Topo topo : {Topo::np, Topo::nn}) {
      const char* tag = (topo == Topo::np ? "np" : "nn");
      auto& vec = fCalcs[topo];
      vec.reserve(fEgrid.size());

      for (double E : fEgrid) {
  const std::string hname = Form("h_weights_map_%s_%0.2fGeV", tag, E);
        TH2D* h = dynamic_cast<TH2D*>(fin.Get(hname.c_str()));
        if (!h) throw std::runtime_error("Missing histogram '" + hname +
                                         "' in file " + fname);

        std::cout << "  Loaded " << fname << " :: " << hname
                  << "  X:[" << h->GetXaxis()->GetXmin() << "," << h->GetXaxis()->GetXmax() << "]"
                  << "  Y:[" << h->GetYaxis()->GetXmin() << "," << h->GetYaxis()->GetXmax() << "]\n";

        auto calc = std::make_unique<MECq0q3ResponseCalc>(h, fWmin, fWmax, mapIsQ3xQ0);
        calc->SetUseNearestBin(useNearestBin);
        calc->SetEdgeClamp(edgeClamp);
        calc->SetOutOfRangeWeight(outOfRangeWeight);
        vec.emplace_back(std::move(calc));
      }
    }
    fin.Close();
  } else {
    // arrays: explicit TH2 names are REQUIRED
    if (np_files.size() != fEgrid.size() || nn_files.size() != fEgrid.size())
      throw std::runtime_error("np_files/nn_files sizes must match EnergyGrid size");

    const std::string hname_np = manifest.get<std::string>("HistNameNP");
    const std::string hname_nn = manifest.get<std::string>("HistNameNN");

    auto load_list = [&](const std::vector<std::string>& files,
                         Topo topo, const std::string& hname) {
      auto& vec = fCalcs[topo];
      vec.reserve(files.size());
      for (size_t i = 0; i < files.size(); ++i) {
        const std::string& fname = files[i];
        TFile fin(fname.c_str(), "READ");
        if (!fin.IsOpen())
          throw std::runtime_error("Cannot open file: " + fname);
        TH2D* h = dynamic_cast<TH2D*>(fin.Get(hname.c_str()));
        if (!h)
          throw std::runtime_error("Missing histogram '" + hname + "' in file " + fname);

        std::cout << "  Loaded " << fname << " :: " << hname
                  << "  X:[" << h->GetXaxis()->GetXmin() << "," << h->GetXaxis()->GetXmax() << "]"
                  << "  Y:[" << h->GetYaxis()->GetXmin() << "," << h->GetYaxis()->GetXmax() << "]\n";

        auto calc = std::make_unique<MECq0q3ResponseCalc>(h, fWmin, fWmax, mapIsQ3xQ0);
        calc->SetUseNearestBin(useNearestBin);
        calc->SetEdgeClamp(edgeClamp);
        calc->SetOutOfRangeWeight(outOfRangeWeight);
        vec.emplace_back(std::move(calc));
        fin.Close();
      }
    };

    load_list(np_files, Topo::np, hname_np);
    load_list(nn_files, Topo::nn, hname_nn);
  }

  std::cout << "[MECq0q3InterpWeighting] SetupResponseCalculator done\n";
  return true;
}

// ---------------------------------------------------------------------------
// Event response
systtools::event_unit_response_t
MECq0q3InterpWeighting::GetEventResponse(genie::EventRecord const& ev)
{
  // classify topology
  const Topo topo = ClassifyEvent(ev);
  if (topo == Topo::unknown)
    return this->GetDefaultEventResponse();

  // compute leptonic transfers
  double q0 = 0.0, q3 = 0.0, Enu = 0.0;
  ComputeQ0Q3(ev, q0, q3, Enu);

  // Determine if Q0 selection is enabled.
  const bool q0SelectionEnabled = (fQ0SelectMax > 0.0);

  // q3/q0 apply-window gate. Outside the apply window the dial has no effect,
  // including in single-dial reweight mode.
  const bool outsideQ3Apply = (q3 <= fQ3ApplyMin + 1e-6 || q3 >= fQ3ApplyMax - 1e-6);
  const bool outsideQ0Apply = (q0 <= fQ0ApplyMin + 1e-6 || q0 >= fQ0ApplyMax - 1e-6);

  if (outsideQ3Apply || outsideQ0Apply)
    return this->GetDefaultEventResponse();

  // NEW: Q0 selection range: weight=1 if outside the selected range (if enabled)
  // This allows fine-grained control: e.g., apply reweight only in 0.05 < q0 < 0.1 GeV
  if (q0SelectionEnabled) {
    if (q0 < fQ0SelectMin - 1e-6 || q0 > fQ0SelectMax + 1e-6) {
      return this->GetDefaultEventResponse();  // weight=1 outside selection range
    }
  }

  // Energy guard: weight=1 outside [fEnuMin, fEnuMax] (after q3/q0 gates)
  if (Enu < fEnuMin - 1e-6 || Enu > fEnuMax + 1e-6)
    return this->GetDefaultEventResponse();

  // find bracketing energies (handles mono grid too)
  auto it_hi = std::lower_bound(fEgrid.begin(), fEgrid.end(), Enu);
  size_t ih = (it_hi == fEgrid.end()) ? fEgrid.size() - 1
                                      : std::distance(fEgrid.begin(), it_hi);
  size_t il = (ih == 0) ? 0 : ih - 1;

  const double Elo = fEgrid[il];
  const double Ehi = fEgrid[ih];

  // Exact-map snapping: if Enu is within EnuSnapTol of a grid point, use it
  double t = 0.0; // interpolation fraction
  if (std::fabs(Enu - Elo) <= fEnuSnapTol) {
    ih = il;          // exact low node
    t  = 0.0;
  } else if (std::fabs(Enu - Ehi) <= fEnuSnapTol) {
    il = ih;          // exact high node
    t  = 0.0;
  } else {
    t = (ih == il || Ehi <= Elo) ? 0.0
        : std::clamp((Enu - Elo) / (Ehi - Elo), 0.0, 1.0);
  }

  // central weights at the two energies (same index if snapped)
  const auto& vec = fCalcs.at(topo);
  const double w_lo = vec[il]->GetCentralWeight(q0, q3);
  const double w_hi = vec[ih]->GetCentralWeight(q0, q3);
  const double w_blend = (1.0 - t) * w_lo + t * w_hi;

  const double w_eff_cv = std::clamp(w_blend, fWmin, fWmax);
  const double one_sigma = (w_eff_cv - 1.0); // for variations

  // Build response vector
  systtools::event_unit_response_t response;
  
  // Q0-binned mode: multiple dials, only one active per event
  if (!fQ0Bins.empty()) {
    const int q0BinIdx = GetQ0BinIndex(q0);
    const auto& smd = this->GetSystMetaData();
    
    // Build response for all dials
    for (size_t dialIdx = 0; dialIdx < smd.size(); ++dialIdx) {
      const auto& hdr = smd[dialIdx];
      systtools::ParamResponses pr;
      pr.pid = hdr.systParamId;
      pr.responses.reserve(hdr.paramVariations.size());
      
      // If this event is in the current dial's q0 bin, use the calculated weight
      // Otherwise, return weight=1.0 (no effect)
      const bool isActiveDial = (static_cast<int>(dialIdx) == q0BinIdx);
      
      if (isActiveDial) {
        // Apply reweighting for this dial
        for (double d : hdr.paramVariations) {
          const double rw = std::clamp(1.0 + d * one_sigma, fWmin, fWmax);
          pr.responses.push_back(rw);
        }
        if (pr.responses.empty()) pr.responses.push_back(w_eff_cv);
      } else {
        // Inactive dial: return weight=1.0
        for (size_t i = 0; i < hdr.paramVariations.size(); ++i) {
          pr.responses.push_back(1.0);
        }
        if (pr.responses.empty()) pr.responses.push_back(1.0);
      }
      
      response.push_back(std::move(pr));
    }
  } 
  // Single dial mode (backward compatible)
  else {
    if (!this->GetSystMetaData().empty()) {
      const auto &hdr = this->GetSystMetaData()[0];
      systtools::ParamResponses pr;
      pr.pid = hdr.systParamId;
      pr.responses.reserve(hdr.paramVariations.size());
      for (double d : hdr.paramVariations) {
        const double rw = std::clamp(1.0 + d * one_sigma, fWmin, fWmax);
        pr.responses.push_back(rw);
      }
      if (pr.responses.empty()) pr.responses.push_back(w_eff_cv);
      response.push_back(std::move(pr));
    } else {
      systtools::ParamResponses pr;
      pr.pid = 0;
      pr.responses = { w_eff_cv };
      response.push_back(std::move(pr));
    }
  }

  return response;
}

/*******************************************************************************
 *  Helpers
 ******************************************************************************/

// q0, q3, Enu from (probe – final-state lepton)
void
MECq0q3InterpWeighting::ComputeQ0Q3(genie::EventRecord const& ev,
                                            double& q0, double& q3, double& Enu)
{
  const TLorentzVector p4nu  = *ev.Probe()->P4();
  const TLorentzVector p4lep = *ev.FinalStatePrimaryLepton()->P4();
  Enu = p4nu.E();

  TLorentzVector qv = p4nu - p4lep;  // four-momentum transfer
  q0 = qv.E();
  q3 = qv.P();

  if (q0 < 0.0 && std::abs(q0) < 1e-6) q0 = 0.0; // numerical safety
}

// classify 2N initial cluster
MECq0q3InterpWeighting::Topo
MECq0q3InterpWeighting::ClassifyEvent(genie::EventRecord const& ev)
{
  for (int i = 0; i < ev.GetEntries(); ++i) {
    const auto* p = ev.Particle(i);
    if (!p) continue;
    if (p->Status() != genie::kIStNucleonTarget) continue;
    const int pdg = p->Pdg();
    if (pdg == genie::kPdgClusterNN) return Topo::nn; // 2n
    if (pdg == genie::kPdgClusterNP) return Topo::np; // 1n+1p
  }
  return Topo::unknown;
}

// Determine which q0 bin (dial index) an event belongs to
// Returns -1 if outside all bins
int
MECq0q3InterpWeighting::GetQ0BinIndex(double q0) const
{
  if (fQ0Bins.empty()) return 0;  // Single dial mode
  
  // Find the bin: [edge_i, edge_{i+1})
  for (size_t i = 0; i < fQ0Bins.size() - 1; ++i) {
    if (q0 >= fQ0Bins[i] && q0 < fQ0Bins[i+1]) {
      return static_cast<int>(i);
    }
  }
  
  // Special case: include upper edge in the last bin
  if (q0 == fQ0Bins.back()) {
    return static_cast<int>(fQ0Bins.size() - 2);
  }
  
  return -1;  // Outside all bins
}
