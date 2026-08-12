/*******************************************************************************
 * MECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "MECq0q3InterpWeighting_tool.hh"

#include <fhiclcpp/ParameterSet.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <utility>

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

NEW_SYSTTOOLS_EXCEPT(invalid_MEC_DataBaseDir_FILEPATH);

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

  auto flavor_pdg = [](const std::string& flavor) -> int {
    if (flavor == "numu") return genie::kPdgNuMu;
    if (flavor == "nue") return genie::kPdgNuE;
    if (flavor == "numubar") return genie::kPdgAntiNuMu;
    if (flavor == "nuebar") return genie::kPdgAntiNuE;
    throw std::runtime_error("Unknown flavor: '" + flavor +
                             "'. Expected 'numu', 'nue', 'numubar', or 'nuebar'");
  };

  if (manifest.has_key("Flavor") && manifest.has_key("Flavors"))
    throw std::runtime_error("Specify either Flavor or Flavors, not both");

  // A zero PDG key preserves the legacy behavior: one map set applies to all
  // incoming flavors. Explicit Flavor/Flavors entries are keyed by probe PDG.
  std::vector<std::pair<std::string, int>> requestedFlavors;
  if (manifest.has_key("Flavors")) {
    const auto flavors = manifest.get<std::vector<std::string>>("Flavors");
    if (flavors.empty())
      throw std::runtime_error("Flavors must be non-empty");
    for (const auto& flavor : flavors) {
      const int pdg = flavor_pdg(flavor);
      const auto duplicate = std::find_if(
          requestedFlavors.begin(), requestedFlavors.end(),
          [pdg](const auto& entry) { return entry.second == pdg; });
      if (duplicate != requestedFlavors.end())
        throw std::runtime_error("Duplicate flavor in Flavors: '" + flavor + "'");
      requestedFlavors.emplace_back(flavor, pdg);
    }
  } else if (manifest.has_key("Flavor")) {
    const std::string flavor = manifest.get<std::string>("Flavor");
    requestedFlavors.emplace_back(flavor, flavor_pdg(flavor));
  } else {
    requestedFlavors.emplace_back("", 0);
  }

  std::vector<double> defaultEnergyGrid;
  if (manifest.has_key("EnergyGrid"))
    defaultEnergyGrid = manifest.get<std::vector<double>>("EnergyGrid");

  fhicl::ParameterSet energyGrids;
  const bool haveFlavorEnergyGrids = manifest.has_key("EnergyGrids");
  if (haveFlavorEnergyGrids)
    energyGrids = manifest.get<fhicl::ParameterSet>("EnergyGrids");

  auto validate_energy_grid = [](const std::vector<double>& grid,
                                 const std::string& label) {
    if (grid.empty())
      throw std::runtime_error(label + " must be non-empty");
    if (!std::is_sorted(grid.begin(), grid.end()) ||
        std::adjacent_find(grid.begin(), grid.end()) != grid.end()) {
      throw std::runtime_error(label + " must be in strictly ascending order");
    }
    if (!std::all_of(grid.begin(), grid.end(),
                     [](double energy) { return std::isfinite(energy); })) {
      throw std::runtime_error(label + " must contain only finite values");
    }
  };

  auto energy_grid_for = [&](const std::string& flavor) {
    if (!flavor.empty() && haveFlavorEnergyGrids && energyGrids.has_key(flavor)) {
      auto grid = energyGrids.get<std::vector<double>>(flavor);
      validate_energy_grid(grid, "EnergyGrids." + flavor);
      return grid;
    }
    if (defaultEnergyGrid.empty()) {
      const std::string suffix = flavor.empty() ? "" : " for flavor '" + flavor + "'";
      throw std::runtime_error("Missing EnergyGrid" + suffix);
    }
    validate_energy_grid(defaultEnergyGrid, "EnergyGrid");
    return defaultEnergyGrid;
  };

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
  const double configuredEnuMin = manifest.get<double>("EnuMin", 0.4);
  const double configuredEnuMax = manifest.get<double>("EnuMax", 2.5);
  fEnuSnapTol = manifest.get<double>("EnuSnapTol", 5e-3); // 5 MeV
  if (!(std::isfinite(configuredEnuMin) && std::isfinite(configuredEnuMax) && configuredEnuMax >= configuredEnuMin))
    throw std::runtime_error("EnuMin/EnuMax must be finite and EnuMax>=EnuMin");
  if (!(std::isfinite(fEnuSnapTol) && fEnuSnapTol >= 0.0))
    throw std::runtime_error("EnuSnapTol must be finite and >= 0");

  std::string weightMapFile3D =
      manifest.get<std::string>("WeightMapFile3D", "");
  weightMapFile3D = systtools::expand_env_vars(weightMapFile3D);
  const bool useWeightMap3D = !weightMapFile3D.empty();

  const bool mapIsQ3xQ0 = useWeightMap3D || manifest.get<bool>("MapIsQ3xQ0", false);
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



  std::cout << "  WeightLimits   : [" << fWmin << ", " << fWmax << "]\n"
            << "  Enu snap tol   : " << fEnuSnapTol << " GeV\n";





  // Histogram sourcing:
  //  (A) Combined q3 x q0 x Enu histograms in WeightMapFile3D
  //  (B) Arrays np_files/nn_files with explicit histogram names: HistNameNP/HistNameNN
  //  (C) Flavor-selected external files under NuisFlatWeightMapBaseDir
  //  (D) Legacy auto-generated file paths based on Model and DataBaseDir

  std::string dataBaseDir =
      manifest.get<std::string>("DataBaseDir", "");
  dataBaseDir = systtools::expand_env_vars(dataBaseDir);

  std::string nuisFlatBaseDir =
      manifest.get<std::string>("NuisFlatWeightMapBaseDir", "");
  nuisFlatBaseDir = systtools::expand_env_vars(nuisFlatBaseDir);

  if (useWeightMap3D) {
    const std::vector<std::string> incompatibleKeys = {
      "NuisFlatWeightMapBaseDir", "DataBaseDir", "np_files", "nn_files",
      "EnergyGrid", "EnergyGrids", "HistNameNP", "HistNameNN"
    };
    for (const auto& key : incompatibleKeys) {
      if (manifest.has_key(key)) {
        throw std::runtime_error(
            "WeightMapFile3D cannot be combined with legacy manifest key '" +
            key + "'");
      }
    }
  }

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

  auto external_model_dir = [&]() -> std::string {
    if (model == "valencia")
      return "Valencia";
    if (model == "martini")
      return "Martini";
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

  const std::string hname_np = manifest.get<std::string>("HistNameNP", "");
  const std::string hname_nn = manifest.get<std::string>("HistNameNN", "");
  fFlavorResponses.clear();

  auto load_response_data = [&](int probePdg, const std::string& flavor,
                                std::vector<double> energyGrid,
                                const std::vector<std::string>& npFiles,
                                const std::vector<std::string>& nnFiles) {
    if (hname_np.empty() || hname_nn.empty())
      throw std::runtime_error("HistNameNP and HistNameNN are required for TH2D weight-map sources");

    if (npFiles.size() != energyGrid.size() || nnFiles.size() != energyGrid.size()) {
      throw std::runtime_error("np_files/nn_files sizes must match the EnergyGrid size" +
                               (flavor.empty() ? std::string{} :
                                " for flavor '" + flavor + "'"));
    }

    auto insertion = fFlavorResponses.emplace(probePdg, FlavorResponseData{});
    if (!insertion.second)
      throw std::runtime_error("More than one MEC response set configured for probe PDG " +
                               std::to_string(probePdg));

    auto& responseData = insertion.first->second;
    responseData.energyGrid = std::move(energyGrid);
    responseData.enuMin = std::max(configuredEnuMin, responseData.energyGrid.front());
    responseData.enuMax = std::min(configuredEnuMax, responseData.energyGrid.back());
    if (responseData.enuMax < responseData.enuMin) {
      throw std::runtime_error("Configured Enu window does not overlap the EnergyGrid" +
                               (flavor.empty() ? std::string{} :
                                " for flavor '" + flavor + "'"));
    }

    std::cout << "  Response maps: "
              << (flavor.empty() ? "legacy flavor-independent" : flavor)
              << " (probe PDG " << probePdg << ")\n"
              << "    EnergyGrid size: " << responseData.energyGrid.size() << "\n"
              << "    Enu window: [" << responseData.enuMin << ", "
              << responseData.enuMax << "] GeV\n";

    auto load_list = [&](const std::vector<std::string>& files,
                         Topo topo, const std::string& hname) {
      auto& calculators = responseData.calcs[topo];
      calculators.reserve(files.size());
      for (const auto& fname : files) {
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
        calculators.emplace_back(std::move(calc));
      }
    };

    load_list(npFiles, Topo::np, hname_np);
    load_list(nnFiles, Topo::nn, hname_nn);
  };

  if (useWeightMap3D) {
    if (model != "valencia" && model != "martini")
      throw std::runtime_error("WeightMapFile3D requires Model to be 'valencia' or 'martini'");

    TFile fin(weightMapFile3D.c_str(), "READ");
    if (!fin.IsOpen() || fin.IsZombie())
      throw std::runtime_error("Cannot open 3D MEC weight-map file: " + weightMapFile3D);

    std::vector<double> referenceXEdges;
    std::vector<double> referenceYEdges;
    std::vector<double> referenceZEdges;

    auto nearly_equal = [](double lhs, double rhs) {
      return std::fabs(lhs - rhs) <= 1e-10;
    };

    auto axis_edges = [](const TAxis& axis) {
      std::vector<double> edges;
      edges.reserve(axis.GetNbins() + 1);
      for (int ibin = 1; ibin <= axis.GetNbins() + 1; ++ibin)
        edges.push_back(axis.GetBinLowEdge(ibin));
      return edges;
    };

    auto require_matching_edges = [&](std::vector<double>& reference,
                                      const TAxis& axis,
                                      const std::string& histName,
                                      const std::string& axisName) {
      auto edges = axis_edges(axis);
      if (reference.empty()) {
        reference = std::move(edges);
        return;
      }
      if (reference.size() != edges.size())
        throw std::runtime_error("Histogram '" + histName +
                                 "' has mismatched " + axisName + " binning");
      for (size_t i = 0; i < edges.size(); ++i) {
        if (!nearly_equal(reference[i], edges[i]))
          throw std::runtime_error("Histogram '" + histName +
                                   "' has mismatched " + axisName + " binning");
      }
    };

    auto validate_3d_histogram = [&](const TH3D& hist,
                                     const std::string& histName) {
      if (hist.GetNbinsX() != 500 || hist.GetNbinsY() != 500 ||
          hist.GetNbinsZ() != 57) {
        throw std::runtime_error("Histogram '" + histName +
                                 "' must have dimensions 500 x 500 x 57");
      }
      if (std::string(hist.GetXaxis()->GetTitle()) != "q3 [GeV]" ||
          std::string(hist.GetYaxis()->GetTitle()) != "q0 [GeV]" ||
          std::string(hist.GetZaxis()->GetTitle()) != "Enu [GeV]") {
        throw std::runtime_error("Histogram '" + histName +
                                 "' must have q3 x q0 x Enu axes in GeV");
      }
      if (!nearly_equal(hist.GetXaxis()->GetXmin(), 0.0) ||
          !nearly_equal(hist.GetXaxis()->GetXmax(), 2.5) ||
          !nearly_equal(hist.GetYaxis()->GetXmin(), 0.0) ||
          !nearly_equal(hist.GetYaxis()->GetXmax(), 2.5)) {
        throw std::runtime_error("Histogram '" + histName +
                                 "' must span 0--2.5 GeV in q3 and q0");
      }
      for (int iz = 1; iz <= hist.GetNbinsZ(); ++iz) {
        const double expected = 0.20 + 0.05 * static_cast<double>(iz - 1);
        if (!nearly_equal(hist.GetZaxis()->GetBinCenter(iz), expected)) {
          throw std::runtime_error("Histogram '" + histName +
                                   "' must have Enu centers from 0.20 to 3.00 GeV in 0.05 GeV steps");
        }
      }
      require_matching_edges(referenceXEdges, *hist.GetXaxis(), histName, "q3");
      require_matching_edges(referenceYEdges, *hist.GetYaxis(), histName, "q0");
      require_matching_edges(referenceZEdges, *hist.GetZaxis(), histName, "Enu");
    };

    std::cout << "  Combined 3D weight-map file: " << weightMapFile3D << "\n";
    for (const auto& flavorEntry : requestedFlavors) {
      const std::string flavor = flavorEntry.first.empty() ? "numu" : flavorEntry.first;
      auto insertion = fFlavorResponses.emplace(flavorEntry.second, FlavorResponseData{});
      if (!insertion.second)
        throw std::runtime_error("More than one MEC response set configured for probe PDG " +
                                 std::to_string(flavorEntry.second));
      auto& responseData = insertion.first->second;

      auto load_3d_histogram = [&](const std::string& topology, Topo topo) {
        const std::string histName = model + "_" + flavor + "_" + topology;
        TObject* object = fin.Get(histName.c_str());
        TH3D* rawHist = dynamic_cast<TH3D*>(object);
        if (!rawHist) {
          const std::string actualType = object ? object->ClassName() : "missing";
          throw std::runtime_error("Expected TH3D histogram '" + histName +
                                   "' in " + weightMapFile3D + "; found " + actualType);
        }

        rawHist->SetDirectory(nullptr);
        std::unique_ptr<TH3D> hist(rawHist);
        validate_3d_histogram(*hist, histName);

        std::vector<double> energyGrid;
        energyGrid.reserve(hist->GetNbinsZ());
        for (int iz = 1; iz <= hist->GetNbinsZ(); ++iz)
          energyGrid.push_back(hist->GetZaxis()->GetBinCenter(iz));
        if (responseData.energyGrid.empty()) {
          responseData.energyGrid = std::move(energyGrid);
        } else if (responseData.energyGrid != energyGrid) {
          throw std::runtime_error("Histogram '" + histName +
                                   "' has a topology-dependent Enu grid");
        }

        auto calc = std::make_unique<MECq0q3ResponseCalc3D>(hist.get(), fWmin, fWmax);
        calc->SetUseNearestBin(useNearestBin);
        calc->SetEdgeClamp(edgeClamp);
        calc->SetOutOfRangeWeight(outOfRangeWeight);
        responseData.calcs3D.emplace(topo, std::move(calc));
        std::cout << "  Loaded " << histName << "  X:["
                  << hist->GetXaxis()->GetXmin() << ","
                  << hist->GetXaxis()->GetXmax() << "]  Y:["
                  << hist->GetYaxis()->GetXmin() << ","
                  << hist->GetYaxis()->GetXmax() << "]  Enu:["
                  << responseData.energyGrid.front() << ","
                  << responseData.energyGrid.back() << "] GeV\n";
      };

      load_3d_histogram("np", Topo::np);
      load_3d_histogram("nn", Topo::nn);
      responseData.enuMin = std::max(configuredEnuMin, responseData.energyGrid.front());
      responseData.enuMax = std::min(configuredEnuMax, responseData.energyGrid.back());
      if (responseData.enuMax < responseData.enuMin) {
        throw std::runtime_error("Configured Enu window does not overlap the 3D map grid for flavor '" +
                                 flavor + "'");
      }
    }
  } else if (manifest.has_key("np_files") || manifest.has_key("nn_files")) {
    if (!(manifest.has_key("np_files") && manifest.has_key("nn_files")))
      throw std::runtime_error("Both np_files and nn_files are required when either is specified");
    if (requestedFlavors.size() != 1)
      throw std::runtime_error("Flavors requires NuisFlatWeightMapBaseDir; explicit file arrays describe only one map set");

    const auto& flavorEntry = requestedFlavors.front();
    load_response_data(
        flavorEntry.second, flavorEntry.first, energy_grid_for(flavorEntry.first),
        manifest.get<std::vector<std::string>>("np_files"),
        manifest.get<std::vector<std::string>>("nn_files"));
  } else if (!nuisFlatBaseDir.empty()) {
    if (model.empty())
      throw std::runtime_error("NuisFlatWeightMapBaseDir requires Model to be 'valencia' or 'martini'");

    const std::string tune = external_model_tune();
    const std::string modelDir = external_model_dir();

    std::cout << "[MECq0q3InterpWeighting] Auto-selecting model: " << model << "\n";
    std::cout << "  External weight map base: " << nuisFlatBaseDir << "\n";

    for (const auto& flavorEntry : requestedFlavors) {
      // Legacy manifests without Flavor continue to use the numu maps under
      // key zero, so the same maps remain applicable to every probe flavor.
      const std::string flavor = flavorEntry.first.empty() ? "numu" : flavorEntry.first;
      const std::string flavorTag = flavor_tag(flavor);
      std::cout << "  Flavor: " << flavor << "\n";

      auto energyGrid = energy_grid_for(flavor);
      std::vector<std::string> npFiles;
      std::vector<std::string> nnFiles;
      npFiles.reserve(energyGrid.size());
      nnFiles.reserve(energyGrid.size());

      auto external_path = [&](const std::string& topo, double E) {
        return nuisFlatBaseDir + "/" + modelDir + "/" + flavor + "/" + topo +
               "/weight_map_ar23_to_" + tune + "_" + flavorTag + "_" +
               topo + "_" + energy_label(E) + "GeV.root";
      };

      for (double E : energyGrid) {
        const std::string np_file = external_path("np", E);
        const std::string nn_file = external_path("nn", E);
        require_existing_map(np_file, flavor, "np", E);
        require_existing_map(nn_file, flavor, "nn", E);
        npFiles.push_back(np_file);
        nnFiles.push_back(nn_file);
        std::cout << "  Generated np file: " << np_file << "\n";
        std::cout << "  Generated nn file: " << nn_file << "\n";
      }
      load_response_data(flavorEntry.second, flavorEntry.first,
                         std::move(energyGrid), npFiles, nnFiles);
    }
  } else if (!model.empty() && !dataBaseDir.empty()) {
    if (requestedFlavors.size() != 1)
      throw std::runtime_error("Flavors requires NuisFlatWeightMapBaseDir; legacy DataBaseDir maps are not flavor-specific");

    std::string modelDir, filePrefix;
    set_legacy_model_info(modelDir, filePrefix);
    const auto& flavorEntry = requestedFlavors.front();
    auto energyGrid = energy_grid_for(flavorEntry.first);
    std::vector<std::string> npFiles;
    std::vector<std::string> nnFiles;
    npFiles.reserve(energyGrid.size());
    nnFiles.reserve(energyGrid.size());

    std::cout << "[MECq0q3InterpWeighting] Auto-selecting model: " << model << "\n";
    std::cout << "  Model directory: " << modelDir << "\n";

    for (double E : energyGrid) {
      const std::string npFile = dataBaseDir + "/" + modelDir + "/" +
          filePrefix + "_np_" + energy_label(E) + "GeV.root";
      const std::string nnFile = dataBaseDir + "/" + modelDir + "/" +
          filePrefix + "_nn_" + energy_label(E) + "GeV.root";
      npFiles.push_back(npFile);
      nnFiles.push_back(nnFile);
      std::cout << "  Generated np file: " << npFile << "\n";
      std::cout << "  Generated nn file: " << nnFile << "\n";
    }
    load_response_data(flavorEntry.second, flavorEntry.first,
                       std::move(energyGrid), npFiles, nnFiles);
  } else {
    throw std::runtime_error("Need WeightMapFile3D, (np_files & nn_files), (NuisFlatWeightMapBaseDir & Model), or (Model & DataBaseDir)");
  }

  std::cout << "[MECq0q3InterpWeighting] SetupResponseCalculator done\n";
  return true;
}

// ---------------------------------------------------------------------------
// Event response
systtools::event_unit_response_t
MECq0q3InterpWeighting::GetEventResponse(genie::EventRecord const& ev)
{
  const auto* probe = ev.Probe();
  if (!probe)
    return this->GetDefaultEventResponse();

  auto flavorResponseIt = fFlavorResponses.find(probe->Pdg());
  if (flavorResponseIt == fFlavorResponses.end())
    flavorResponseIt = fFlavorResponses.find(0);
  if (flavorResponseIt == fFlavorResponses.end())
    return this->GetDefaultEventResponse();

  const auto& flavorResponse = flavorResponseIt->second;
  const auto& energyGrid = flavorResponse.energyGrid;

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

  // Energy guard: weight=1 outside this flavor's available map grid.
  if (Enu < flavorResponse.enuMin - 1e-6 || Enu > flavorResponse.enuMax + 1e-6)
    return this->GetDefaultEventResponse();

  double w_blend = 1.0;
  const auto calc3DIt = flavorResponse.calcs3D.find(topo);
  if (calc3DIt != flavorResponse.calcs3D.end()) {
    w_blend = calc3DIt->second->GetCentralWeight(q0, q3, Enu, fEnuSnapTol);
  } else {
    // Legacy TH2D mode: find and blend the two bracketing energy maps.
    auto it_hi = std::lower_bound(energyGrid.begin(), energyGrid.end(), Enu);
    size_t ih = (it_hi == energyGrid.end()) ? energyGrid.size() - 1
                                            : std::distance(energyGrid.begin(), it_hi);
    size_t il = (ih == 0) ? 0 : ih - 1;

    const double Elo = energyGrid[il];
    const double Ehi = energyGrid[ih];
    double t = 0.0;
    if (std::fabs(Enu - Elo) <= fEnuSnapTol) {
      ih = il;
    } else if (std::fabs(Enu - Ehi) <= fEnuSnapTol) {
      il = ih;
    } else {
      t = (ih == il || Ehi <= Elo)
          ? 0.0
          : std::clamp((Enu - Elo) / (Ehi - Elo), 0.0, 1.0);
    }

    const auto& vec = flavorResponse.calcs.at(topo);
    const double w_lo = vec[il]->GetCentralWeight(q0, q3);
    const double w_hi = vec[ih]->GetCentralWeight(q0, q3);
    w_blend = (1.0 - t) * w_lo + t * w_hi;
  }

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
