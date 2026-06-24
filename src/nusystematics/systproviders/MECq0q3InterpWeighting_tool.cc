/*******************************************************************************
 * MECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "MECq0q3InterpWeighting_tool.hh"


#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2.h>
#include <TString.h>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

// SystematicsTools helper
#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"
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
    const YAML::Node& p)
  : IGENIESystProvider_tool(p) {}

// ---------------------------------------------------------------------------
// Build metadata (standard NuSyst header parsing)
SystMetaData
MECq0q3InterpWeighting::BuildSystMetaData(YAML::Node const &yamlnd,
                                                  systtools::paramId_t firstId) {

  std::cout << "[MECq0q3InterpWeighting::BuildSystMetaData] Called\n";

  SystMetaData smd;
  
  // Check if Q0 binning is enabled
  auto man = yamlnd["MECResponse_input_manifest"];
  std::vector<double> q0Bins;
  if (man && man["Q0Bins"]) {
    q0Bins = man["Q0Bins"].as<std::vector<double>>();
  }
  
  // If Q0 binning is enabled, create multiple dials (one per bin)
  if (!q0Bins.empty() && q0Bins.size() >= 2) {
    std::cout << "  Q0-binned mode: creating " << (q0Bins.size() - 1) << " dials\n";
    
    for (size_t i = 0; i < q0Bins.size() - 1; ++i) {
      std::string dialName = "MECResponse_q0bin" + std::to_string(i);
      SystParamHeader phdr;
      if (ParseYAMLToolConfigurationParameter(yamlnd, dialName, phdr, firstId)) {
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
    if (ParseYAMLToolConfigurationParameter(yamlnd, "MECResponse", phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  // stash manifest for SetupResponseCalculator
  tool_options["MECResponse_input_manifest"] = man;

  return smd;
}

// ---------------------------------------------------------------------------
// Read manifest and build calculators
bool
MECq0q3InterpWeighting::SetupResponseCalculator(YAML::Node const &tool_opts)
{
  std::cout << "[MECq0q3InterpWeighting] SetupResponseCalculator begin\n";

  const auto manifest =
      tool_opts["MECResponse_input_manifest"];

  // Energy grid (required)
  if (!manifest || !manifest["EnergyGrid"])
    throw std::runtime_error("Missing EnergyGrid");
  fEgrid = manifest["EnergyGrid"].as<std::vector<double>>();
  if (fEgrid.empty())
    throw std::runtime_error("EnergyGrid must be non-empty");

  // Weight clamp (optional)
  if (manifest["WeightLimits"]) {
    auto wl = manifest["WeightLimits"].as<std::vector<double>>();
    if (wl.size() >= 2) {
      fWmin = std::min(wl.front(), wl.back());
      fWmax = std::max(wl.front(), wl.back());
    }
  }
  if (!(std::isfinite(fWmin) && std::isfinite(fWmax) && fWmin >= 0.0 && fWmax >= fWmin))
    throw std::runtime_error("Invalid WeightLimits");


  // NEW: q0 apply window (optional, defaults to [0, +inf))
  fQ0ApplyMin = manifest["Q0ApplyMin"] ? manifest["Q0ApplyMin"].as<double>() : 0.0;
  if (manifest["Q0ApplyMax"]) {
    fQ0ApplyMax = manifest["Q0ApplyMax"].as<double>();
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
  fQ3ApplyMin = manifest["Q3ApplyMin"] ? manifest["Q3ApplyMin"].as<double>() : 0.0;
  if (manifest["Q3ApplyMax"]) {
    fQ3ApplyMax = manifest["Q3ApplyMax"].as<double>();
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
  if (manifest["Q0Bins"]) {
    fQ0Bins = manifest["Q0Bins"].as<std::vector<double>>();
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
    if (manifest["Q0SelectMin"] || manifest["Q0SelectMax"]) {
      std::cout << "  WARNING: Q0SelectMin/Max ignored when Q0Bins is set\n";
    }
  } else {
    // Single dial mode - legacy Q0 selection range still works
    fQ0Bins.clear();
  }

  // NEW: Q0 selection range (optional, defaults to disabled)
  fQ0SelectMin = manifest["Q0SelectMin"] ? manifest["Q0SelectMin"].as<double>() : 0.0;
  fQ0SelectMax = manifest["Q0SelectMax"] ? manifest["Q0SelectMax"].as<double>() : 0.0;
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
  fEnuMin     = manifest["EnuMin"] ? manifest["EnuMin"].as<double>() : 0.4;
  fEnuMax     = manifest["EnuMax"] ? manifest["EnuMax"].as<double>() : 2.5;
  fEnuSnapTol = manifest["EnuSnapTol"] ? manifest["EnuSnapTol"].as<double>() : 5e-3; // 5 MeV
  if (!(std::isfinite(fEnuMin) && std::isfinite(fEnuMax) && fEnuMax >= fEnuMin))
    throw std::runtime_error("EnuMin/EnuMax must be finite and EnuMax>=EnuMin");
  if (!(std::isfinite(fEnuSnapTol) && fEnuSnapTol >= 0.0))
    throw std::runtime_error("EnuSnapTol must be finite and >= 0");

  const bool mapIsQ3xQ0 = manifest["MapIsQ3xQ0"] ? manifest["MapIsQ3xQ0"].as<bool>() : false;
  const bool useNearestBin = manifest["UseNearestBin"] ? manifest["UseNearestBin"].as<bool>() : true;  // turn ON new behavior
  const bool edgeClamp     = manifest["EdgeClamp"] ? manifest["EdgeClamp"].as<bool>() : true;  // clamp OOR to edge bin
  
  // Read Model parameter early to determine default out-of-range behavior
  std::string model = manifest["Model"] ? manifest["Model"].as<std::string>() : "";
  
  // Out-of-range weight: default depends on model
  // Valencia: 0.0 (suppress beyond q3~1.2 GeV, q0~??? GeV)
  // Martini:  0.0 (suppress beyond q0~0.995 GeV)
  // Custom:   configurable via OutOfRangeWeight parameter
  double outOfRangeWeight = 1.0;  // backward compatible default
  if (manifest["OutOfRangeWeight"]) {
    outOfRangeWeight = manifest["OutOfRangeWeight"].as<double>();
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
  //  (C) Auto-generate file paths based on Model parameter
  
  // Read DataBaseDir for auto-generation (model already read earlier for outOfRangeWeight)
  std::string dataBaseDir = manifest["DataBaseDir"] ? manifest["DataBaseDir"].as<std::string>() : "";
  dataBaseDir = systtools::expand_env_vars(dataBaseDir);

  std::vector<std::string> np_files, nn_files;
  
  // If Model is specified, auto-generate file paths
  if (!model.empty() && !dataBaseDir.empty()) {
    std::string modelDir, filePrefix;
    
    if (model == "valencia") {
      modelDir = "ValenciaMECq0q3";
      filePrefix = "reweight_data_SuSAv2_to_valencia";
    } else if (model == "martini") {
      modelDir = "martini_2p2h_weights";
      filePrefix = "reweight_data_SuSAv2_to_martini";
    } else {
      throw std::runtime_error("Unknown Model: '" + model + "'. Expected 'valencia' or 'martini'");
    }
    
    std::cout << "[MECq0q3InterpWeighting] Auto-selecting model: " << model << "\n";
    std::cout << "  Model directory: " << modelDir << "\n";
    
    // Generate file paths for each energy point
    for (double E : fEgrid) {
  std::string np_file = dataBaseDir + "/" + modelDir + "/" + 
           filePrefix + "_np_" + Form("%0.2f", E) + "GeV.root";
  std::string nn_file = dataBaseDir + "/" + modelDir + "/" + 
           filePrefix + "_nn_" + Form("%0.2f", E) + "GeV.root";
      np_files.push_back(np_file);
      nn_files.push_back(nn_file);
      std::cout << "  Generated np file: " << np_file << "\n";
      std::cout << "  Generated nn file: " << nn_file << "\n";
    }
  } else {
    // Fallback to explicitly provided file lists
    if (manifest["np_files"] && manifest["nn_files"]) {
      np_files = manifest["np_files"].as<std::vector<std::string>>();
      nn_files = manifest["nn_files"].as<std::vector<std::string>>();
    }
  }
  
  const bool haveSingle = manifest["WeightFile"] ? true : false;
  const bool haveArrays = !np_files.empty() && !nn_files.empty();
  if (!haveSingle && !haveArrays)
    throw std::runtime_error("Need either WeightFile, (np_files & nn_files), or (Model & DataBaseDir)");

  fCalcs.clear();

  if (haveSingle) {
    const std::string fname = manifest["WeightFile"].as<std::string>();
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

    const std::string hname_np = manifest["HistNameNP"].as<std::string>();
    const std::string hname_nn = manifest["HistNameNN"].as<std::string>();

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

  // Helper lambda to create zero-weight response (suppress events)
  // Same as setting w_eff_cv = 0 or one_sigma = -1.0
  auto GetZeroWeightResponse = [this]() {
    auto const& smd = this->GetSystMetaData();
    systtools::event_unit_response_t resp;
    resp.reserve(smd.size());
    for(auto const& sph : smd) {
      if (sph.isCorrection) {
        const double this_rw = std::clamp(1.0 + (sph.centralParamValue) * (-1.), 0., fWmax);
        resp.push_back({sph.systParamId, std::vector<double>{this_rw}});
      } else {
        std::vector<double> arr_rw;
        arr_rw.reserve(sph.paramVariations.size());
        for (double d : sph.paramVariations) {
          const double this_rw = std::clamp(1.0 + d * (-1.), 0., fWmax);
          arr_rw.push_back(this_rw);
        }
        resp.push_back({sph.systParamId, arr_rw});
      }
    }
    return resp;
  };

  // Determine if Q0 selection/binning is enabled
  // Either Q0SelectMin/Max OR Q0Bins means we want weight=1 outside apply windows
  const bool q0SelectionEnabled = (fQ0SelectMax > 0.0);
  const bool q0BinningEnabled = (!fQ0Bins.empty());

  // NEW: q3/q0 apply window gate logic
  // If Q0 selection OR Q0 binning is ENABLED: return weight=1 if outside apply windows
  // If BOTH are DISABLED: return weight=0 if outside apply windows
  const bool outsideQ3Apply = (q3 <= fQ3ApplyMin + 1e-6 || q3 >= fQ3ApplyMax - 1e-6);
  const bool outsideQ0Apply = (q0 <= fQ0ApplyMin + 1e-6 || q0 >= fQ0ApplyMax - 1e-6);
  
  if (outsideQ3Apply || outsideQ0Apply) {
    if (q0SelectionEnabled || q0BinningEnabled) {
      // When Q0 selection or binning is enabled, return weight=1 outside apply region
      return this->GetDefaultEventResponse();
    } else {
      // When both are disabled, return weight=0 outside apply region
      return GetZeroWeightResponse();
    }
  }

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
