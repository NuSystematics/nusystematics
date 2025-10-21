/*******************************************************************************
 * ValenciaMECq0q3InterpWeighting_tool.cc
 ******************************************************************************/
#include "ValenciaMECq0q3InterpWeighting_tool.hh"

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
#include <limits>

// SystematicsTools helper
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

// GENIE
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/InitialState.h"

using namespace nusyst;
using namespace systtools;

ValenciaMECq0q3InterpWeighting::ValenciaMECq0q3InterpWeighting(
    const fhicl::ParameterSet& p)
  : IGENIESystProvider_tool(p) {}

// ---------------------------------------------------------------------------
// Build metadata (standard NuSyst header parsing)
SystMetaData
ValenciaMECq0q3InterpWeighting::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                                  systtools::paramId_t firstId) {

  std::cout << "[ValenciaMECq0q3InterpWeighting::BuildSystMetaData] Called\n";

  SystMetaData smd;
  SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(ps,
                                           "ValenciaMECResponse",
                                           phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  // stash manifest for SetupResponseCalculator
  auto man = ps.get<fhicl::ParameterSet>("ValenciaMECResponse_input_manifest");
  tool_options.put("ValenciaMECResponse_input_manifest", man);

  return smd;
}

// ---------------------------------------------------------------------------
// Read manifest and build calculators
bool
ValenciaMECq0q3InterpWeighting::SetupResponseCalculator(fhicl::ParameterSet const &tool_opts)
{
  std::cout << "[ValenciaMECq0q3InterpWeighting] SetupResponseCalculator begin\n";

  const auto manifest =
      tool_opts.get<fhicl::ParameterSet>("ValenciaMECResponse_input_manifest");

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





  // q0 ramp guard (optional)
  fQ0GuardMin = manifest.get<double>("Q0GuardMin", 0.0);
  fQ0GuardMax = manifest.get<double>("Q0GuardMax", fQ0GuardMin);
  if (!(std::isfinite(fQ0GuardMin) && fQ0GuardMin >= 0.0))
    throw std::runtime_error("Q0GuardMin must be finite and >= 0");
  if (!(std::isfinite(fQ0GuardMax) && fQ0GuardMax >= fQ0GuardMin))
    throw std::runtime_error("Q0GuardMax must be finite and >= Q0GuardMin");

  


  // Energy window & snap tolerance (defaults implement your request)
  fEnuMin     = manifest.get<double>("EnuMin",     0.4);
  fEnuMax     = manifest.get<double>("EnuMax",     2.5);
  fEnuSnapTol = manifest.get<double>("EnuSnapTol", 5e-3); // 5 MeV
  if (!(std::isfinite(fEnuMin) && std::isfinite(fEnuMax) && fEnuMax >= fEnuMin))
    throw std::runtime_error("EnuMin/EnuMax must be finite and EnuMax>=EnuMin");
  if (!(std::isfinite(fEnuSnapTol) && fEnuSnapTol >= 0.0))
    throw std::runtime_error("EnuSnapTol must be finite and >= 0");

  const bool mapIsQ3xQ0 = manifest.get<bool>("MapIsQ3xQ0", false);
  std::cout << "  MapIsQ3xQ0     : " << (mapIsQ3xQ0 ? "true" : "false") << "\n";



  std::cout << "  EnergyGrid size: " << fEgrid.size() << "\n"
            << "  WeightLimits   : [" << fWmin << ", " << fWmax << "]\n"
            << "  Q0Guard Min/Max: " << fQ0GuardMin << " / " << fQ0GuardMax << " GeV\n"
            << "  Enu window     : [" << fEnuMin << ", " << fEnuMax << "] GeV\n"
            << "  Enu snap tol   : " << fEnuSnapTol << " GeV\n";





  // Histogram sourcing:
  //  (A) WeightFile with conventional names: h_weights_map_{np|nn}_<E>GeV
  //  (B) Arrays np_files/nn_files with explicit histogram names: HistNameNP/HistNameNN
  const bool haveSingle = manifest.has_key("WeightFile");
  const bool haveArrays = manifest.has_key("np_files") && manifest.has_key("nn_files");
  if (!haveSingle && !haveArrays)
    throw std::runtime_error("Need either WeightFile or (np_files & nn_files)");

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
        const std::string hname = Form("h_weights_map_%s_%0.1fGeV", tag, E);
        TH2D* h = dynamic_cast<TH2D*>(fin.Get(hname.c_str()));
        if (!h) throw std::runtime_error("Missing histogram '" + hname +
                                         "' in file " + fname);

        std::cout << "  Loaded " << fname << " :: " << hname
                  << "  X:[" << h->GetXaxis()->GetXmin() << "," << h->GetXaxis()->GetXmax() << "]"
                  << "  Y:[" << h->GetYaxis()->GetXmin() << "," << h->GetYaxis()->GetXmax() << "]\n";

        vec.emplace_back(std::make_unique<ValenciaMECq0q3ResponseCalc>(h, fWmin, fWmax, mapIsQ3xQ0));
      }
    }
    fin.Close();
  } else {
    // arrays: explicit TH2 names are REQUIRED
    auto np_files = manifest.get<std::vector<std::string>>("np_files");
    auto nn_files = manifest.get<std::vector<std::string>>("nn_files");
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

        vec.emplace_back(std::make_unique<ValenciaMECq0q3ResponseCalc>(h, fWmin, fWmax, mapIsQ3xQ0));
        fin.Close();
      }
    };

    load_list(np_files, Topo::np, hname_np);
    load_list(nn_files, Topo::nn, hname_nn);
  }

  std::cout << "[ValenciaMECq0q3InterpWeighting] SetupResponseCalculator done\n";
  return true;
}

// ---------------------------------------------------------------------------
// Event response
systtools::event_unit_response_t
ValenciaMECq0q3InterpWeighting::GetEventResponse(genie::EventRecord const& ev)
{
  // classify topology
  const Topo topo = ClassifyEvent(ev);
  if (topo == Topo::unknown)
    return this->GetDefaultEventResponse();

  // compute leptonic transfers
  double q0 = 0.0, q3 = 0.0, Enu = 0.0;
  ComputeQ0Q3(ev, q0, q3, Enu);

  // Energy guard: weight=1 outside [fEnuMin, fEnuMax]
  if (Enu < fEnuMin - 1e-6 || Enu > fEnuMax + 1e-6)
    return this->GetDefaultEventResponse();

  // NEW: q0 apply window gate: unity if outside (q0min, q0max)
  if (q0 <= fQ0ApplyMin + 1e-6 || q0 >= fQ0ApplyMax - 1e-6)
    return this->GetDefaultEventResponse();


  // q0 guard: unity below Min (inclusive with epsilon)
  if (q0 <= fQ0GuardMin + 1e-6)
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

  // ramp in q0: w_eff = 1 + alpha * (w_blend - 1)
  double alpha = 1.0;
  if (fQ0GuardMax > fQ0GuardMin && q0 < fQ0GuardMax)
    alpha = std::clamp((q0 - fQ0GuardMin) / (fQ0GuardMax - fQ0GuardMin), 0.0, 1.0);

  const double w_eff_cv = std::clamp(1.0 + alpha * (w_blend - 1.0), fWmin, fWmax);
  const double one_sigma = (std::clamp(w_blend, fWmin, fWmax) - 1.0); // for variations

  // Build response vector
  systtools::event_unit_response_t response;
  if (!this->GetSystMetaData().empty()) {
    const auto &hdr = this->GetSystMetaData()[0];
    systtools::ParamResponses pr;
    pr.pid = hdr.systParamId;
    pr.responses.reserve(hdr.paramVariations.size());
    for (double d : hdr.paramVariations) {
      const double rw = std::clamp(1.0 + d * alpha * one_sigma, fWmin, fWmax);
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

  return response;
}

/*******************************************************************************
 *  Helpers
 ******************************************************************************/

// q0, q3, Enu from (probe – final-state lepton)
void
ValenciaMECq0q3InterpWeighting::ComputeQ0Q3(genie::EventRecord const& ev,
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
ValenciaMECq0q3InterpWeighting::Topo
ValenciaMECq0q3InterpWeighting::ClassifyEvent(genie::EventRecord const& ev)
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
