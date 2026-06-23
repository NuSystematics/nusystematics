// DumpDialResponseNuSyst.cxx
// For each configured dial, sample N events per interaction channel and
// plot the weight-vs-parameter-value curve per event. Reveals the shape of
// the dial response (linearity, monotonicity, asymmetry, event dependence).
//
// Uses IGENIESystProvider_tool::GetEventResponse(evt, paramVals) to force
// arbitrary grid evaluation via OverrideVariations.

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystParamHeader.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/response_helper.hh"
#include "nusystematics/utility/silence_genie.hh"

#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"

#include <cstdlib>
#include <unistd.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace systtools;
using namespace nusyst;

// ===== CLI =================================================================
namespace cliopts {
std::string fclname;
std::string input_file;
std::string output_base = "dial_response";
std::string genie_branch_name = "gmcrec";
std::vector<std::string> parameters;
size_t NMax = std::numeric_limits<size_t>::max();
int n_per_channel = 5;
int nsteps = 0; // 0 = use paramVariations as-is
bool make_pdf = true;
bool make_root = true;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Required:\n"
    "    -c <config.fcl>      FHiCL config file\n"
    "    -i <ghep.root>       GENIE ghep input\n\n"
    "  Optional:\n"
    "    -p <par1,par2,...>    Filter dials by name substring\n"
    "    -n <N>                Events per channel (default 5)\n"
    "    --nsteps <N>          Resample paramVariations range with N points\n"
    "                          (default: use paramVariations as-is)\n"
    "    -o <basename>         Output basename (default: dial_response)\n"
    "    -b <branch>           NtpMCEventRecord branch (default: gmcrec)\n"
    "    -N <NMax>             Max events to scan for sampling\n"
    "    --no-pdf / --no-root  Skip output\n"
    "    -?|--help             Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "-c") cliopts::fclname = argv[++opt];
    else if (s == "-i") cliopts::input_file = argv[++opt];
    else if (s == "-o") {
      // -o is a *base* name; the tool appends .pdf and .root. If the user
      // passed e.g. `-o foo.pdf`, strip the extension so we don't produce
      // `foo.pdf.pdf` / `foo.pdf.root`.
      std::string v = argv[++opt];
      for (auto const &ext : {std::string(".pdf"), std::string(".root")}) {
        if (v.size() >= ext.size() &&
            v.compare(v.size() - ext.size(), ext.size(), ext) == 0) {
          v.erase(v.size() - ext.size());
          break;
        }
      }
      cliopts::output_base = v;
    }
    else if (s == "-b") cliopts::genie_branch_name = argv[++opt];
    else if (s == "-n") cliopts::n_per_channel = str2T<int>(argv[++opt]);
    else if (s == "-N") cliopts::NMax = str2T<size_t>(argv[++opt]);
    else if (s == "--nsteps") cliopts::nsteps = str2T<int>(argv[++opt]);
    else if (s == "--no-pdf") cliopts::make_pdf = false;
    else if (s == "--no-root") cliopts::make_root = false;
    else if (s == "-p") {
      std::string tok; std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::parameters.push_back(tok);
    } else {
      std::cout << "[ERROR]: Unknown option: " << s << std::endl;
      SayUsage(argv); exit(1);
    }
    opt++;
  }
}

bool ParamSelected(const std::string &name) {
  if (cliopts::parameters.empty()) return true;
  for (auto &sel : cliopts::parameters)
    if (name.find(sel) != std::string::npos) return true;
  return false;
}

// ===== Channel classification (copied from PlotSystVariationsNuSyst) =======
std::string NuName(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "numu" : "numubar";
    case 12: return pdg > 0 ? "nue"  : "nuebar";
    case 16: return pdg > 0 ? "nutau" : "nutaubar";
    default: return Form("nu%d", pdg);
  }
}

std::string TargetName(int A) {
  switch (A) {
    case 40: return "Ar40";
    case 12: return "C12";
    case 16: return "O16";
    case 56: return "Fe56";
    case 1:  return "H";
    default: return Form("A%d", A);
  }
}

std::string InteractionType(const genie::EventRecord &ev) {
  auto const &proc = ev.Summary()->ProcInfo();
  if (proc.IsQuasiElastic())    return "QE";
  if (proc.IsMEC())             return "MEC";
  if (proc.IsResonant())        return "RES";
  if (proc.IsDeepInelastic())   return "DIS";
  return "Other";
}

std::string MakeChannelKey(const genie::EventRecord &ev) {
  std::string cc = ev.Summary()->ProcInfo().IsWeakCC() ? "CC" : "NC";
  int nu_pdg = ev.Probe() ? ev.Probe()->Pdg() : 0;
  int tgt_A = ev.TargetNucleus() ? ev.TargetNucleus()->A() : 1;
  return cc + "_" + NuName(nu_pdg) + "_" + TargetName(tgt_A) + "_" + InteractionType(ev);
}

// ===== Color palette =======================================================
void GetVarColor(int ivar, int &color, int &style) {
  static const int palette[] = {
    kBlue+2, kBlue, kAzure+7, kGreen+3, kSpring+4,
    kOrange+7, kRed, kRed+2, kMagenta+2, kViolet+1
  };
  static const int npal = sizeof(palette) / sizeof(palette[0]);
  color = palette[ivar % npal];
  style = 1 + (ivar / npal);
}

// ===== Data containers =====================================================
struct DialInfo {
  paramId_t pid;
  size_t provider_idx;
  std::string provider_name;
  std::string name;
  std::string fullname;
  double central;
  bool is_correction;
  bool is_natural;
  bool is_splineable;
  bool is_weight_syst;
  std::array<double, 2> oneSigmaShifts;     // {down, up}
  std::array<double, 2> validityRange;      // {min, max}
  std::vector<double> paramVariations;       // original from fhicl
  std::vector<double> grid;                  // x-values where we evaluate
};

struct EventRef {
  Long64_t index;
  std::string channel;
};

// Per (dial, event) → vector of weights aligned with dial.grid
using ResponseMap = std::map<std::string, std::map<std::string, std::vector<std::pair<Long64_t, std::vector<double>>>>>;
// ResponseMap[dial_fullname][channel] -> list of (event_index, weights)

// ===== Cache resolution + provider-level filter (shared with nusyst tweaks) =
// If `-c` is omitted while `-p` is given, resolve fclname against the cached
// kitchen-sink path: $NUSYST_INVENTORY_FCL, then $nusystematics_ROOT/fcl/
// nusyst_inventory.fcl, then $NUSYST/fcl/..., then /tmp/...
// Auto-generating via `nusyst config` on first use.
constexpr const char *kInventoryEnvVar = "NUSYST_INVENTORY_FCL";

std::string InventoryDefaultPath() {
#ifdef NUSYST_INSTALL_PREFIX
  return std::string(NUSYST_INSTALL_PREFIX) + "/fcl/nusyst_inventory.fcl";
#else
  for (char const *var : {"nusystematics_ROOT", "NUSYST"}) {
    char const *val = std::getenv(var);
    if (val && *val) return std::string(val) + "/fcl/nusyst_inventory.fcl";
  }
  return "/tmp/nusyst_inventory.fcl";
#endif
}

// Returns true if `name` matches any substring in cliopts::parameters, or if
// the filter is empty (no filter = keep everything).
bool DialMatchesPFilter(const std::string &name) {
  if (cliopts::parameters.empty()) return true;
  for (auto const &sub : cliopts::parameters)
    if (name.find(sub) != std::string::npos) return true;
  return false;
}

// ===== Main ================================================================
int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  // Cache fallback for the common `-p DialName -i events.root -o curves`
  // invocation. Same resolution + auto-generation as nusyst inventory / tweaks.
  if (cliopts::fclname.empty() && !cliopts::parameters.empty()) {
    char const *env = std::getenv(kInventoryEnvVar);
    cliopts::fclname = (env && *env) ? std::string(env) : InventoryDefaultPath();
    if (::access(cliopts::fclname.c_str(), R_OK) != 0) {
      std::cerr << "[INFO]: -p was given without -c; auto-generating "
                << cliopts::fclname << " via `nusyst config --mode all`.\n";
      std::string cmd = "GenerateAllDialsConfigNuSyst --mode all -o " +
                        cliopts::fclname + " > /dev/null 2>&1";
      if (std::system(cmd.c_str()) != 0) {
        std::cerr << "[ERROR]: Auto-generation failed; re-running visibly:\n";
        std::system(("GenerateAllDialsConfigNuSyst --mode all -o " +
                     cliopts::fclname).c_str());
        return 3;
      }
    }
  }

  if (cliopts::fclname.empty() || cliopts::input_file.empty()) {
    std::cout << "[ERROR]: -c and -i are required.\n"
                 "         (Pass -p <dial,...> to auto-resolve -c from the "
                 "cached kitchen sink.)" << std::endl;
    SayUsage(argv);
    return 1;
  }

  // Load providers -- silence GENIE banner + per-provider chatter.
  // Two-step load so a -p filter can drop providers whose dials all miss
  // before any provider is instantiated.
  nusyst::quiet::SetGlobalQuiet();
  response_helper phh;
  {
    nusyst::quiet::StdoutSink _quiet;
    genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
        "Messenger_whisper.xml");

    fhicl::ParameterSet raw_ps;
    {
      std::unique_ptr<cet::filepath_maker> fm =
          std::make_unique<cet::filepath_lookup_nonabsolute>("FHICL_FILE_PATH");
      raw_ps = fhicl::ParameterSet::make(cliopts::fclname, *fm);
    }
    fhicl::ParameterSet gen_ps = raw_ps.get<fhicl::ParameterSet>(
        "generated_systematic_provider_configuration");

    if (!cliopts::parameters.empty()) {
      auto provider_names =
          gen_ps.get<std::vector<std::string>>("syst_providers");
      std::vector<std::string> kept;
      for (auto const &pname : provider_names) {
        fhicl::ParameterSet prov;
        try { prov = gen_ps.get<fhicl::ParameterSet>(pname); }
        catch (...) { continue; }
        bool any_match = false;
        for (auto const &key : prov.get_names()) {
          if (!prov.is_key_to_table(key)) continue;
          try {
            fhicl::ParameterSet sub = prov.get<fhicl::ParameterSet>(key);
            if (!sub.has_key("prettyName")) continue;
            std::string pretty = sub.get<std::string>("prettyName");
            if (DialMatchesPFilter(pretty)) { any_match = true; break; }
          } catch (...) {}
        }
        if (any_match) kept.push_back(pname);
        else           gen_ps.erase(pname);
      }
      gen_ps.put_or_replace<std::vector<std::string>>("syst_providers", kept);
      std::cerr << "[INFO]: -p filter kept " << kept.size() << " of "
                << provider_names.size() << " providers ("
                << (provider_names.size() - kept.size())
                << " skipped -- neither constructed nor evaluated).\n";
    }
    phh.LoadProvidersAndHeaders(gen_ps);
  }

  // Collect dial info
  std::vector<DialInfo> dials;
  for (paramId_t pid : phh.GetParameters()) {
    SystParamHeader const &hdr = phh.GetHeader(pid);
    if (hdr.isResponselessParam) continue;
    if (!ParamSelected(hdr.prettyName)) continue;

    // Find provider
    size_t prov_idx = 0;
    std::string provname = "?";
    for (size_t si = 0; si < phh.GetSystProvider().size(); ++si) {
      if (phh.GetSystProvider()[si]->ParamIsHandled(pid)) {
        prov_idx = si;
        provname = phh.GetSystProvider()[si]->GetFullyQualifiedName();
        break;
      }
    }

    DialInfo di;
    di.pid = pid;
    di.provider_idx = prov_idx;
    di.provider_name = provname;
    di.name = hdr.prettyName;
    di.fullname = provname + "_" + hdr.prettyName;
    di.central = hdr.centralParamValue;
    di.is_correction = hdr.isCorrection;
    di.is_natural = hdr.unitsAreNatural;
    di.is_splineable = hdr.isSplineable;
    di.is_weight_syst = hdr.isWeightSystematicVariation;
    di.oneSigmaShifts = hdr.oneSigmaShifts;
    di.validityRange = hdr.paramValidityRange;
    di.paramVariations = hdr.paramVariations;

    if (hdr.paramVariations.empty()) {
      std::cout << "[INFO]: Skipping, " << di.fullname
                << " has no variations." << std::endl;
      continue;
    }

    // Build evaluation grid
    if (cliopts::nsteps > 0) {
      double vmin = *std::min_element(hdr.paramVariations.begin(), hdr.paramVariations.end());
      double vmax = *std::max_element(hdr.paramVariations.begin(), hdr.paramVariations.end());
      di.grid.reserve(cliopts::nsteps);
      for (int i = 0; i < cliopts::nsteps; ++i) {
        double t = (cliopts::nsteps == 1) ? 0.5 : double(i) / (cliopts::nsteps - 1);
        di.grid.push_back(vmin + t * (vmax - vmin));
      }
    } else {
      di.grid = hdr.paramVariations;
    }

    dials.push_back(std::move(di));
  }

  if (dials.empty()) {
    std::cout << "[ERROR]: No dials to process." << std::endl;
    return 2;
  }

  std::cout << "Configured " << dials.size() << " dials." << std::endl;

  // Open ghep input
  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(cliopts::input_file.c_str())) {
    std::cout << "[ERROR]: Failed to add " << cliopts::input_file << std::endl;
    return 3;
  }
  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  gevs->SetBranchAddress(cliopts::genie_branch_name.c_str(), &GenieNtpl);

  // Pass 1: scan events to pick N per channel
  Long64_t nentries = gevs->GetEntries(); if (cliopts::NMax < (size_t)nentries) nentries = (Long64_t)cliopts::NMax;
  std::map<std::string, std::vector<Long64_t>> channel_events;

  std::cout << "Scanning " << nentries << " events to select samples..." << std::endl;
  for (Long64_t ev = 0; ev < nentries; ++ev) {
    gevs->GetEntry(ev);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;
    std::string chkey = MakeChannelKey(GenieGHep);
    auto &vec = channel_events[chkey];
    if ((int)vec.size() < cliopts::n_per_channel) {
      vec.push_back(ev);
    }
    GenieNtpl->Clear();
  }

  size_t total_selected = 0;
  for (auto &[ch, evs] : channel_events) total_selected += evs.size();
  std::cout << "Selected " << total_selected << " events across "
            << channel_events.size() << " channels." << std::endl;

  // Pass 2: for each selected event, compute dial response curves
  // Storage: responses[dial_fullname][channel] -> list of (event_index, weights)
  ResponseMap responses;

  std::cout << "Computing dial responses..." << std::endl;
  size_t processed = 0;
  for (auto &[chkey, event_list] : channel_events) {
    for (Long64_t ev_idx : event_list) {
      gevs->GetEntry(ev_idx);
      genie::EventRecord const &GenieGHep = *GenieNtpl->event;

      for (auto &di : dials) {
        IGENIESystProvider_tool *prov = phh.GetSystProvider()[di.provider_idx].get();

        // Build paramVals for this single dial with its grid
        std::vector<std::pair<paramId_t, std::vector<double>>> paramVals = {
            {di.pid, di.grid}};

        std::vector<double> weights(di.grid.size(), 1.0);
        try {
          auto resp = prov->GetEventResponse(GenieGHep, paramVals);
          // Find our pid in the response
          for (auto const &pr : resp) {
            if (pr.pid == di.pid) {
              for (size_t i = 0; i < pr.responses.size() && i < weights.size(); ++i)
                weights[i] = pr.responses[i];
              break;
            }
          }
        } catch (std::exception &e) {
          // leave weights at 1
        }

        responses[di.fullname][chkey].push_back({ev_idx, std::move(weights)});
      }

      processed++;
      GenieNtpl->Clear();
    }
  }
  std::cout << "Processed " << processed << " events." << std::endl;

  // Write ROOT output
  if (cliopts::make_root) {
    std::string rootname = cliopts::output_base + ".root";
    TFile *fout = new TFile(rootname.c_str(), "RECREATE");
    for (auto &di : dials) {
      TDirectory *ddir = fout->mkdir(di.fullname.c_str());
      auto dit = responses.find(di.fullname);
      if (dit == responses.end()) continue;
      for (auto &[chkey, ev_list] : dit->second) {
        TDirectory *cdir = ddir->mkdir(chkey.c_str());
        cdir->cd();
        for (size_t i = 0; i < ev_list.size(); ++i) {
          auto &[ev_idx, w] = ev_list[i];
          TGraph *g = new TGraph((int)di.grid.size());
          for (size_t k = 0; k < di.grid.size(); ++k)
            g->SetPoint(k, di.grid[k], w[k]);
          g->SetName(Form("event_%lld", ev_idx));
          g->SetTitle(Form("%s | %s | event %lld;%s;weight",
                            di.fullname.c_str(), chkey.c_str(), ev_idx,
                            di.is_natural ? "param value (natural)" : "param value (sigma)"));
          g->Write();
        }
      }
    }
    fout->Close();
    delete fout;
    std::cout << "ROOT file written to " << rootname << std::endl;
  }

  // ---- Determine which (dial, channel) pairs are "active" --------------
  // Active = at least one event has at least one weight differing from 1
  // by more than the threshold. Inactive (dial, channel) pairs are skipped
  // in the PDF.
  const double activity_threshold = 1e-4;
  // active_channels[dial_fullname] = ordered list of channels that show variation
  std::map<std::string, std::vector<std::string>> active_channels;
  for (auto &di : dials) {
    auto dit = responses.find(di.fullname);
    if (dit == responses.end()) continue;
    for (auto &[chkey, ev_list] : dit->second) {
      bool active = false;
      for (auto &[idx, w] : ev_list) {
        for (double v : w) {
          if (std::isfinite(v) && std::abs(v - 1.0) > activity_threshold) {
            active = true; break;
          }
        }
        if (active) break;
      }
      if (active) active_channels[di.fullname].push_back(chkey);
    }
  }

  // Write PDF
  if (cliopts::make_pdf) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::string pdfname = cliopts::output_base + ".pdf";
    TCanvas *c = new TCanvas("c", "", 1200, 900);
    c->Print((pdfname + "[").c_str());

    auto print_with_bookmark = [&](const std::string &bookmark) {
      // ROOT TPDF uses "Title:..." in the Print options string to set a
      // PDF bookmark for the current page; PDF viewers expose this as a
      // clickable sidebar / outline pane for navigation.
      c->Print(pdfname.c_str(), Form("Title:%s", bookmark.c_str()));
    };

    // ---- Page 1: Summary --------------------------------------------------
    c->Clear();
    {
      TPad *sumpad = new TPad("sumpad", "", 0, 0, 1, 1);
      sumpad->SetFillColor(kWhite);
      sumpad->Draw(); sumpad->cd();
      TLatex tex;
      tex.SetNDC(); tex.SetTextFont(42);
      tex.SetTextSize(0.05); tex.SetTextFont(62);
      tex.DrawLatex(0.05, 0.93, "DumpDialResponseNuSyst Summary");
      tex.SetTextFont(42); tex.SetTextSize(0.028);
      double y = 0.86;
      tex.DrawLatex(0.05, y, Form("Config: %s", cliopts::fclname.c_str())); y -= 0.03;
      tex.DrawLatex(0.05, y, Form("Input:  %s", cliopts::input_file.c_str())); y -= 0.03;
      tex.DrawLatex(0.05, y, Form("Dials:  %zu   Active dials: %zu   Channels: %zu   Events/channel: %d",
                                   dials.size(), active_channels.size(),
                                   channel_events.size(), cliopts::n_per_channel));
      y -= 0.04;
      tex.SetTextFont(62);
      tex.DrawLatex(0.05, y, "Sampled events per channel:"); y -= 0.03;
      tex.SetTextFont(42);
      for (auto &[chkey, evs] : channel_events) {
        std::ostringstream ss;
        ss << chkey << ":  ";
        for (size_t i = 0; i < evs.size(); ++i) {
          if (i > 0) ss << ", ";
          ss << evs[i];
        }
        tex.DrawLatex(0.07, y, ss.str().c_str());
        y -= 0.025;
        if (y < 0.05) break;
      }
    }
    print_with_bookmark("Summary");

    // ---- Page 2+: Table of contents --------------------------------------
    // Track the page number each dial starts on. Page 1 = Summary; TOC
    // pages start at 2 and have unknown length, so first dial starts at
    // 2 + ceil(active_dials / lines_per_toc_page) + 1.
    const size_t lines_per_toc_page = 35;
    size_t n_active = active_channels.size();
    size_t n_toc_pages = (n_active + lines_per_toc_page - 1) / lines_per_toc_page;
    if (n_toc_pages == 0) n_toc_pages = 1;
    // First dial summary page starts after summary + TOC pages.
    size_t first_dial_page = 1 + n_toc_pages + 1; // 1 (summary) + n_toc + 1
    // But the index above counts the *next* page after the TOC, which is the
    // first dial-summary page. PDF page numbers are 1-based externally.

    // Compute dial → starting page map
    std::map<std::string, size_t> dial_page;
    {
      size_t cur = first_dial_page;
      for (auto &di : dials) {
        auto it = active_channels.find(di.fullname);
        if (it == active_channels.end()) continue;
        dial_page[di.fullname] = cur;
        cur += 1 /* metadata page */ + it->second.size();
      }
    }

    // Render TOC across n_toc_pages pages
    {
      std::vector<std::pair<std::string, size_t>> entries; // (dial fullname, page)
      for (auto &di : dials) {
        auto pit = dial_page.find(di.fullname);
        if (pit == dial_page.end()) continue;
        entries.emplace_back(di.fullname, pit->second);
      }
      size_t n_pages = (entries.size() + lines_per_toc_page - 1) / lines_per_toc_page;
      if (n_pages == 0) n_pages = 1;
      for (size_t pg = 0; pg < n_pages; ++pg) {
        c->Clear();
        TPad *tp = new TPad(Form("tocpad_%zu", pg), "", 0, 0, 1, 1);
        tp->SetFillColor(kWhite);
        tp->Draw(); tp->cd();
        TLatex tex;
        tex.SetNDC(); tex.SetTextFont(62); tex.SetTextSize(0.045);
        tex.DrawLatex(0.05, 0.93,
                       Form("Table of contents (%zu/%zu) -- %zu active dials",
                             pg + 1, n_pages, entries.size()));
        tex.SetTextFont(42); tex.SetTextSize(0.022);
        tex.DrawLatex(0.05, 0.89,
                       "Use your PDF viewer's outline / bookmark sidebar "
                       "for clickable navigation.");
        double y = 0.85;
        size_t i_start = pg * lines_per_toc_page;
        size_t i_end = std::min(i_start + lines_per_toc_page, entries.size());
        for (size_t i = i_start; i < i_end; ++i) {
          tex.DrawLatex(0.06, y,
                         Form("%-50s . . . . . . . . p. %zu",
                               entries[i].first.c_str(), entries[i].second));
          y -= 0.022;
        }
        std::string bm = (n_pages == 1) ? "Table of contents"
                                        : Form("Table of contents (%zu/%zu)", pg + 1, n_pages);
        print_with_bookmark(bm);
      }
    }

    // ---- For each active dial: metadata page + per-channel response pages
    for (auto &di : dials) {
      auto it = active_channels.find(di.fullname);
      if (it == active_channels.end()) continue;
      auto &chans = it->second;

      // Metadata page
      c->Clear();
      {
        TPad *p = new TPad("metap", "", 0, 0, 1, 1);
        p->SetFillColor(kWhite);
        p->Draw(); p->cd();
        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(62); tex.SetTextSize(0.05);
        tex.DrawLatex(0.05, 0.93, di.fullname.c_str());
        tex.SetTextFont(42); tex.SetTextSize(0.030);
        double y = 0.86;
        auto fmt_dbl = [](double v) -> std::string {
          if (v == kDefaultDouble) return "(unset)";
          std::ostringstream ss; ss << v; return ss.str();
        };
        auto fmt_vec = [](const std::vector<double> &v) -> std::string {
          if (v.empty()) return "(empty)";
          std::ostringstream ss;
          for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) ss << ", ";
            ss << v[i];
          }
          return ss.str();
        };
        auto line = [&](const std::string &k, const std::string &v) {
          tex.SetTextFont(62); tex.DrawLatex(0.06, y, Form("%s:", k.c_str()));
          tex.SetTextFont(42); tex.DrawLatex(0.32, y, v.c_str());
          y -= 0.038;
        };
        line("Provider", di.provider_name);
        line("paramId", std::to_string(di.pid));
        line("centralParamValue", fmt_dbl(di.central));
        line("paramVariations (fhicl)", fmt_vec(di.paramVariations));
        line("Eval grid", fmt_vec(di.grid));
        line("oneSigma down / up",
              fmt_dbl(di.oneSigmaShifts[0]) + " / " + fmt_dbl(di.oneSigmaShifts[1]));
        line("validity range",
              fmt_dbl(di.validityRange[0]) + " ... " + fmt_dbl(di.validityRange[1]));
        line("isSplineable", di.is_splineable ? "yes" : "no");
        line("isCorrection", di.is_correction ? "yes" : "no");
        line("unitsAreNatural", di.is_natural ? "yes (natural)" : "no (sigma)");
        line("isWeightSystematic", di.is_weight_syst ? "yes" : "no (property shift)");
        // Active channels
        tex.SetTextFont(62); tex.DrawLatex(0.06, y, "Acts on channels:");
        y -= 0.034; tex.SetTextFont(42);
        for (auto &ch : chans) {
          tex.DrawLatex(0.10, y, ch.c_str());
          y -= 0.026;
          if (y < 0.05) break;
        }
      }
      print_with_bookmark(di.fullname);

      // Response curve page per active channel
      auto dit = responses.find(di.fullname);
      for (auto &chkey : chans) {
        auto evlit = dit->second.find(chkey);
        if (evlit == dit->second.end()) continue;
        auto &ev_list = evlit->second;
        if (ev_list.empty()) continue;

        c->Clear();
        c->cd();

        TPad *tp = new TPad("tp", "", 0, 0.93, 1, 1.0);
        tp->SetFillColor(kWhite);
        tp->Draw(); tp->cd();
        TLatex title;
        title.SetNDC(); title.SetTextSize(0.45); title.SetTextAlign(22); title.SetTextFont(62);
        std::string ptitle = di.fullname + "  |  " + chkey;
        title.DrawLatex(0.5, 0.5, ptitle.c_str());

        c->cd();
        TPad *mp = new TPad("mp", "", 0.08, 0.04, 0.98, 0.92);
        mp->SetFillColor(kWhite);
        mp->SetLeftMargin(0.13); mp->SetBottomMargin(0.13);
        mp->SetTopMargin(0.05); mp->SetRightMargin(0.05);
        mp->Draw(); mp->cd();

        double wmin = 1, wmax = 1;
        for (auto &[idx, w] : ev_list) {
          for (double v : w) {
            if (std::isfinite(v)) {
              wmin = std::min(wmin, v);
              wmax = std::max(wmax, v);
            }
          }
        }
        double pad = 0.1 * (wmax - wmin);
        if (pad == 0) pad = 0.1;
        double gmin = *std::min_element(di.grid.begin(), di.grid.end());
        double gmax = *std::max_element(di.grid.begin(), di.grid.end());

        TH1D *frame = new TH1D(Form("frame_%s_%s", di.fullname.c_str(), chkey.c_str()),
                                 "", 100, gmin, gmax);
        frame->SetMinimum(wmin - pad);
        frame->SetMaximum(wmax + pad);
        std::string xaxis = di.is_natural ? "parameter value (natural units)"
                                           : "parameter value (#sigma)";
        frame->GetXaxis()->SetTitle(xaxis.c_str());
        frame->GetYaxis()->SetTitle("weight");
        frame->GetXaxis()->SetTitleSize(0.045);
        frame->GetYaxis()->SetTitleSize(0.045);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->Draw("AXIS");

        TLine *unity = new TLine(gmin, 1, gmax, 1);
        unity->SetLineStyle(2); unity->SetLineColor(kGray + 1);
        unity->Draw();

        if (di.central != kDefaultDouble && di.central >= gmin && di.central <= gmax) {
          TLine *cvline = new TLine(di.central, wmin - pad, di.central, wmax + pad);
          cvline->SetLineStyle(3); cvline->SetLineColor(kGray + 2);
          cvline->Draw();
        }

        TLegend *leg = new TLegend(0.65, 0.65, 0.94, 0.93);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->SetTextSize(0.032);
        leg->SetHeader("event index");

        for (size_t i = 0; i < ev_list.size(); ++i) {
          auto &[ev_idx, w] = ev_list[i];
          TGraph *g = new TGraph((int)di.grid.size());
          for (size_t k = 0; k < di.grid.size(); ++k)
            g->SetPoint(k, di.grid[k], w[k]);
          int col, sty;
          GetVarColor((int)i, col, sty);
          g->SetLineColor(col); g->SetLineWidth(2); g->SetLineStyle(sty);
          g->SetMarkerColor(col); g->SetMarkerStyle(20 + (int)(i % 8)); g->SetMarkerSize(0.9);
          g->Draw("LP SAME");
          leg->AddEntry(g, Form("%lld", ev_idx), "lp");
        }
        leg->Draw();

        print_with_bookmark(di.fullname + " | " + chkey);
      }
    }
    c->Print((pdfname + "]").c_str());
    std::cout << "PDF written to " << pdfname << std::endl;
    std::cout << "Active dials: " << active_channels.size() << " of " << dials.size() << std::endl;
    delete c;
  }

  delete gevs;
  return 0;
}
