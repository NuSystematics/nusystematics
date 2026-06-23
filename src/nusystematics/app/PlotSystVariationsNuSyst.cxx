// PlotSystVariationsNuSyst.cxx
// Plots systematic variations for nusystematics parameters,
// separated by interaction channel (QE, MEC, RES, DIS, Other)
// and by CC/NC, neutrino species, and target.
// Produces differential cross section plots with ratio panels,
// a summary page, and a shared legend.

#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/KinVarUtils.hh"
#include "nusystematics/utility/channel_classification.hh"
#include "nusystematics/utility/response_helper.hh"
#include "nusystematics/utility/silence_genie.hh"

#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Conventions/Units.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObjString.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"
#include "TExec.h"
#include "TTreeFormula.h"

#include <algorithm>
#include <cstdlib>
#include <unistd.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace systtools;
using namespace nusyst;

// ===== Predefined variables ================================================
struct VarDef {
  std::string branch;     // unique key
  // Quantity label only (no units). Used as both the x-axis quantity name
  // and as the numerator in the auto-generated y-axis "d#sigma/d<label>".
  // Examples: "Q^{2}", "p_{p}^{lead}", "cos#theta_{lep}".
  std::string label;
  // Units appended in parentheses on the x-axis ("<label> (<units>)") and as
  // a per-unit denominator on the y-axis ("(cm^{2}/<units>/nucleus)"). Leave
  // empty for dimensionless quantities (counts, angles); the y-axis then
  // reads "(cm^{2}/nucleus)".
  std::string units;
  // Optional explicit override for the dsigma/d<x> numerator. If empty,
  // auto-generated as "d#sigma/d" + label.
  std::string diffLabel;
  double xmin, xmax;
  // If non-empty, evaluate this TTreeFormula expression per event instead of
  // reading a single branch via GetVarValue. Set from a -v spec like
  //   -v "name:Max$(sqrt(fsparticles_px*fsparticles_px+...)):0,2[,nbins]"
  // Formula vars work only in flat-tree mode (TTreeFormula needs a TTree).
  std::string formula;
  // Per-variable bin count override; 0 = use cliopts::nbins.
  int         nbins = 0;
};

// X-axis title built from label + optional units.
inline std::string BuildXAxisTitle(VarDef const &vd) {
  if (vd.units.empty()) return vd.label;
  return vd.label + " (" + vd.units + ")";
}

// Y-axis title for the diff-xsec panel. Uses explicit diffLabel if set, else
// auto-generates "d#sigma/d<label>". Appends a per-unit denominator that
// reflects how ScaleToDiffXsec normalises the histogram (per nucleus, in cm
// squared, per bin width in x-axis units).
inline std::string BuildYAxisTitle(VarDef const &vd) {
  std::string numerator = vd.diffLabel.empty()
                            ? std::string("d#sigma/d") + vd.label
                            : vd.diffLabel;
  std::string denom = vd.units.empty()
                        ? std::string("cm^{2}/nucleus")
                        : std::string("cm^{2}/") + vd.units + "/nucleus";
  return numerator + " (" + denom + ")";
}

// A 2D plotting pair, used in `--dim 2` mode. Each pair just names two
// existing 1D entries from the variable registry; binning, labels and
// units inherit from those.
struct Var2DPair {
  std::string key;     // unique identifier and figure tag, e.g. "q0_q3"
  std::string x_var;   // x-axis variable name (matches kPredefinedVars key)
  std::string y_var;   // y-axis variable name
};

// Hardcoded 2D defaults. Same role as kPredefinedVars but for pairs.
// Convention: x first, y second. q3 on x / q0 on y matches the canonical
// nuclear-response (q3, q0) plane orientation.
static const std::vector<Var2DPair> kPredefined2DPairs = {
  {"q3_q0", "q3", "q0"},
  {"Q2_W",  "Q2", "W"},
  {"x_Q2",  "x",  "Q2"},
  {"pL_pT", "pL", "pT"},
};

// Built-in variable registry. The runtime registry `g_vars_registry` starts
// from a copy of this map and may be extended / overridden / disabled by
// `--plot-config <file.fcl>` (see LoadPlotConfig).
static const std::map<std::string, VarDef> kPredefinedVars = {
  // key            branch                label              units      diffLabel (empty -> auto from label)   xmin xmax
  {"Enu",         {"Enu_true",          "E_{#nu}",           "GeV",     "",                                     0,  10}},
  {"q0",          {"q0",                "q_{0}",             "GeV",     "",                                     0,   5}},
  {"Q2",          {"Q2",                "Q^{2}",             "GeV^{2}", "",                                     0,   5}},
  {"q3",          {"q3",                "|#vec{q}|",         "GeV",     "",                                     0,   5}},
  {"W",           {"w",                 "W",                 "GeV",     "",                                     0,   3}},
  {"plep",        {"plep",              "p_{lep}",           "GeV",     "",                                     0,   5}},
  {"coslep",      {"coslep",            "cos#theta_{lep}",   "",        "",                                    -1,   1}},
  {"EAvail",      {"EAvail_GeV",        "E_{avail}",         "GeV",     "",                                     0,   3}},
  {"Tp",          {"leading_proton_KE", "T_{p}",             "GeV",     "",                                     0, 1.5}},
  {"pLeadProton", {"leading_proton_p",  "p_{p}^{lead}",      "GeV",     "",                                     0, 2.0}},
  {"x",           {"Bjorken_x",         "x",                 "",        "",                                     0, 2.0}},
  {"pL",          {"plep_L",            "p_{L}^{lep}",       "GeV",     "",                                     0, 5.0}},
  {"pT",          {"plep_T",            "p_{T}^{lep}",       "GeV",     "",                                     0, 2.0}},
  {"Ereco",       {"Ereco_cal",         "E_{reco}^{cal}",    "GeV",     "",                                     0,  10}},
  {"nproton",     {"nproton",           "N_{proton}",        "",        "",                                     0,   6}},
  {"npip",        {"npip",              "N_{#pi^{+}}",       "",        "",                                     0,   5}},
  {"npim",        {"npim",              "N_{#pi^{-}}",       "",        "",                                     0,   5}},
  {"npi0",        {"npi0",              "N_{#pi^{0}}",       "",        "",                                     0,   5}},
  {"nneutron",    {"nneutron",          "N_{neutron}",       "",        "",                                     0,   6}},
};

// Forward-declared cliopts fields written by LoadPlotConfig (which is
// defined earlier in the file than the cliopts namespace block).
namespace cliopts {
  extern int plots_per_page;
  extern double ratio_min_cv_count;
}

// Runtime variable registry. Initialised lazily from kPredefinedVars and
// extended / overridden by --plot-config. Names with `enabled: false` in the
// config are erased from this map.
std::map<std::string, VarDef> g_vars_registry;

void EnsureRegistryInitialised() {
  if (g_vars_registry.empty())
    g_vars_registry = kPredefinedVars;
}

void PrintAvailableVars(std::ostream &os) {
  EnsureRegistryInitialised();
  os << "Available -v variables (case-sensitive):\n";
  for (auto const &[key, vd] : g_vars_registry) {
    os << "  " << std::left << std::setw(14) << key
       << "  range=[" << vd.xmin << ", " << vd.xmax << "]"
       << (vd.formula.empty()
              ? std::string("  branch=") + vd.branch
              : std::string("  formula"))
       << "\n";
  }
  os << "\nInline formula form (one-off, no config file):\n"
     << "  -v \"<name>:<TTreeFormula expr>:<xmin>,<xmax>[,<nbins>]\"\n";
}

// Parse "name:formula:xmin,xmax[,nbins]" into a VarDef. Splits on the FIRST
// colon (between name and the rest) and the LAST colon (between formula and
// the range triple) so that the formula itself may contain a colon (ternary
// `a?b:c`). Returns true if `spec` looks like the formula form.
bool TryParseInlineFormulaSpec(const std::string &spec, VarDef &out_vd) {
  size_t first = spec.find(':');
  size_t last  = spec.rfind(':');
  if (first == std::string::npos || first == last) return false;

  std::string name    = spec.substr(0, first);
  std::string formula = spec.substr(first + 1, last - first - 1);
  std::string range   = spec.substr(last + 1);

  std::vector<std::string> parts;
  std::istringstream ss(range);
  std::string tok;
  while (std::getline(ss, tok, ',')) parts.push_back(tok);
  if (parts.size() < 2) return false;

  try {
    out_vd = VarDef{};
    out_vd.branch     = name;     // unique key / tag
    out_vd.label      = name;     // no nice TLatex form available inline
    out_vd.units      = "";       // unknown from inline spec; use config for units
    out_vd.diffLabel  = "";       // auto from label
    out_vd.xmin       = std::stod(parts[0]);
    out_vd.xmax       = std::stod(parts[1]);
    out_vd.formula    = formula;
    out_vd.nbins      = (parts.size() >= 3) ? std::stoi(parts[2]) : 0;
  } catch (std::exception const &) {
    return false;
  }
  return true;
}

VarDef ResolveVar(const std::string &spec) {
  EnsureRegistryInitialised();

  // Inline formula form (contains a colon)
  VarDef vd;
  if (TryParseInlineFormulaSpec(spec, vd)) return vd;

  // Lookup in registry (built-in + config overrides)
  auto it = g_vars_registry.find(spec);
  if (it != g_vars_registry.end()) return it->second;

  std::cerr << "[ERROR]: Unknown variable '" << spec << "'.\n\n";
  PrintAvailableVars(std::cerr);
  std::cerr << "\nTo add new variables, write a plot-config fhicl and pass it\n"
               "via --plot-config <file.fcl>.  See doc/PlotConfig.example.fcl\n"
               "for the schema." << std::endl;
  std::exit(4);
}

// Load a plot-config fhicl file and merge its `plot_config` table into the
// global variable registry. Entries with `enabled: false` are erased. Entries
// whose names match an existing built-in override that built-in's fields
// (binning / labels / formula); fields omitted in the config inherit from the
// built-in.
void LoadPlotConfig(const std::string &path) {
  EnsureRegistryInitialised();

  // Try the literal path first (handles absolute paths and paths relative
  // to cwd -- what most users will type), fall back to FHICL_FILE_PATH lookup
  // (handles bare basenames like `MyPlotConfig.fcl` resolved from the search
  // path). Without the literal-first attempt, fhicl's nonabsolute lookup
  // refuses a path like `doc/Plot.fcl` because it interprets the slash as
  // "must be findable verbatim in FHICL_FILE_PATH".
  fhicl::ParameterSet raw;
  std::string literal_err, lookup_err;
  try {
    cet::filepath_maker fm;
    raw = fhicl::ParameterSet::make(path, fm);
  } catch (std::exception const &e) {
    literal_err = e.what();
    try {
      cet::filepath_lookup_nonabsolute fm("FHICL_FILE_PATH");
      raw = fhicl::ParameterSet::make(path, fm);
    } catch (std::exception const &e2) {
      lookup_err = e2.what();
    }
  }
  if (raw.is_empty()) {
    std::cerr << "[ERROR]: Failed to parse plot-config fhicl '" << path << "'.\n"
              << "         As literal path: " << literal_err << "\n"
              << "         Via FHICL_FILE_PATH: " << lookup_err << std::endl;
    std::exit(5);
  }

  fhicl::ParameterSet plot_cfg;
  try { plot_cfg = raw.get<fhicl::ParameterSet>("plot_config"); }
  catch (std::exception const &) {
    std::cerr << "[ERROR]: '" << path
              << "' has no top-level 'plot_config' table." << std::endl;
    std::exit(5);
  }

  // Top-level scalar knobs in plot_config (not per-variable settings).
  // The per-variable loop below skips these because is_key_to_table is false.
  cliopts::plots_per_page      = plot_cfg.get<int>("plots_per_page",
                                                    cliopts::plots_per_page);
  cliopts::ratio_min_cv_count  = plot_cfg.get<double>("ratio_min_cv_count",
                                                       cliopts::ratio_min_cv_count);

  for (auto const &name : plot_cfg.get_names()) {
    if (!plot_cfg.is_key_to_table(name)) continue;
    fhicl::ParameterSet entry;
    try { entry = plot_cfg.get<fhicl::ParameterSet>(name); }
    catch (std::exception const &) { continue; }

    if (!entry.get<bool>("enabled", true)) {
      g_vars_registry.erase(name);
      continue;
    }

    // Start from existing entry if overriding a built-in, else a blank VarDef
    // seeded with sensible defaults.
    auto existing = g_vars_registry.find(name);
    VarDef vd = (existing != g_vars_registry.end()) ? existing->second : VarDef{};
    if (existing == g_vars_registry.end()) {
      vd.branch    = name;
      vd.label     = name;  // user should override with a TLatex-clean version
      vd.units     = "";
      vd.diffLabel = "";    // auto = "d#sigma/d" + label
    }
    vd.branch    = entry.get<std::string>("branch", vd.branch);
    vd.label     = entry.get<std::string>("label", vd.label);
    vd.units     = entry.get<std::string>("units", vd.units);
    vd.diffLabel = entry.get<std::string>("diff_label", vd.diffLabel);
    vd.xmin      = entry.get<double>("xmin", vd.xmin);
    vd.xmax      = entry.get<double>("xmax", vd.xmax);
    vd.nbins     = entry.get<int>("nbins", vd.nbins);
    vd.formula   = entry.get<std::string>("formula", vd.formula);
    g_vars_registry[name] = vd;
  }

  std::cerr << "[INFO]: plot-config '" << path
            << "' applied -- registry has " << g_vars_registry.size()
            << " variables." << std::endl;
}

// ===== Channel classification ==============================================
struct EventVars {
  double Enu, q0, Q2, q3, w, plep, coslep, EAvail, leading_proton_KE, xsec;
  // Magnitude of highest-|p| proton's 3-momentum; -999 if no FS proton.
  double leading_proton_p;
  // Bjorken x and lepton momentum decomposition (along / perpendicular to
  // the incoming neutrino direction). Used as 1D vars and as axes for the
  // built-in 2D pairs.
  double Bjorken_x;
  double plep_L;
  double plep_T;
  double Emiss, pmiss, Ereco_cal;
  int nproton, npip, npim, npi0, nneutron;
  bool is_cc, is_qe, is_mec, is_res, is_dis;
  int nu_pdg, tgt_A, tgt_Z;
  bool has_coslep, has_leading_proton;
};

std::string NuName(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "numu" : "numubar";
    case 12: return pdg > 0 ? "nue"  : "nuebar";
    case 16: return pdg > 0 ? "nutau" : "nutaubar";
    default: return Form("nu%d", pdg);
  }
}

std::string NuLabel(int pdg) {
  switch (std::abs(pdg)) {
    case 14: return pdg > 0 ? "#nu_{#mu}" : "#bar{#nu}_{#mu}";
    case 12: return pdg > 0 ? "#nu_{e}"   : "#bar{#nu}_{e}";
    case 16: return pdg > 0 ? "#nu_{#tau}": "#bar{#nu}_{#tau}";
    default: return Form("#nu_{%d}", pdg);
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

std::string TargetLabel(int A) {
  switch (A) {
    case 40: return "^{40}Ar";
    case 12: return "^{12}C";
    case 16: return "^{16}O";
    case 56: return "^{56}Fe";
    case 1:  return "H";
    default: return Form("A=%d", A);
  }
}

std::string InteractionType(const EventVars &ev) {
  if (ev.is_qe)  return "QE";
  if (ev.is_mec) return "MEC";
  if (ev.is_res) return "RES";
  if (ev.is_dis) return "DIS";
  return "Other";
}

std::string InteractionLabel(const std::string &type) {
  if (type == "MEC") return "2p2h";
  return type;
}

// Sort key: CC before NC, numu before nue, QE/MEC/RES/DIS/Other order
std::string ChannelSortKey(const std::string &channel) {
  std::string key;
  key += (channel.find("CC_") == 0) ? "0" : "1";
  if (channel.find("numu_")    != std::string::npos) key += "0";
  else if (channel.find("nue_") != std::string::npos) key += "1";
  else if (channel.find("nutau_") != std::string::npos) key += "2";
  else key += "3";
  if (channel.find("bar_") != std::string::npos) key += "1"; else key += "0";
  if (channel.find("_QE") != std::string::npos) key += "0";
  else if (channel.find("_MEC") != std::string::npos) key += "1";
  else if (channel.find("_RES") != std::string::npos) key += "2";
  else if (channel.find("_DIS") != std::string::npos) key += "3";
  else key += "4";
  return key + "_" + channel;
}

std::string MakeChannelKey(const EventVars &ev) {
  std::string cc = ev.is_cc ? "CC" : "NC";
  return cc + "_" + NuName(ev.nu_pdg) + "_" + TargetName(ev.tgt_A) + "_" + InteractionType(ev);
}

std::string MakeChannelLabel(const std::string &key) {
  if (key == "Total") return "Total (all channels)";
  if (key == "Total_CC") return "Total CC";
  if (key == "Total_NC") return "Total NC";
  // Parse key: CC_numu_Ar40_QE
  std::istringstream ss(key);
  std::string cc, nu, tgt, inttype;
  std::getline(ss, cc, '_');
  // nu name may contain "bar" so read until target
  std::string rest;
  std::getline(ss, rest);
  // Find last two underscore-separated tokens as target and interaction
  size_t last_us = rest.rfind('_');
  inttype = rest.substr(last_us + 1);
  rest = rest.substr(0, last_us);
  size_t second_last = rest.rfind('_');
  tgt = rest.substr(second_last + 1);
  std::string nuname = rest.substr(0, second_last);

  // Build label
  int pdg = 0;
  if (nuname == "numu") pdg = 14;
  else if (nuname == "numubar") pdg = -14;
  else if (nuname == "nue") pdg = 12;
  else if (nuname == "nuebar") pdg = -12;
  else if (nuname == "nutau") pdg = 16;
  else if (nuname == "nutaubar") pdg = -16;

  int A = 0;
  if (tgt == "Ar40") A = 40;
  else if (tgt == "C12") A = 12;
  else if (tgt == "O16") A = 16;
  else if (tgt == "Fe56") A = 56;
  else if (tgt == "H") A = 1;
  else sscanf(tgt.c_str(), "A%d", &A);

  return cc + " " + (pdg ? NuLabel(pdg) : nuname) + " " +
         (A ? TargetLabel(A) : tgt) + " " + InteractionLabel(inttype);
}

// ===== Histogram bookkeeping ===============================================
struct SystHists {
  TH1D *h_cv = nullptr;
  std::vector<TH1D *> h_var;
  std::vector<double> tweakvals;
};

// 2D analogue of SystHists. Same per-(channel, dial) ownership pattern; the
// outer key in Hist2DMap is the Var2DPair::key (e.g. "q0_q3") rather than
// the 1D var.branch.
struct Syst2DHists {
  TH2D *h_cv = nullptr;
  std::vector<TH2D *> h_var;
  std::vector<double> tweakvals;
};

struct ParamMeta {
  std::string fullname;
  int ntweaks;
  std::vector<double> tweakvalues;
};

// ===== CLI =================================================================
namespace cliopts {
std::string input_file;
std::string fclname;
std::string fhicl_key = "generated_systematic_provider_configuration";
std::string output_base = "syst_variations";
std::string treename = "events";
std::string genie_branch_name = "gmcrec";
std::vector<std::string> variables;
std::vector<std::string> parameters;
size_t NMax = std::numeric_limits<size_t>::max();
int nbins = 50;
bool make_pdf = true;
bool make_root = true;
// 1 = 1D differential xsec mode (default); 2 = 2D pair mode with
// -1sigma / CV / +1sigma per row. Higher dimensions are not supported.
int dim = 1;
// Page layout: how many (variable, channel) panels per PDF page. 6 keeps the
// 3×2 default grid that fits 8.5×11 nicely; overridable via plot_config.
int plots_per_page = 6;
// Ratio-panel low-stats suppression: any (var/CV - 1) bin whose CV count is
// below this many events is blanked out so a near-zero denominator doesn't
// produce a spike that dominates the panel's y-axis range. NIWG calls the
// equivalent knob `min_event_rate`. Default 0 = no suppression; overridable
// via plot_config `ratio_min_cv_count`.
double ratio_min_cv_count = 0.0;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Required:\n"
    "    -i <file.root>       Input ROOT file (ghep or flat tree)\n\n"
    "  Optional:\n"
    "    -c <config.fcl>      FHiCL config (required for GHEP mode; if given\n"
    "                         in flat-tree mode, applies_to_channels patterns\n"
    "                         are read from it to skip dial x channel pages\n"
    "                         where the dial does not apply)\n"
    "    -k <fhicl_key>       Top-level fhicl key (default: generated_systematic_provider_configuration)\n"
    "    -o <output>          Output basename (default: syst_variations)\n"
    "    -v <var1,var2,...>    Variables to plot (default: all predefined).\n"
    "                         Pass --list-vars to see what's available.\n"
    "    -p <par1,par2,...>    Systematic parameters to plot (default: all)\n"
    "    -t <treename>        Tree name (default: events)\n"
    "    -b <branch>          NtpMCEventRecord branch (default: gmcrec)\n"
    "    -N <NMax>            Max events to process\n"
    "    --nbins <n>          Number of bins (default: 50)\n"
    "    --no-pdf / --no-root Skip output\n"
    "    --list-vars          Print the table of available -v variables and exit\n"
    "    --plot-config <f>    fhicl file overriding / adding plot variables.\n"
    "                         See config/PlotConfig.example.fcl for the schema.\n"
    "    --dim <1|2>          Plot dimensionality. Default 1 (differential xsec\n"
    "                         per kinematic variable). 2 produces 2D pair plots:\n"
    "                         per (channel, dial) page, two rows of three\n"
    "                         panels (-1sigma ratio | CV 2D | +1sigma ratio).\n"
    "                         Default pairs: q0:q3, Q2:W, x:Q2, pL:pT.\n"
    "    -?|--help            Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "--list-vars") {
      PrintAvailableVars(std::cout);
      exit(0);
    }
    else if (s == "--plot-config") {
      // Load immediately so subsequent --list-vars sees the augmented set,
      // and so -v lookups can hit config-only entries.
      LoadPlotConfig(argv[++opt]);
    }
    else if (s == "--dim") {
      cliopts::dim = str2T<int>(argv[++opt]);
      if (cliopts::dim != 1 && cliopts::dim != 2) {
        std::cerr << "[ERROR]: --dim must be 1 or 2 (got " << cliopts::dim
                  << ")" << std::endl;
        std::exit(2);
      }
    }
    else if (s == "-i") cliopts::input_file = argv[++opt];
    else if (s == "-c") cliopts::fclname = argv[++opt];
    else if (s == "-k") cliopts::fhicl_key = argv[++opt];
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
    else if (s == "-t") cliopts::treename = argv[++opt];
    else if (s == "-b") cliopts::genie_branch_name = argv[++opt];
    else if (s == "-N") cliopts::NMax = str2T<size_t>(argv[++opt]);
    else if (s == "--nbins") cliopts::nbins = str2T<int>(argv[++opt]);
    else if (s == "--no-pdf") cliopts::make_pdf = false;
    else if (s == "--no-root") cliopts::make_root = false;
    else if (s == "-v") {
      // Two accepted forms:
      //   1. Plain comma-separated list of variable names: -v Q2,W,Enu
      //   2. Single inline-formula spec:                  -v "name:Form:xmin,xmax[,nbins]"
      // Form 2 contains commas inside the range triple; splitting the whole
      // arg on commas would shred it. Detect form 2 by the presence of a
      // colon and treat the whole arg as a single spec. Multiple inline
      // formulas: pass -v multiple times.
      std::string arg = argv[++opt];
      if (arg.find(':') != std::string::npos) {
        cliopts::variables.push_back(arg);
      } else {
        std::string tok; std::istringstream ss(arg);
        while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::variables.push_back(tok);
      }
    } else if (s == "-p") {
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

// ===== Event counting for summary ==========================================
struct EventCounts {
  size_t total = 0;
  std::map<int, size_t> by_nu_pdg;
  std::map<int, size_t> by_tgt_A;
  std::map<std::string, size_t> by_channel; // full channel key
  std::map<std::string, size_t> by_current; // "CC" or "NC"
  std::map<std::string, size_t> by_inttype; // "CC_QE", "NC_RES", etc.

  void Count(const EventVars &ev) {
    total++;
    by_nu_pdg[ev.nu_pdg]++;
    by_tgt_A[ev.tgt_A]++;
    std::string cc = ev.is_cc ? "CC" : "NC";
    by_current[cc]++;
    by_inttype[cc + "_" + InteractionType(ev)]++;
    by_channel[MakeChannelKey(ev)]++;
  }
};

// ===== Helpers =============================================================
std::pair<double, double> AutoRange(TTree *tree, const std::string &branch) {
  double val;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus(branch.c_str(), 1);
  if (!tree->GetBranch(branch.c_str())) { tree->SetBranchStatus("*", 1); return {0, 1}; }
  tree->SetBranchAddress(branch.c_str(), &val);
  double vmin = 1e30, vmax = -1e30;
  Long64_t n = std::min((Long64_t)10000, tree->GetEntries());
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    if (std::isfinite(val)) { vmin = std::min(vmin, val); vmax = std::max(vmax, val); }
  }
  tree->SetBranchStatus("*", 1);
  tree->ResetBranchAddresses();
  if (vmin >= vmax) return {0, 1};
  double pad = 0.05 * (vmax - vmin);
  return {vmin - pad, vmax + pad};
}

void BookHistograms(const std::string &channel, const ParamMeta &pm,
                    const std::vector<VarDef> &vars,
                    std::map<std::string, SystHists> &var_hists,
                    TTree *tree_for_autorange = nullptr) {
  for (auto &vd : vars) {
    double xmin = vd.xmin, xmax = vd.xmax;
    if (xmin == 0 && xmax == 0 && tree_for_autorange && vd.formula.empty()) {
      auto r = AutoRange(tree_for_autorange, vd.branch);
      xmin = r.first; xmax = r.second;
    }
    if (xmin == 0 && xmax == 0) { xmin = 0; xmax = 5; }

    // Per-var nbins override (set by `--plot-config` or inline `-v` spec);
    // 0 = use the global --nbins.
    int nbins = vd.nbins > 0 ? vd.nbins : cliopts::nbins;

    std::string tag = channel + "_" + pm.fullname + "_" + vd.branch;
    SystHists sh;
    sh.h_cv = new TH1D(("h_cv_" + tag).c_str(), "", nbins, xmin, xmax);
    sh.h_cv->Sumw2();
    sh.tweakvals = pm.tweakvalues;
    for (int i = 0; i < pm.ntweaks; ++i) {
      TH1D *hv = new TH1D(Form("h_var%d_%s", i, tag.c_str()), "", nbins, xmin, xmax);
      hv->Sumw2();
      sh.h_var.push_back(hv);
    }
    var_hists[vd.branch] = sh;
  }
}

double GetVarValue(const EventVars &ev, const std::string &branch) {
  if (branch == "Enu_true") return ev.Enu;
  if (branch == "q0") return ev.q0;
  if (branch == "Q2") return ev.Q2;
  if (branch == "q3") return ev.q3;
  if (branch == "w") return ev.w;
  if (branch == "plep") return ev.plep;
  if (branch == "coslep") return ev.has_coslep ? ev.coslep : -999;
  if (branch == "EAvail_GeV") return ev.EAvail;
  if (branch == "leading_proton_KE") return ev.has_leading_proton ? ev.leading_proton_KE : -999;
  if (branch == "leading_proton_p")  return ev.leading_proton_p;
  if (branch == "Bjorken_x")         return ev.Bjorken_x;
  if (branch == "plep_L")            return ev.plep_L;
  if (branch == "plep_T")            return ev.plep_T;
  if (branch == "Emiss") return ev.Emiss;
  if (branch == "pmiss") return ev.pmiss;
  if (branch == "Ereco_cal") return ev.Ereco_cal;
  if (branch == "nproton") return (double)ev.nproton;
  if (branch == "npip") return (double)ev.npip;
  if (branch == "npim") return (double)ev.npim;
  if (branch == "npi0") return (double)ev.npi0;
  if (branch == "nneutron") return (double)ev.nneutron;
  return -999;
}

// ===== Histogram map type ==================================================
// allhists[channel][param][var] -> SystHists
using HistMap = std::map<std::string, std::map<std::string, std::map<std::string, SystHists>>>;
// Same shape, but the innermost key is a Var2DPair::key (e.g. "q0_q3") and
// each entry holds TH2Ds instead of TH1Ds.
using Hist2DMap = std::map<std::string, std::map<std::string, std::map<std::string, Syst2DHists>>>;

// ===== applies_to_channels loader (mirrors DumpConfiguredTweaksNuSyst) ======
// Re-parse the parameter-headers fhicl to extract per-provider
// `applies_to_channels` patterns. Returns an empty map if no provider declares
// the key, in which case no channel-aware skipping happens and behaviour
// matches the pre-filter baseline.
std::map<std::string, std::vector<std::string>>
LoadAppliesToChannelsMap(const std::string &fclname,
                        const std::string &top_key) {
  std::map<std::string, std::vector<std::string>> out;
  auto try_parse = [&](cet::filepath_maker &fm) -> fhicl::ParameterSet {
    return fhicl::ParameterSet::make(fclname, fm);
  };
  fhicl::ParameterSet raw;
  bool parsed = false;
  // 1) Treat the path as-is (handles absolute paths and paths relative to cwd).
  try {
    cet::filepath_maker fm;
    raw = try_parse(fm);
    parsed = true;
  } catch (std::exception const &) {}
  // 2) Fall back to FHICL_FILE_PATH lookup (handles bare basenames).
  if (!parsed) {
    try {
      cet::filepath_lookup_nonabsolute fm("FHICL_FILE_PATH");
      raw = try_parse(fm);
      parsed = true;
    } catch (std::exception const &e) {
      std::cerr << "[WARN]: Failed to load applies_to_channels map from "
                << fclname << ": " << e.what()
                << ". Falling back to no channel filtering." << std::endl;
      return {};
    }
  }
  try {
    fhicl::ParameterSet gen = raw.get<fhicl::ParameterSet>(top_key);
    auto provider_names = gen.get<std::vector<std::string>>("syst_providers");
    for (auto const &name : provider_names) {
      fhicl::ParameterSet prov =
          gen.get<fhicl::ParameterSet>(name, fhicl::ParameterSet{});
      fhicl::ParameterSet topts =
          prov.get<fhicl::ParameterSet>("tool_options", fhicl::ParameterSet{});
      auto patterns =
          topts.get<std::vector<std::string>>("applies_to_channels", {});
      if (!patterns.empty()) out[name] = patterns;
    }
  } catch (std::exception const &e) {
    std::cerr << "[WARN]: applies_to_channels lookup in " << fclname
              << " failed: " << e.what()
              << ". Falling back to no channel filtering." << std::endl;
    return {};
  }
  return out;
}

// Compute Bjorken x / pL / pT from base branches when the dedicated scalar
// branches aren't present in the input tree (older flat trees pre-dating
// nusyst tweaks shipping Bjorken_x/plep_L/plep_T). Mutates the relevant
// EventVars fields in place; leaves them untouched if the inputs needed
// for the fallback are also missing.
//
// Required inputs:
//   * Bjorken_x: just Q2 and q0 (always present in the tweaks tree).
//   * pL / pT  : isparticles_* and fsparticles_* vectors (per-particle
//                4-vectors). We pick the IS neutrino's direction (|PDG| in
//                {12, 14, 16}) and the leading-|p| FS lepton (|PDG| in
//                {11, 12, 13, 14, 15, 16}). If the vectors are absent, pL/pT
//                stay -999 and the consuming pair (pL_pT) silently skips.
inline void ApplyDerivedFallbacks(
    EventVars &evars,
    bool has_x, bool has_pL, bool has_pT,
    std::vector<int>    const *is_pdg, std::vector<double> const *is_px,
    std::vector<double> const *is_py,  std::vector<double> const *is_pz,
    std::vector<int>    const *fs_pdg, std::vector<double> const *fs_px,
    std::vector<double> const *fs_py,  std::vector<double> const *fs_pz) {

  // -- Bjorken x ---------------------------------------------------------
  if (!has_x) {
    static const double M_N = 0.93827203;  // GeV
    double safe_q0 = (evars.q0 > 1e-9) ? evars.q0 : 1e-9;
    evars.Bjorken_x = evars.Q2 / (2.0 * M_N * safe_q0);
  }

  // -- pL / pT -----------------------------------------------------------
  if (has_pL && has_pT) return;  // nothing to do

  // Preferred path: full per-particle vectors. Pick the IS neutrino's
  // direction and the leading-|p| FS lepton; project.
  bool from_vectors = false;
  if (is_pdg && fs_pdg && !is_pdg->empty() && !fs_pdg->empty()) {
    double nu_px = 0, nu_py = 0, nu_pz = 0, nu_mag2 = 0;
    for (size_t i = 0; i < is_pdg->size(); ++i) {
      int a = std::abs((*is_pdg)[i]);
      if (a == 12 || a == 14 || a == 16) {
        nu_px = (*is_px)[i]; nu_py = (*is_py)[i]; nu_pz = (*is_pz)[i];
        nu_mag2 = nu_px*nu_px + nu_py*nu_py + nu_pz*nu_pz;
        break;
      }
    }
    if (nu_mag2 > 1e-18) {
      double nu_mag = std::sqrt(nu_mag2);
      double nx = nu_px/nu_mag, ny = nu_py/nu_mag, nz = nu_pz/nu_mag;
      double best_p = -1, best_px = 0, best_py = 0, best_pz = 0;
      for (size_t i = 0; i < fs_pdg->size(); ++i) {
        int a = std::abs((*fs_pdg)[i]);
        if (a == 11 || a == 12 || a == 13 || a == 14 || a == 15 || a == 16) {
          double px = (*fs_px)[i], py = (*fs_py)[i], pz = (*fs_pz)[i];
          double pmag = std::sqrt(px*px + py*py + pz*pz);
          if (pmag > best_p) { best_p = pmag; best_px = px; best_py = py; best_pz = pz; }
        }
      }
      if (best_p >= 0) {
        double pL = best_px * nx + best_py * ny + best_pz * nz;
        double pT2 = std::max(0.0, best_p * best_p - pL * pL);
        if (!has_pL) evars.plep_L = pL;
        if (!has_pT) evars.plep_T = std::sqrt(pT2);
        from_vectors = true;
      }
    }
  }
  if (from_vectors) return;

  // Older-tree fallback: per-particle vectors are missing too, but plep,
  // Enu and q3 are always available. Solve for pL under the convention
  // that the incoming neutrino is along +z (true for typical GENIE flat
  // events): q3^2 = pT^2 + (Enu - pL)^2, plep^2 = pT^2 + pL^2  =>
  //   pL = (plep^2 + Enu^2 - q3^2) / (2 Enu)
  //   pT = sqrt(max(0, plep^2 - pL^2))
  // Reproduces the per-vector result when the neutrino is along z (the
  // typical case); a generic non-z beam direction would need the per-
  // particle vectors and bypasses this branch.
  if (evars.Enu > 1e-9 && evars.plep > 0) {
    double pL = (evars.plep * evars.plep + evars.Enu * evars.Enu - evars.q3 * evars.q3)
                / (2.0 * evars.Enu);
    double pT2 = std::max(0.0, evars.plep * evars.plep - pL * pL);
    if (!has_pL) evars.plep_L = pL;
    if (!has_pT) evars.plep_T = std::sqrt(pT2);
  }
}

// Resolve which provider declared a parameter by longest-prefix match on the
// param's full name (e.g. "GENIEReWeight_CCQE_MaCCQE" -> "GENIEReWeight_CCQE").
// Returns nullptr if no provider in the map is a prefix.
const std::vector<std::string> *
FindAppliesPatternsForFullname(
    const std::map<std::string, std::vector<std::string>> &map,
    const std::string &fullname) {
  const std::vector<std::string> *best = nullptr;
  size_t best_len = 0;
  for (auto const &kv : map) {
    auto const &provname = kv.first;
    if (fullname.size() > provname.size() + 1 &&
        fullname.compare(0, provname.size(), provname) == 0 &&
        fullname[provname.size()] == '_' &&
        provname.size() > best_len) {
      best = &kv.second;
      best_len = provname.size();
    }
  }
  return best;
}

// ===== FLAT TREE MODE ======================================================
void FillFromFlatTree(
    TTree *tree, TTree *meta,
    const std::vector<VarDef> &vars,
    const std::vector<ParamMeta> &params,
    const std::map<std::string, std::vector<std::string>> &applies_to_channels,
    HistMap &allhists,
    EventCounts &counts) {

  // Branch addresses for event classification
  double xsec = 0;
  int nu_pdg = 0, tgt_A = 0, tgt_Z = 0;
  bool is_cc = false, is_qe = false, is_mec = false, is_res = false, is_dis = false;
  double q0 = 0, Q2 = 0, q3 = 0, w = 0, plep = 0, Enu_true = 0, EAvail_GeV = 0;
  double Emiss = 0, pmiss = 0;
  std::vector<double> *fsprotons_KE = nullptr;
  std::vector<int> *fsi_pdgs = nullptr;

  // Build one TTreeFormula per formula-backed VarDef. nullptr for plain
  // built-in vars; we just check formulas[i] != nullptr at eval time.
  std::vector<std::unique_ptr<TTreeFormula>> formulas(vars.size());
  for (size_t i = 0; i < vars.size(); ++i) {
    if (vars[i].formula.empty()) continue;
    formulas[i] = std::make_unique<TTreeFormula>(
        ("v_" + vars[i].branch).c_str(), vars[i].formula.c_str(), tree);
    if (formulas[i]->GetNdim() == 0) {
      std::cerr << "[ERROR]: TTreeFormula for variable '" << vars[i].branch
                << "' failed to compile (expression: " << vars[i].formula
                << "). Aborting." << std::endl;
      std::exit(5);
    }
  }

  tree->SetBranchAddress("xsec", &xsec);
  tree->SetBranchAddress("nu_pdg", &nu_pdg);
  tree->SetBranchAddress("tgt_A", &tgt_A);
  tree->SetBranchAddress("tgt_Z", &tgt_Z);
  tree->SetBranchAddress("is_cc", &is_cc);
  tree->SetBranchAddress("is_qe", &is_qe);
  tree->SetBranchAddress("is_mec", &is_mec);
  tree->SetBranchAddress("is_res", &is_res);
  tree->SetBranchAddress("is_dis", &is_dis);
  tree->SetBranchAddress("q0", &q0);
  tree->SetBranchAddress("Q2", &Q2);
  tree->SetBranchAddress("q3", &q3);
  tree->SetBranchAddress("w", &w);
  tree->SetBranchAddress("plep", &plep);
  tree->SetBranchAddress("Enu_true", &Enu_true);
  tree->SetBranchAddress("EAvail_GeV", &EAvail_GeV);

  bool has_Emiss = tree->GetBranch("Emiss") != nullptr;
  bool has_pmiss = tree->GetBranch("pmiss") != nullptr;
  if (has_Emiss) tree->SetBranchAddress("Emiss", &Emiss);
  if (has_pmiss) tree->SetBranchAddress("pmiss", &pmiss);

  // Scalar branches added to the events tree by DumpConfiguredTweaksNuSyst
  // for compatibility with the kPredefinedVars list. has_<name> defaults to
  // false so older trees (pre-dating these branches) silently fall back to
  // the -999 sentinel and the per-event loop skips those vars.
  double Ereco_cal = 0;
  int nproton = 0, npip = 0, npim = 0, npi0 = 0, nneutron = 0;
  bool has_Ereco    = tree->GetBranch("Ereco_cal") != nullptr;
  bool has_nproton  = tree->GetBranch("nproton")   != nullptr;
  bool has_npip     = tree->GetBranch("npip")      != nullptr;
  bool has_npim     = tree->GetBranch("npim")      != nullptr;
  bool has_npi0     = tree->GetBranch("npi0")      != nullptr;
  bool has_nneutron = tree->GetBranch("nneutron")  != nullptr;
  if (has_Ereco)    tree->SetBranchAddress("Ereco_cal", &Ereco_cal);
  if (has_nproton)  tree->SetBranchAddress("nproton",   &nproton);
  if (has_npip)     tree->SetBranchAddress("npip",      &npip);
  if (has_npim)     tree->SetBranchAddress("npim",      &npim);
  if (has_npi0)     tree->SetBranchAddress("npi0",      &npi0);
  if (has_nneutron) tree->SetBranchAddress("nneutron",  &nneutron);

  bool has_proton_branch = tree->GetBranch("fsprotons_KE") != nullptr;
  if (has_proton_branch) tree->SetBranchAddress("fsprotons_KE", &fsprotons_KE);
  double leading_proton_p = -999;
  bool has_leading_proton_p = tree->GetBranch("leading_proton_p") != nullptr;
  if (has_leading_proton_p)
    tree->SetBranchAddress("leading_proton_p", &leading_proton_p);
  // x / pL / pT scalars (added to nusyst tweaks for the 2D pair defaults).
  double Bjorken_x = -999, plep_L = -999, plep_T = -999;
  bool has_x  = tree->GetBranch("Bjorken_x") != nullptr;
  bool has_pL = tree->GetBranch("plep_L")    != nullptr;
  bool has_pT = tree->GetBranch("plep_T")    != nullptr;
  if (has_x)  tree->SetBranchAddress("Bjorken_x", &Bjorken_x);
  if (has_pL) tree->SetBranchAddress("plep_L",    &plep_L);
  if (has_pT) tree->SetBranchAddress("plep_T",    &plep_T);
  // Per-particle 4-vectors for fallback computation of x/pL/pT when the
  // dedicated scalar branches are absent (older flat trees).
  std::vector<int>    *is_pdg = nullptr, *fs_pdg = nullptr;
  std::vector<double> *is_px = nullptr, *is_py = nullptr, *is_pz = nullptr;
  std::vector<double> *fs_px = nullptr, *fs_py = nullptr, *fs_pz = nullptr;
  if (tree->GetBranch("isparticles_pdg")) tree->SetBranchAddress("isparticles_pdg", &is_pdg);
  if (tree->GetBranch("isparticles_px"))  tree->SetBranchAddress("isparticles_px",  &is_px);
  if (tree->GetBranch("isparticles_py"))  tree->SetBranchAddress("isparticles_py",  &is_py);
  if (tree->GetBranch("isparticles_pz"))  tree->SetBranchAddress("isparticles_pz",  &is_pz);
  if (tree->GetBranch("fsparticles_pdg")) tree->SetBranchAddress("fsparticles_pdg", &fs_pdg);
  if (tree->GetBranch("fsparticles_px"))  tree->SetBranchAddress("fsparticles_px",  &fs_px);
  if (tree->GetBranch("fsparticles_py"))  tree->SetBranchAddress("fsparticles_py",  &fs_py);
  if (tree->GetBranch("fsparticles_pz"))  tree->SetBranchAddress("fsparticles_pz",  &fs_pz);

  bool has_fsi_pdgs = tree->GetBranch("fsi_pdgs") != nullptr;
  if (has_fsi_pdgs) tree->SetBranchAddress("fsi_pdgs", &fsi_pdgs);

  // Tweak response branches
  struct PBranch {
    std::string fullname;
    int ntweaks;
    std::vector<double> responses;
    double cv_weight;
    std::vector<std::string> applies_patterns; // empty -> applies to all channels
  };
  std::vector<PBranch> pbranches;
  size_t n_filtered = 0;
  for (auto &pm : params) {
    PBranch pb;
    pb.fullname = pm.fullname;
    pb.ntweaks = pm.ntweaks;
    pb.responses.resize(pm.ntweaks, 1.0);
    pb.cv_weight = 1.0;
    if (tree->GetBranch(("tweak_responses_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("tweak_responses_" + pm.fullname).c_str(), pb.responses.data());
    if (tree->GetBranch(("paramCVWeight_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("paramCVWeight_" + pm.fullname).c_str(), &pb.cv_weight);
    if (auto *pats = FindAppliesPatternsForFullname(applies_to_channels, pm.fullname)) {
      pb.applies_patterns = *pats;
      ++n_filtered;
    }
    pbranches.push_back(std::move(pb));
  }
  if (!applies_to_channels.empty()) {
    std::cout << "[INFO]: " << n_filtered << " of " << pbranches.size()
              << " parameters have applies_to_channels patterns; channel x dial "
                 "combinations outside those patterns will be skipped."
              << std::endl;
  }

  Long64_t nentries = tree->GetEntries(); if (cliopts::NMax < (size_t)nentries) nentries = (Long64_t)cliopts::NMax;
  Long64_t nshout = std::max(nentries / 20, (Long64_t)1);

  for (Long64_t ev = 0; ev < nentries; ++ev) {
    tree->GetEntry(ev);
    if (!(ev % nshout))
      std::cout << "\rProcessing event " << ev << "/" << nentries << std::flush;

    EventVars evars;
    evars.Enu = Enu_true; evars.q0 = q0; evars.Q2 = Q2; evars.q3 = q3;
    evars.w = w; evars.plep = plep; evars.EAvail = EAvail_GeV; evars.xsec = xsec;
    evars.is_cc = is_cc; evars.is_qe = is_qe; evars.is_mec = is_mec;
    evars.is_res = is_res; evars.is_dis = is_dis;
    evars.nu_pdg = nu_pdg; evars.tgt_A = tgt_A; evars.tgt_Z = tgt_Z;
    evars.has_coslep = false; evars.coslep = -999;
    evars.has_leading_proton = has_proton_branch && fsprotons_KE && !fsprotons_KE->empty();
    evars.leading_proton_KE = evars.has_leading_proton ? (*fsprotons_KE)[0] : -999;
    evars.Emiss = has_Emiss ? Emiss : -999;
    evars.pmiss = has_pmiss ? pmiss : -999;

    // Calorimetric reco energy + multiplicities: prefer the new dedicated
    // scalar branches written by DumpConfiguredTweaksNuSyst; if the tree
    // pre-dates them (has_<name> false), fall back to legacy derivations.
    if (has_Ereco) {
      evars.Ereco_cal = Ereco_cal;
    } else {
      double Elep = std::sqrt(plep * plep + 0.10566 * 0.10566); // assume muon
      evars.Ereco_cal = EAvail_GeV + Elep;
    }
    evars.leading_proton_p = has_leading_proton_p ? leading_proton_p : -999;
    evars.Bjorken_x        = has_x ? Bjorken_x : -999;
    evars.plep_L           = has_pL ? plep_L   : -999;
    evars.plep_T           = has_pT ? plep_T   : -999;
    ApplyDerivedFallbacks(evars, has_x, has_pL, has_pT,
                          is_pdg, is_px, is_py, is_pz,
                          fs_pdg, fs_px, fs_py, fs_pz);
    evars.nproton  = has_nproton  ? nproton
                                  : (has_proton_branch && fsprotons_KE
                                        ? (int)fsprotons_KE->size() : 0);
    evars.npip     = has_npip     ? npip     : 0;
    evars.npim     = has_npim     ? npim     : 0;
    evars.npi0     = has_npi0     ? npi0     : 0;
    evars.nneutron = has_nneutron ? nneutron : 0;

    counts.Count(evars);
    std::string chkey = MakeChannelKey(evars);

    for (size_t ip = 0; ip < pbranches.size(); ++ip) {
      auto &pb = pbranches[ip];
      auto &pm = params[ip];

      // Skip booking + filling for (param, channel) combinations the dial
      // does not apply to. Empty patterns -> dial applies everywhere.
      if (!pb.applies_patterns.empty() &&
          !nusyst::channel::MatchesAny(chkey, pb.applies_patterns)) {
        continue;
      }

      // Book histograms on first encounter of this channel
      if (allhists[chkey].find(pm.fullname) == allhists[chkey].end()) {
        BookHistograms(chkey, pm, vars, allhists[chkey][pm.fullname], tree);
      }

      auto &var_hists = allhists[chkey][pm.fullname];
      for (size_t vi = 0; vi < vars.size(); ++vi) {
        auto const &vd = vars[vi];
        double x;
        if (formulas[vi]) {
          // TTreeFormula::EvalInstance(0) returns the value for the current
          // tree entry; for Max$/Sum$ aggregate forms it folds across the
          // referenced vectors and returns a single scalar per event.
          x = formulas[vi]->EvalInstance(0);
        } else {
          x = GetVarValue(evars, vd.branch);
          if (x <= -999) continue;
        }
        auto vit = var_hists.find(vd.branch);
        if (vit == var_hists.end()) continue;
        auto &sh = vit->second;

        sh.h_cv->Fill(x, xsec * pb.cv_weight);
        for (int i = 0; i < pb.ntweaks && i < (int)sh.h_var.size(); ++i)
          sh.h_var[i]->Fill(x, xsec * pb.cv_weight * pb.responses[i]);
      }
    }
  }
  std::cout << "\rProcessed " << nentries << " events.          " << std::endl;
  tree->ResetBranchAddresses();
}

// ===== FLAT TREE MODE -- 2D ==================================================
// Mirrors FillFromFlatTree but fills TH2Ds per (channel, dial, pair) instead
// of TH1Ds per (channel, dial, var). 1D event reading + channel/applies
// logic is re-used by calling GetVarValue() per axis variable.
//
// Histograms are booked lazily on first event for a given (channel, dial)
// using the binning carried by each pair's x and y VarDef entries in the
// runtime registry.
void FillFromFlatTree2D(
    TTree *tree, TTree *meta,
    const std::vector<Var2DPair> &pairs,
    const std::vector<ParamMeta> &params,
    const std::map<std::string, std::vector<std::string>> &applies_to_channels,
    Hist2DMap &allhists2d,
    EventCounts &counts) {

  // Resolve VarDefs for each pair up front; abort early if any axis var
  // isn't in the registry.
  struct PairResolved {
    Var2DPair pair;
    VarDef    x_vd;
    VarDef    y_vd;
  };
  std::vector<PairResolved> resolved;
  for (auto const &p : pairs) {
    auto it_x = g_vars_registry.find(p.x_var);
    auto it_y = g_vars_registry.find(p.y_var);
    if (it_x == g_vars_registry.end() || it_y == g_vars_registry.end()) {
      std::cerr << "[ERROR]: 2D pair '" << p.key << "' references unknown var ("
                << p.x_var << " or " << p.y_var << "). Skipping." << std::endl;
      continue;
    }
    resolved.push_back({p, it_x->second, it_y->second});
  }

  // Same branch reading as FillFromFlatTree.
  double xsec = 0;
  int nu_pdg = 0, tgt_A = 0, tgt_Z = 0;
  bool is_cc = false, is_qe = false, is_mec = false, is_res = false, is_dis = false;
  double q0 = 0, Q2 = 0, q3 = 0, w = 0, plep = 0, Enu_true = 0, EAvail_GeV = 0;
  double Emiss = 0, pmiss = 0;
  std::vector<double> *fsprotons_KE = nullptr;
  std::vector<int> *fsi_pdgs = nullptr;
  tree->SetBranchAddress("xsec",  &xsec);
  tree->SetBranchAddress("nu_pdg", &nu_pdg);
  tree->SetBranchAddress("tgt_A",  &tgt_A);
  tree->SetBranchAddress("tgt_Z",  &tgt_Z);
  tree->SetBranchAddress("is_cc",  &is_cc);
  tree->SetBranchAddress("is_qe",  &is_qe);
  tree->SetBranchAddress("is_mec", &is_mec);
  tree->SetBranchAddress("is_res", &is_res);
  tree->SetBranchAddress("is_dis", &is_dis);
  tree->SetBranchAddress("q0", &q0);
  tree->SetBranchAddress("Q2", &Q2);
  tree->SetBranchAddress("q3", &q3);
  tree->SetBranchAddress("w",  &w);
  tree->SetBranchAddress("plep",       &plep);
  tree->SetBranchAddress("Enu_true",   &Enu_true);
  tree->SetBranchAddress("EAvail_GeV", &EAvail_GeV);
  bool has_Emiss = tree->GetBranch("Emiss") != nullptr;
  bool has_pmiss = tree->GetBranch("pmiss") != nullptr;
  if (has_Emiss) tree->SetBranchAddress("Emiss", &Emiss);
  if (has_pmiss) tree->SetBranchAddress("pmiss", &pmiss);
  double Ereco_cal = 0;
  int nproton = 0, npip = 0, npim = 0, npi0 = 0, nneutron = 0;
  bool has_Ereco    = tree->GetBranch("Ereco_cal") != nullptr;
  bool has_nproton  = tree->GetBranch("nproton")   != nullptr;
  bool has_npip     = tree->GetBranch("npip")      != nullptr;
  bool has_npim     = tree->GetBranch("npim")      != nullptr;
  bool has_npi0     = tree->GetBranch("npi0")      != nullptr;
  bool has_nneutron = tree->GetBranch("nneutron")  != nullptr;
  if (has_Ereco)    tree->SetBranchAddress("Ereco_cal", &Ereco_cal);
  if (has_nproton)  tree->SetBranchAddress("nproton",   &nproton);
  if (has_npip)     tree->SetBranchAddress("npip",      &npip);
  if (has_npim)     tree->SetBranchAddress("npim",      &npim);
  if (has_npi0)     tree->SetBranchAddress("npi0",      &npi0);
  if (has_nneutron) tree->SetBranchAddress("nneutron",  &nneutron);
  bool has_proton_branch = tree->GetBranch("fsprotons_KE") != nullptr;
  if (has_proton_branch) tree->SetBranchAddress("fsprotons_KE", &fsprotons_KE);
  double leading_proton_p = -999;
  bool has_leading_proton_p = tree->GetBranch("leading_proton_p") != nullptr;
  if (has_leading_proton_p)
    tree->SetBranchAddress("leading_proton_p", &leading_proton_p);
  double Bjorken_x = -999, plep_L = -999, plep_T = -999;
  bool has_x  = tree->GetBranch("Bjorken_x") != nullptr;
  bool has_pL = tree->GetBranch("plep_L")    != nullptr;
  bool has_pT = tree->GetBranch("plep_T")    != nullptr;
  if (has_x)  tree->SetBranchAddress("Bjorken_x", &Bjorken_x);
  if (has_pL) tree->SetBranchAddress("plep_L",    &plep_L);
  if (has_pT) tree->SetBranchAddress("plep_T",    &plep_T);
  // Per-particle 4-vectors for fallback computation of x/pL/pT when the
  // dedicated scalar branches are absent (older flat trees).
  std::vector<int>    *is_pdg = nullptr, *fs_pdg = nullptr;
  std::vector<double> *is_px = nullptr, *is_py = nullptr, *is_pz = nullptr;
  std::vector<double> *fs_px = nullptr, *fs_py = nullptr, *fs_pz = nullptr;
  if (tree->GetBranch("isparticles_pdg")) tree->SetBranchAddress("isparticles_pdg", &is_pdg);
  if (tree->GetBranch("isparticles_px"))  tree->SetBranchAddress("isparticles_px",  &is_px);
  if (tree->GetBranch("isparticles_py"))  tree->SetBranchAddress("isparticles_py",  &is_py);
  if (tree->GetBranch("isparticles_pz"))  tree->SetBranchAddress("isparticles_pz",  &is_pz);
  if (tree->GetBranch("fsparticles_pdg")) tree->SetBranchAddress("fsparticles_pdg", &fs_pdg);
  if (tree->GetBranch("fsparticles_px"))  tree->SetBranchAddress("fsparticles_px",  &fs_px);
  if (tree->GetBranch("fsparticles_py"))  tree->SetBranchAddress("fsparticles_py",  &fs_py);
  if (tree->GetBranch("fsparticles_pz"))  tree->SetBranchAddress("fsparticles_pz",  &fs_pz);
  bool has_fsi_pdgs = tree->GetBranch("fsi_pdgs") != nullptr;
  if (has_fsi_pdgs) tree->SetBranchAddress("fsi_pdgs", &fsi_pdgs);

  // Per-dial tweak response branch setup (same as 1D path).
  struct PBranch {
    std::string fullname;
    int ntweaks;
    std::vector<double> responses;
    double cv_weight;
    std::vector<std::string> applies_patterns;
  };
  std::vector<PBranch> pbranches;
  for (auto &pm : params) {
    PBranch pb;
    pb.fullname = pm.fullname;
    pb.ntweaks = pm.ntweaks;
    pb.responses.resize(pm.ntweaks, 1.0);
    pb.cv_weight = 1.0;
    if (tree->GetBranch(("tweak_responses_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("tweak_responses_" + pm.fullname).c_str(),
                              pb.responses.data());
    if (tree->GetBranch(("paramCVWeight_" + pm.fullname).c_str()))
      tree->SetBranchAddress(("paramCVWeight_" + pm.fullname).c_str(),
                              &pb.cv_weight);
    if (auto *pats = FindAppliesPatternsForFullname(applies_to_channels, pm.fullname))
      pb.applies_patterns = *pats;
    pbranches.push_back(std::move(pb));
  }

  Long64_t nentries = tree->GetEntries();
  if (cliopts::NMax < (size_t)nentries) nentries = (Long64_t)cliopts::NMax;
  Long64_t nshout = std::max(nentries / 20, (Long64_t)1);

  for (Long64_t ev = 0; ev < nentries; ++ev) {
    if (ev % nshout == 0)
      std::cout << "\rProcessing event " << ev << "/" << nentries << std::flush;
    tree->GetEntry(ev);

    EventVars evars;
    evars.Enu = Enu_true; evars.q0 = q0; evars.Q2 = Q2; evars.q3 = q3;
    evars.w = w; evars.plep = plep; evars.EAvail = EAvail_GeV; evars.xsec = xsec;
    evars.is_cc = is_cc; evars.is_qe = is_qe; evars.is_mec = is_mec;
    evars.is_res = is_res; evars.is_dis = is_dis;
    evars.nu_pdg = nu_pdg; evars.tgt_A = tgt_A; evars.tgt_Z = tgt_Z;
    evars.has_coslep = false; evars.coslep = -999;
    evars.has_leading_proton = has_proton_branch && fsprotons_KE && !fsprotons_KE->empty();
    evars.leading_proton_KE = evars.has_leading_proton ? (*fsprotons_KE)[0] : -999;
    evars.Emiss = has_Emiss ? Emiss : -999;
    evars.pmiss = has_pmiss ? pmiss : -999;
    evars.Ereco_cal = has_Ereco ? Ereco_cal : -999;
    evars.leading_proton_p = has_leading_proton_p ? leading_proton_p : -999;
    evars.Bjorken_x = has_x  ? Bjorken_x : -999;
    evars.plep_L    = has_pL ? plep_L    : -999;
    evars.plep_T    = has_pT ? plep_T    : -999;
    ApplyDerivedFallbacks(evars, has_x, has_pL, has_pT,
                          is_pdg, is_px, is_py, is_pz,
                          fs_pdg, fs_px, fs_py, fs_pz);
    evars.nproton  = has_nproton  ? nproton  : 0;
    evars.npip     = has_npip     ? npip     : 0;
    evars.npim     = has_npim     ? npim     : 0;
    evars.npi0     = has_npi0     ? npi0     : 0;
    evars.nneutron = has_nneutron ? nneutron : 0;

    counts.Count(evars);
    std::string chkey = MakeChannelKey(evars);

    for (size_t ip = 0; ip < pbranches.size(); ++ip) {
      auto &pb = pbranches[ip];
      auto &pm = params[ip];
      if (!pb.applies_patterns.empty() &&
          !nusyst::channel::MatchesAny(chkey, pb.applies_patterns))
        continue;

      // Lazy-book TH2Ds for this (channel, dial) on first encounter.
      if (allhists2d[chkey].find(pm.fullname) == allhists2d[chkey].end()) {
        auto &pair_map = allhists2d[chkey][pm.fullname];
        for (auto const &pr : resolved) {
          std::string tag = chkey + "_" + pm.fullname + "_" + pr.pair.key;
          int xn = pr.x_vd.nbins > 0 ? pr.x_vd.nbins : cliopts::nbins;
          int yn = pr.y_vd.nbins > 0 ? pr.y_vd.nbins : cliopts::nbins;
          double xlo = pr.x_vd.xmin, xhi = pr.x_vd.xmax;
          double ylo = pr.y_vd.xmin, yhi = pr.y_vd.xmax;
          if (xlo == 0 && xhi == 0) { xlo = 0; xhi = 5; }
          if (ylo == 0 && yhi == 0) { ylo = 0; yhi = 5; }
          Syst2DHists sh;
          sh.tweakvals = pm.tweakvalues;
          sh.h_cv = new TH2D(("h2_cv_" + tag).c_str(), "",
                              xn, xlo, xhi, yn, ylo, yhi);
          sh.h_cv->Sumw2();
          for (int i = 0; i < pm.ntweaks; ++i) {
            TH2D *hv = new TH2D(Form("h2_var%d_%s", i, tag.c_str()), "",
                                xn, xlo, xhi, yn, ylo, yhi);
            hv->Sumw2();
            sh.h_var.push_back(hv);
          }
          pair_map[pr.pair.key] = sh;
        }
      }

      auto &pair_map = allhists2d[chkey][pm.fullname];
      for (auto const &pr : resolved) {
        double x = GetVarValue(evars, pr.x_vd.branch);
        double y = GetVarValue(evars, pr.y_vd.branch);
        if (x <= -999 || y <= -999) continue;
        auto pit = pair_map.find(pr.pair.key);
        if (pit == pair_map.end()) continue;
        auto &sh = pit->second;
        sh.h_cv->Fill(x, y, xsec * pb.cv_weight);
        for (int i = 0; i < pb.ntweaks && i < (int)sh.h_var.size(); ++i)
          sh.h_var[i]->Fill(x, y, xsec * pb.cv_weight * pb.responses[i]);
      }
    }
  }
  std::cout << "\rProcessed " << nentries << " events (2D).          " << std::endl;
  tree->ResetBranchAddresses();
}

// ===== GHEP MODE ===========================================================
void FillFromGHEP(
    const std::string &input, const std::string &fclname,
    const std::vector<VarDef> &vars,
    const std::map<std::string, std::vector<std::string>> &applies_to_channels,
    HistMap &allhists,
    EventCounts &counts) {

  // Two-step load so a -p filter can drop providers whose dials all miss
  // before any provider is instantiated. Mirrors nusyst tweaks / response.
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
      raw_ps = fhicl::ParameterSet::make(fclname, *fm);
    }
    fhicl::ParameterSet gen_ps =
        raw_ps.get<fhicl::ParameterSet>(cliopts::fhicl_key);

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
            // Reuse the existing ParamSelected() helper (matches the
            // provider-qualified `<provider>_<dial>` form as well as the
            // bare dial name, so the user can pass either).
            if (ParamSelected(pretty) ||
                ParamSelected(pname + "_" + pretty)) {
              any_match = true; break;
            }
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

  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(input.c_str())) {
    std::cout << "[ERROR]: No gtree found in " << input << std::endl;
    return;
  }
  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  gevs->SetBranchAddress(cliopts::genie_branch_name.c_str(), &GenieNtpl);

  // Build param metadata
  struct GPInfo {
    paramId_t pid;
    ParamMeta meta;
    std::vector<std::string> applies_patterns; // empty -> applies to all channels
  };
  std::vector<GPInfo> gparams;
  size_t n_filtered = 0;
  for (paramId_t pid : phh.GetParameters()) {
    SystParamHeader const &hdr = phh.GetHeader(pid);
    if (hdr.isResponselessParam) continue;
    std::string provname;
    for (size_t si = 0; si < phh.GetSystProvider().size(); ++si)
      if (phh.GetSystProvider()[si]->ParamIsHandled(pid))
        { provname = phh.GetSystProvider()[si]->GetFullyQualifiedName(); break; }
    std::string fullname = provname + "_" + hdr.prettyName;
    if (!ParamSelected(fullname)) continue;
    GPInfo gi;
    gi.pid = pid;
    gi.meta.fullname = fullname;
    gi.meta.ntweaks = hdr.isCorrection ? 1 : (int)hdr.paramVariations.size();
    gi.meta.tweakvalues = hdr.paramVariations;
    auto it = applies_to_channels.find(provname);
    if (it != applies_to_channels.end()) {
      gi.applies_patterns = it->second;
      ++n_filtered;
    }
    gparams.push_back(gi);
  }
  if (gparams.empty()) { std::cout << "[WARN]: No matching parameters." << std::endl; delete gevs; return; }
  if (!applies_to_channels.empty()) {
    std::cout << "[INFO]: " << n_filtered << " of " << gparams.size()
              << " parameters have applies_to_channels patterns; channel x dial "
                 "combinations outside those patterns will be skipped."
              << std::endl;
  }

  Long64_t nentries = gevs->GetEntries(); if (cliopts::NMax < (size_t)nentries) nentries = (Long64_t)cliopts::NMax;
  Long64_t nshout = std::max(nentries / 20, (Long64_t)1);

  for (Long64_t ev_it = 0; ev_it < nentries; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;
    if (!(ev_it % nshout))
      std::cout << "\rGHEP event " << ev_it << "/" << nentries << std::flush;

    genie::GHepParticle *FSLep = GenieGHep.FinalStatePrimaryLepton();
    genie::GHepParticle *ISLep = GenieGHep.Probe();
    if (!FSLep || !ISLep) { GenieNtpl->Clear(); continue; }

    TLorentzVector FSLepP4 = *FSLep->P4();
    TLorentzVector ISLepP4 = *ISLep->P4();
    TLorentzVector emTransfer = ISLepP4 - FSLepP4;

    EventVars evars;
    evars.xsec = GenieGHep.XSec() / genie::units::cm2;
    evars.Enu = ISLepP4.E();
    evars.q0 = emTransfer.E();
    evars.Q2 = -emTransfer.Mag2();
    evars.q3 = emTransfer.Vect().Mag();
    double MN = 0.93827203;
    double W2 = MN*MN + 2.0*MN*evars.q0 - evars.Q2;
    evars.w = (W2 > 0) ? std::sqrt(W2) : -1;
    evars.plep = FSLepP4.Vect().Mag();
    evars.has_coslep = true;
    evars.coslep = FSLepP4.CosTheta();
    evars.EAvail = GetErecoil_MINERvA_LowRecoil(GenieGHep);

    // Bjorken x with q0 -> 0 guard.
    {
      double safe_q0 = (evars.q0 > 1e-9) ? evars.q0 : 1e-9;
      evars.Bjorken_x = evars.Q2 / (2.0 * MN * safe_q0);
    }
    // pL / pT relative to the IS-neutrino direction.
    {
      TVector3 nuDir = ISLepP4.Vect();
      double nu_mag = nuDir.Mag();
      if (nu_mag > 1e-9) {
        nuDir *= 1.0 / nu_mag;
        TVector3 lp = FSLepP4.Vect();
        double pL = lp.Dot(nuDir);
        double pmag = lp.Mag();
        double pT2 = std::max(0.0, pmag * pmag - pL * pL);
        evars.plep_L = pL;
        evars.plep_T = std::sqrt(pT2);
      } else {
        evars.plep_L = evars.plep;
        evars.plep_T = 0.0;
      }
    }

    // FS particle loop: proton KE & |p|, multiplicities
    evars.has_leading_proton = false;
    evars.leading_proton_KE = -999;
    evars.leading_proton_p  = -999;
    evars.nproton = 0; evars.npip = 0; evars.npim = 0; evars.npi0 = 0; evars.nneutron = 0;
    genie::GHepParticle *p = nullptr;
    TIter iter(&GenieGHep);
    double max_pKE = -1;
    double max_pp  = -1;
    while ((p = dynamic_cast<genie::GHepParticle*>(iter.Next()))) {
      if (p->Status() != genie::kIStStableFinalState) continue;
      int pdgc = p->Pdg();
      if (pdgc == 2212) {
        evars.nproton++;
        double ke = p->KinE();
        if (ke > max_pKE) { max_pKE = ke; evars.has_leading_proton = true; evars.leading_proton_KE = ke; }
        double pmag = p->P4()->Vect().Mag();
        if (pmag > max_pp) { max_pp = pmag; evars.leading_proton_p = pmag; }
      }
      else if (pdgc ==  211) evars.npip++;
      else if (pdgc == -211) evars.npim++;
      else if (pdgc ==  111) evars.npi0++;
      else if (pdgc == 2112) evars.nneutron++;
    }

    evars.Emiss = GetEmiss(GenieGHep, false);
    evars.pmiss = GetPmiss(GenieGHep, false);
    double Elep = std::sqrt(evars.plep * evars.plep + 0.10566 * 0.10566);
    evars.Ereco_cal = evars.EAvail + Elep;

    evars.is_cc = GenieGHep.Summary()->ProcInfo().IsWeakCC();
    evars.is_qe = GenieGHep.Summary()->ProcInfo().IsQuasiElastic();
    evars.is_mec = GenieGHep.Summary()->ProcInfo().IsMEC();
    evars.is_res = GenieGHep.Summary()->ProcInfo().IsResonant();
    evars.is_dis = GenieGHep.Summary()->ProcInfo().IsDeepInelastic();
    evars.nu_pdg = ISLep->Pdg();
    evars.tgt_A = GenieGHep.TargetNucleus() ? GenieGHep.TargetNucleus()->A() : 1;
    evars.tgt_Z = GenieGHep.TargetNucleus() ? GenieGHep.TargetNucleus()->Z() : 1;

    counts.Count(evars);
    std::string chkey = MakeChannelKey(evars);

    event_unit_response_w_cv_t resp = phh.GetEventVariationAndCVResponse(GenieGHep);

    for (auto &gp : gparams) {
      // Skip booking + filling for (param, channel) combinations the dial
      // does not apply to. Empty patterns -> dial applies everywhere.
      if (!gp.applies_patterns.empty() &&
          !nusyst::channel::MatchesAny(chkey, gp.applies_patterns)) {
        continue;
      }
      if (allhists[chkey].find(gp.meta.fullname) == allhists[chkey].end())
        BookHistograms(chkey, gp.meta, vars, allhists[chkey][gp.meta.fullname]);

      size_t ridx = GetParamContainerIndex(resp, gp.pid);
      double cv_w = 1.0;
      std::vector<double> responses;
      if (ridx != kParamUnhandled<size_t>) {
        cv_w = resp[ridx].CV_response;
        responses = resp[ridx].responses;
      }

      auto &var_hists = allhists[chkey][gp.meta.fullname];
      for (auto &vd : vars) {
        // Formula vars require a TTree behind them; GHEP mode reads
        // genie::EventRecord directly and has no tree to bind to. Skip with
        // a one-time warning (issued from main before the loop starts).
        if (!vd.formula.empty()) continue;
        double x = GetVarValue(evars, vd.branch);
        if (x <= -999) continue;
        auto vit = var_hists.find(vd.branch);
        if (vit == var_hists.end()) continue;
        auto &sh = vit->second;
        sh.h_cv->Fill(x, evars.xsec * cv_w);
        for (int i = 0; i < (int)responses.size() && i < (int)sh.h_var.size(); ++i)
          sh.h_var[i]->Fill(x, evars.xsec * cv_w * responses[i]);
      }
    }
    GenieNtpl->Clear();
  }
  std::cout << "\rProcessed " << nentries << " GHEP events.          " << std::endl;
  delete gevs;
}

// ===== Scale histograms to differential xsec ===============================
// Per-event xsec values were summed into each bin during the fill pass.
// To turn that into the canonical differential cross section per nucleus,
// divide by:
//   * the total number of events processed (so the average per event is the
//     true sigma_tot for the channel) and
//   * the bin width (so dsigma/dx is a density).
// The factor sums to:  bin = (1/N_events) * Sum(xsec * w) / dx,  units
// cm^2 / (x-axis unit) per nucleus -- which is what BuildYAxisTitle reports.
void ScaleToDiffXsec(HistMap &allhists, size_t n_events_processed) {
  double inv_n = (n_events_processed > 0) ? 1.0 / double(n_events_processed) : 1.0;
  for (auto &[ch, param_map] : allhists)
    for (auto &[par, var_map] : param_map)
      for (auto &[var, sh] : var_map) {
        sh.h_cv->Scale(inv_n, "width");
        for (auto *h : sh.h_var) h->Scale(inv_n, "width");
      }
}

// ===== Remove empty channels ===============================================
void PruneEmpty(HistMap &allhists) {
  std::vector<std::string> empty_channels;
  for (auto &[ch, param_map] : allhists) {
    bool has_entries = false;
    for (auto &[par, var_map] : param_map)
      for (auto &[var, sh] : var_map)
        if (sh.h_cv->GetEntries() > 0) has_entries = true;
    if (!has_entries) empty_channels.push_back(ch);
  }
  for (auto &ch : empty_channels) allhists.erase(ch);
}

// Drop (channel, param) pairs whose variation histograms are bin-for-bin
// identical to the CV histogram across every variable. These arise when a
// dial's applies_to_channels topology matches but its internal GENIE cuts
// (struck-nucleon PDG, FS pion multiplicity, W cuts, FS-hadron presence for
// FSI dials, etc.) leave every event in the channel with response = 1.
void PruneTrivialParams(HistMap &allhists) {
  size_t pruned = 0, total = 0;
  for (auto &chan_kv : allhists) {
    auto &param_map = chan_kv.second;
    for (auto it = param_map.begin(); it != param_map.end(); ) {
      ++total;
      bool trivial = true;
      for (auto const &var_kv : it->second) {
        auto const &sh = var_kv.second;
        int nb = sh.h_cv->GetNbinsX();
        for (auto *hv : sh.h_var) {
          for (int b = 1; b <= nb && trivial; ++b) {
            double cv = sh.h_cv->GetBinContent(b);
            double v  = hv->GetBinContent(b);
            // Bit-exact compare is enough: when a dial doesn't fire, each
            // event's CV-weighted and var-weighted fills are identical
            // (response is exactly 1.0), so the bin sums match exactly.
            if (cv != v) trivial = false;
          }
          if (!trivial) break;
        }
        if (!trivial) break;
      }
      if (trivial) {
        for (auto &var_kv : it->second) {
          auto &sh = var_kv.second;
          delete sh.h_cv;
          for (auto *h : sh.h_var) delete h;
        }
        it = param_map.erase(it);
        ++pruned;
      } else {
        ++it;
      }
    }
  }
  std::cout << "[INFO]: Pruned " << pruned << " of " << total
            << " (channel, dial) combinations with response == CV in every bin."
            << std::endl;
}

// ===== Color palette =======================================================
// Diverging blue (negative variations) → grey (CV / 0) → red (positive),
// chosen to be readable for protanopia / deuteranopia (blue vs red is
// preserved under both). Shades deepen with |variation|; values beyond
// ±3 saturate at the extreme shade.
void GetVarColor(double var_value, int &color, int &style) {
  style = 1;  // solid line -- color is the only differentiator
  if (std::abs(var_value) < 1e-9) {
    color = TColor::GetColor("#404040");  // dark grey for the 0-σ sample
    return;
  }
  // Three-class blue / red shades, lightest at |1σ|, darkest at |3σ|.
  // The hex values are from a ColorBrewer RdBu 7-class diverging palette.
  static const int blue_shades[] = {
    TColor::GetColor("#4393C3"),   // |var| ≈ 1
    TColor::GetColor("#2166AC"),   // |var| ≈ 2
    TColor::GetColor("#053061"),   // |var| ≥ 3
  };
  static const int red_shades[] = {
    TColor::GetColor("#F4A582"),   // |var| ≈ 1
    TColor::GetColor("#D6604D"),   // |var| ≈ 2
    TColor::GetColor("#B2182B"),   // |var| ≥ 3
  };
  double a = std::abs(var_value);
  int idx = (a >= 3.0) ? 2 : (a >= 2.0) ? 1 : 0;
  color = (var_value < 0) ? blue_shades[idx] : red_shades[idx];
}

// ===== Summary page ========================================================
void MakeSummaryPage(TCanvas *c, const EventCounts &counts) {
  c->Clear();
  c->cd();
  TPad *pad = new TPad("summary", "", 0, 0, 1, 1);
  pad->SetFillColor(kWhite);
  pad->Draw();
  pad->cd();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);

  double y = 0.92;
  double dy = 0.035;

  tex.SetTextSize(0.05);
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Event Summary"); y -= 1.5 * dy;
  tex.SetTextFont(42);
  tex.SetTextSize(0.032);

  tex.DrawLatex(0.05, y, Form("Total events processed: %zu", counts.total)); y -= 1.2 * dy;

  // Neutrino species
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Neutrino species:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[pdg, n] : counts.by_nu_pdg) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", NuName(pdg).c_str(), n));
    y -= dy;
  }

  // Targets
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "Targets:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[A, n] : counts.by_tgt_A) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", TargetName(A).c_str(), n));
    y -= dy;
  }

  // By current
  y -= 0.5 * dy;
  tex.SetTextFont(62);
  tex.DrawLatex(0.05, y, "By current:"); y -= dy;
  tex.SetTextFont(42);
  for (auto &[cur, n] : counts.by_current) {
    tex.DrawLatex(0.08, y, Form("%-12s %zu", cur.c_str(), n));
    y -= dy;
  }

  // Interaction types
  double y_right = 0.92 - 1.5 * dy; // start on right column
  tex.SetTextFont(62);
  tex.DrawLatex(0.50, y_right, "Interaction breakdown:"); y_right -= dy;
  tex.SetTextFont(42);

  // Sort by key for consistent ordering
  std::vector<std::pair<std::string, size_t>> sorted_int(
      counts.by_inttype.begin(), counts.by_inttype.end());
  std::sort(sorted_int.begin(), sorted_int.end());
  for (auto &[key, n] : sorted_int) {
    tex.DrawLatex(0.53, y_right, Form("%-16s %zu", key.c_str(), n));
    y_right -= dy;
  }
}

// gStyle->SetPalette is GLOBAL, looked up by COLZ at canvas-print time —
// not at hist->Draw() time. For pages with multiple pads using different
// palettes, all pads end up using the LAST palette set unless we ship a
// per-pad TExec that re-sets the palette as the pad is rendered. The
// strings below are Cling-evaluated by TExec at that moment.
//
// TExec command strings (Cling-evaluated at pad render time). Initially try
// custom gradient tables; if those don't switch correctly per-pad on this
// ROOT build, fall back to stock palettes that Cling definitely knows.
inline const char *kPaletteBWRCmd =
  "double s[7]={0.00,0.17,0.33,0.50,0.67,0.83,1.00};"
  "double r[7]={0.129,0.263,0.573,0.969,0.957,0.839,0.698};"
  "double g[7]={0.400,0.576,0.773,0.969,0.647,0.376,0.094};"
  "double b[7]={0.674,0.764,0.871,0.969,0.510,0.302,0.169};"
  "TColor::CreateGradientColorTable(7,s,r,g,b,255);";

inline const char *kPaletteMagmaCmd =
  "double s[6]={0.00,0.20,0.40,0.60,0.80,1.00};"
  "double r[6]={0.001,0.234,0.550,0.866,0.987,0.987};"
  "double g[6]={0.000,0.060,0.149,0.219,0.471,0.991};"
  "double b[6]={0.014,0.402,0.506,0.420,0.299,0.750};"
  "TColor::CreateGradientColorTable(6,s,r,g,b,255);";

// Sequential single-hue (very light -> deep navy blue) for the CV panel.
// White-on-blue intensity scale; no second color, no hue rotation.
inline const char *kPaletteBlueSequentialCmd =
  "double s[5]={0.00,0.25,0.50,0.75,1.00};"
  "double r[5]={0.969,0.776,0.486,0.192,0.031};"
  "double g[5]={0.984,0.859,0.722,0.510,0.188};"
  "double b[5]={1.000,0.937,0.847,0.741,0.420};"
  "TColor::CreateGradientColorTable(5,s,r,g,b,255);";

// Compile-time helpers used outside TExec (e.g. one-time init in main).
// 7 stops, divergent blue->white->red. The midpoint stop is pure white
// (1.0, 1.0, 1.0) -- not the ColorBrewer 247/247/247 light gray -- so the
// central band on a 99-contour palette renders as true white and the
// colourbar tick at 0 visibly sits on white, not on faint pink.
inline void UseBlueWhiteRedPalette() {
  static Double_t stops[7]  = {0.00,  0.17,  0.33,  0.50, 0.67,  0.83,  1.00};
  static Double_t reds[7]   = {0.129, 0.263, 0.573, 1.0,  0.957, 0.839, 0.698};
  static Double_t greens[7] = {0.400, 0.576, 0.773, 1.0,  0.647, 0.376, 0.094};
  static Double_t blues[7]  = {0.674, 0.764, 0.871, 1.0,  0.510, 0.302, 0.169};
  TColor::CreateGradientColorTable(7, stops, reds, greens, blues, 255);
}
inline void UseMagmaPalette() {
  static Double_t stops[6]  = {0.00, 0.20, 0.40, 0.60, 0.80, 1.00};
  static Double_t reds[6]   = {0.001, 0.234, 0.550, 0.866, 0.987, 0.987};
  static Double_t greens[6] = {0.000, 0.060, 0.149, 0.219, 0.471, 0.991};
  static Double_t blues[6]  = {0.014, 0.402, 0.506, 0.420, 0.299, 0.750};
  TColor::CreateGradientColorTable(6, stops, reds, greens, blues, 255);
}

// Pick the index in `tweakvals` closest to `target` (e.g. +1 or -1).
// Returns -1 if the vector is empty.
inline int IndexClosestTo(const std::vector<double> &tweakvals, double target) {
  if (tweakvals.empty()) return -1;
  int best = 0;
  double best_d = std::abs(tweakvals[0] - target);
  for (int i = 1; i < (int)tweakvals.size(); ++i) {
    double d = std::abs(tweakvals[i] - target);
    if (d < best_d) { best_d = d; best = i; }
  }
  return best;
}

// ===== PDF plotting -- 2D mode =============================================
// Per (channel, dial): pairs are emitted two per page, each row is
//   [ (var_-1sigma - CV)/CV  |  CV 2D  |  (var_+1sigma - CV)/CV ]
// drawn with a diverging palette centered at 0 for the ratio panels and
// the default palette for the CV. Snaps to the nearest tweakval if exact
// +-1 sigma is not available.
void Make2DPlots(const Hist2DMap &allhists2d,
                 const std::vector<Var2DPair> &pairs,
                 size_t n_events_processed) {
  // Same monotonic-counter trick as MakePlots: avoids gDirectory name
  // collisions across (channel, dial, page, row) iterations that would
  // otherwise cause Clone() to silently reuse a stale histogram.
  static int g_panel_2d_uid = 0;
  std::string pdfname = cliopts::output_base + ".pdf";
  TCanvas *c = new TCanvas("c2dplots", "", 1400, 1000);
  // Suppress the default ROOT stat box on every histogram drawn from here
  // on. It otherwise sits over the upper-right corner and obscures the
  // y-axis labels in tight panels.
  gStyle->SetOptStat(0);

  // Force an ODD number of contour bands so the central band sits exactly
  // on palette position 0.5. ROOT's default 20 bands (even) makes the band
  // straddling 0 fetch its colour from position 10.5/20 = 0.525, which is
  // already past the white midpoint and on the salmon side -- the colourbar
  // tick at "0" then visibly lines up with light pink instead of white.
  // 99 bands gives a smooth gradient and a clean central band at 49.5/99 =
  // 0.500 = the BWR palette's pure-white midpoint.
  gStyle->SetNumberContours(99);

  // Same channel ordering helper as 1D.
  auto channel_priority = [](const std::string &s) -> std::string {
    if (s == "Total") return "000";
    if (s == "Total_CC") return "001";
    if (s == "Total_NC") return "002";
    return "100_" + s;
  };
  std::vector<std::string> channels;
  for (auto &[ch, _] : allhists2d) channels.push_back(ch);
  std::sort(channels.begin(), channels.end(),
            [&](const std::string &a, const std::string &b) {
              return channel_priority(a) < channel_priority(b);
            });

  c->Print((pdfname + "[").c_str());

  // Resolve x and y VarDefs for each pair (registry lookup) to build axis
  // titles in the rendering loop.
  auto resolve_pair_vds = [](Var2DPair const &p, VarDef &x_vd, VarDef &y_vd) {
    auto ix = g_vars_registry.find(p.x_var);
    auto iy = g_vars_registry.find(p.y_var);
    if (ix == g_vars_registry.end() || iy == g_vars_registry.end()) return false;
    x_vd = ix->second; y_vd = iy->second; return true;
  };

  double inv_n = n_events_processed > 0 ? 1.0 / double(n_events_processed) : 1.0;
  const int per_page = 2;  // pairs per page (2 rows of 3 panels)

  for (auto &ch : channels) {
    auto pit_ch = allhists2d.find(ch);
    if (pit_ch == allhists2d.end()) continue;
    for (auto &[parname, pair_map] : pit_ch->second) {
      int n_pairs = (int)pairs.size();
      int n_pages = (n_pairs + per_page - 1) / per_page;
      for (int ipage = 0; ipage < n_pages; ++ipage) {
        int pair_start = ipage * per_page;
        int pair_end   = std::min(pair_start + per_page, n_pairs);
        int nrows = pair_end - pair_start;

        c->Clear();
        c->cd();

        TPad *titlepad = new TPad(Form("tp2_%d", ipage), "", 0, 0.93, 1, 1.0);
        titlepad->SetFillColor(kWhite);
        titlepad->Draw(); titlepad->cd();
        TLatex title; title.SetNDC();
        title.SetTextSize(0.42); title.SetTextAlign(22); title.SetTextFont(62);
        std::string page_title = parname + "  |  " + MakeChannelLabel(ch);
        if (n_pages > 1) page_title += Form("  (%d/%d)", ipage + 1, n_pages);
        title.DrawLatex(0.5, 0.45, page_title.c_str());

        c->cd();
        double grid_top = 0.92, grid_bot = 0.04;
        double row_h = (grid_top - grid_bot) / per_page;
        // ncols hardcoded to 3 (sigma- | CV | sigma+)
        double col_w[3] = {1.0/3, 1.0/3, 1.0/3};

        for (int ip = pair_start; ip < pair_end; ++ip) {
          int row = ip - pair_start;
          double y_top = grid_top - row * row_h;
          double y_bot = y_top - row_h;
          // Leave small inner margin between top of pad and ceiling
          double pad_top = y_top - 0.005, pad_bot = y_bot + 0.005;

          auto const &pr = pairs[ip];
          auto sh_it = pair_map.find(pr.key);
          if (sh_it == pair_map.end()) continue;
          auto const &sh = sh_it->second;
          if (!sh.h_cv) continue;

          VarDef x_vd, y_vd;
          if (!resolve_pair_vds(pr, x_vd, y_vd)) continue;

          // Find indices closest to -1 and +1 sigma.
          int im = IndexClosestTo(sh.tweakvals, -1.0);
          int ip_sig = IndexClosestTo(sh.tweakvals, +1.0);
          double tw_m = (im >= 0 && im < (int)sh.tweakvals.size()) ? sh.tweakvals[im] : -1;
          double tw_p = (ip_sig >= 0 && ip_sig < (int)sh.tweakvals.size())
                          ? sh.tweakvals[ip_sig] : +1;

          // ---- left panel: -sigma ratio --------------------------------------
          double x1 = 0.0, x2 = col_w[0];
          c->cd();
          TPad *pL = new TPad(Form("p2L_%d_%d", ipage, row), "",
                              x1, pad_bot, x2, pad_top);
          pL->SetTopMargin(0.10); pL->SetBottomMargin(0.18);
          pL->SetLeftMargin(0.17); pL->SetRightMargin(0.16);
          pL->Draw(); pL->cd();
          if (im >= 0 && im < (int)sh.h_var.size() && sh.h_var[im]) {
            int uid2d = ++g_panel_2d_uid;
            TH2D *rL = (TH2D*)sh.h_var[im]->Clone(Form("rL_%d", uid2d));
            rL->Add(sh.h_cv, -1.0);
            rL->Divide(sh.h_cv);
            // Zero bins where either:
            //   - the CV denominator has no events (ratio meaningless;
            //     TH2::Divide can leak a coloured speckle on empty regions),
            //   - the CV is below ratio_min_cv_count * peak_CV (low-stats noise),
            //   - the ratio itself is non-finite (NaN / Inf safety net).
            // ratio_min_cv_count is a FRACTION OF PEAK (scale-invariant).
            {
              double _peak = sh.h_cv->GetMaximum();
              double _thresh = _peak * cliopts::ratio_min_cv_count;
              for (int bx = 1; bx <= rL->GetNbinsX(); ++bx) {
                for (int by = 1; by <= rL->GetNbinsY(); ++by) {
                  double _cv = sh.h_cv->GetBinContent(bx, by);
                  double _r  = rL->GetBinContent(bx, by);
                  if (_cv <= 0 || _cv < _thresh || !std::isfinite(_r)) {
                    rL->SetBinContent(bx, by, 0);
                    rL->SetBinError(bx, by, 0);
                  }
                }
              }
            }
            // Clamp the colourbar range to [0.05, 5.0] so a single low-stats
            // bin with a huge ratio doesn't drive zmax to e.g. 100 and squash
            // every other bin to near-white. Same logic as the 1D ratio panel.
            double zmax = std::max(std::abs(rL->GetMinimum()),
                                   std::abs(rL->GetMaximum()));
            zmax = std::min(std::max(zmax, 0.05), 5.0);
            rL->SetMinimum(-zmax); rL->SetMaximum(+zmax);
            rL->SetTitle(Form("(%.2f#sigma - CV)/CV", tw_m));
            rL->GetXaxis()->SetTitle(BuildXAxisTitle(x_vd).c_str());
            rL->GetYaxis()->SetTitle(BuildXAxisTitle(y_vd).c_str());
            rL->GetXaxis()->SetTitleSize(0.06); rL->GetYaxis()->SetTitleSize(0.06);
            rL->GetXaxis()->SetLabelSize(0.05); rL->GetYaxis()->SetLabelSize(0.05);
            rL->GetXaxis()->SetTitleOffset(1.1); rL->GetYaxis()->SetTitleOffset(1.2);
            UseBlueWhiteRedPalette();
            rL->Draw("COLZ0");
            // Overlay the axis frame on top of the COLZ fill -- on this ROOT
            // build, COLZ0 with SetMinimum/SetMaximum sometimes skips the
            // frame redraw, leaving the data region without its black
            // rectangular border. AXIS SAME re-paints axes + frame.
            rL->Draw("AXIS SAME");
          }

          // ---- middle panel: CV 2D ------------------------------------------
          x1 = col_w[0]; x2 = col_w[0] + col_w[1];
          c->cd();
          TPad *pM = new TPad(Form("p2M_%d_%d", ipage, row), "",
                              x1, pad_bot, x2, pad_top);
          pM->SetTopMargin(0.10); pM->SetBottomMargin(0.18);
          pM->SetLeftMargin(0.17); pM->SetRightMargin(0.16);
          pM->Draw(); pM->cd();
          int uid2d_cv = ++g_panel_2d_uid;
          TH2D *cv_scaled = (TH2D*)sh.h_cv->Clone(Form("cv_%d", uid2d_cv));
          // Match the 1D ScaleToDiffXsec convention: divide by N_events and
          // by bin area (Scale "width" handles both x and y widths).
          cv_scaled->Scale(inv_n, "width");
          cv_scaled->SetTitle("CV");
          cv_scaled->GetXaxis()->SetTitle(BuildXAxisTitle(x_vd).c_str());
          cv_scaled->GetYaxis()->SetTitle(BuildXAxisTitle(y_vd).c_str());
          cv_scaled->GetXaxis()->SetTitleSize(0.06);
          cv_scaled->GetYaxis()->SetTitleSize(0.06);
          cv_scaled->GetXaxis()->SetLabelSize(0.05);
          cv_scaled->GetYaxis()->SetLabelSize(0.05);
          cv_scaled->GetXaxis()->SetTitleOffset(1.1);
          cv_scaled->GetYaxis()->SetTitleOffset(1.2);
          // All three panels share the BWR palette: ROOT's multi-palette
          // TExec pattern doesn't pin per-pad palettes on this build, so
          // any second palette set later (e.g. for the right panel) would
          // overwrite the middle panel's at canvas-print time.
          UseBlueWhiteRedPalette();
          cv_scaled->Draw("COLZ0");

          // ---- right panel: +sigma ratio ------------------------------------
          x1 = col_w[0] + col_w[1]; x2 = 1.0;
          c->cd();
          TPad *pR = new TPad(Form("p2R_%d_%d", ipage, row), "",
                              x1, pad_bot, x2, pad_top);
          pR->SetTopMargin(0.10); pR->SetBottomMargin(0.18);
          pR->SetLeftMargin(0.17); pR->SetRightMargin(0.16);
          pR->Draw(); pR->cd();
          if (ip_sig >= 0 && ip_sig < (int)sh.h_var.size() && sh.h_var[ip_sig]) {
            int uid2d_R = ++g_panel_2d_uid;
            TH2D *rR = (TH2D*)sh.h_var[ip_sig]->Clone(Form("rR_%d", uid2d_R));
            rR->Add(sh.h_cv, -1.0);
            rR->Divide(sh.h_cv);
            // See the left-panel block above for rationale.
            for (int bx = 1; bx <= rR->GetNbinsX(); ++bx) {
              for (int by = 1; by <= rR->GetNbinsY(); ++by) {
                double _cv = sh.h_cv->GetBinContent(bx, by);
                double _r  = rR->GetBinContent(bx, by);
                double _thresh_R = sh.h_cv->GetMaximum() * cliopts::ratio_min_cv_count;
                if (_cv <= 0 || _cv < _thresh_R || !std::isfinite(_r)) {
                  rR->SetBinContent(bx, by, 0);
                  rR->SetBinError(bx, by, 0);
                }
              }
            }
            double zmax = std::max(std::abs(rR->GetMinimum()),
                                   std::abs(rR->GetMaximum()));
            zmax = std::min(std::max(zmax, 0.05), 5.0);
            rR->SetMinimum(-zmax); rR->SetMaximum(+zmax);
            rR->SetTitle(Form("(+%.2f#sigma - CV)/CV", tw_p));
            rR->GetXaxis()->SetTitle(BuildXAxisTitle(x_vd).c_str());
            rR->GetYaxis()->SetTitle(BuildXAxisTitle(y_vd).c_str());
            rR->GetXaxis()->SetTitleSize(0.06); rR->GetYaxis()->SetTitleSize(0.06);
            rR->GetXaxis()->SetLabelSize(0.05); rR->GetYaxis()->SetLabelSize(0.05);
            rR->GetXaxis()->SetTitleOffset(1.1); rR->GetYaxis()->SetTitleOffset(1.2);
            UseBlueWhiteRedPalette();
            rR->Draw("COLZ0");
            rR->Draw("AXIS SAME");
          }
        }

        c->Print(pdfname.c_str());
      }
    }
  }

  c->Print((pdfname + "]").c_str());
  delete c;
  std::cout << "PDF written to " << pdfname << std::endl;
}

// ===== PDF plotting (appends pages to already-open PDF) ====================
void MakePlots(const HistMap &allhists, const std::vector<VarDef> &vars) {
  std::string pdfname = cliopts::output_base + ".pdf";
  TCanvas *c = new TCanvas("cplots", "", 1400, 1000);

  // Sort channels: Total first, then Total_CC, Total_NC, then per-channel
  std::vector<std::string> channels;
  for (auto &[ch, _] : allhists) channels.push_back(ch);
  std::sort(channels.begin(), channels.end(),
            [](const std::string &a, const std::string &b) {
              // Priority: Total=0, Total_CC=1, Total_NC=2, rest sorted
              auto pri = [](const std::string &s) -> std::string {
                if (s == "Total") return "000";
                if (s == "Total_CC") return "001";
                if (s == "Total_NC") return "002";
                return "1" + ChannelSortKey(s);
              };
              return pri(a) < pri(b);
            });

  // Collect parameter names (ordered)
  std::set<std::string> param_set;
  for (auto &[ch, pm] : allhists)
    for (auto &[pname, _] : pm) param_set.insert(pname);
  std::vector<std::string> param_names(param_set.begin(), param_set.end());

  // One page per (parameter, channel)
  for (auto &parname : param_names) {
    for (auto &chkey : channels) {
      auto ch_it = allhists.find(chkey);
      if (ch_it == allhists.end()) continue;
      auto par_it = ch_it->second.find(parname);
      if (par_it == ch_it->second.end()) continue;
      auto &var_map = par_it->second;

      // Count valid (non-empty) variables for this page
      std::vector<const VarDef *> page_vars;
      for (auto &vd : vars) {
        auto vit = var_map.find(vd.branch);
        if (vit != var_map.end() && vit->second.h_cv->GetEntries() > 0)
          page_vars.push_back(&vd);
      }
      if (page_vars.empty()) continue;

      // Split into pages of max 6 plots each
      const int max_per_page = cliopts::plots_per_page;
      int total_vars = (int)page_vars.size();
      int npages = (total_vars + max_per_page - 1) / max_per_page;

      for (int ipage = 0; ipage < npages; ++ipage) {
      int var_start = ipage * max_per_page;
      int var_end = std::min(var_start + max_per_page, total_vars);
      int nvars_this_page = var_end - var_start;

      int ncols = std::min(nvars_this_page, 3);
      int nrows = (nvars_this_page + ncols - 1) / ncols;

      c->Clear();
      c->SetFillColor(kWhite);

      // Title pad
      c->cd();
      TPad *titlepad = new TPad(Form("tp_%d", ipage), "", 0, 0.93, 1, 1.0);
      titlepad->SetFillColor(kWhite);
      titlepad->Draw(); titlepad->cd();
      TLatex title;
      title.SetNDC(); title.SetTextSize(0.42); title.SetTextAlign(22); title.SetTextFont(62);
      std::string page_title = parname + "  |  " + MakeChannelLabel(chkey);
      if (npages > 1) page_title += Form("  (%d/%d)", ipage + 1, npages);
      title.DrawLatex(0.5, 0.45, page_title.c_str());

      // Legend pad (horizontal, below title)
      c->cd();
      static int legpad_id = 0;
      TPad *legpad = new TPad(Form("lp_%d", legpad_id++), "", 0, 0.88, 1, 0.93);
      legpad->SetFillColor(kWhite);
      legpad->Draw(); legpad->cd();

      // Get variation values from the first variable's SystHists
      const SystHists &sh_first = var_map.begin()->second;
      int nvar = (int)sh_first.h_var.size();
      int nleg = nvar + 1;

      TLegend *leg = new TLegend(0.01, 0.05, 0.99, 0.95);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      leg->SetTextSize(0.55); leg->SetNColumns(std::min(nleg, 10));

      // Dummy histograms for legend (unique names)
      static int hleg_uid = 0;
      TH1D *hleg_cv = new TH1D(Form("hleg_cv_%d", hleg_uid), "", 1, 0, 1);
      hleg_cv->SetLineColor(kBlack); hleg_cv->SetLineWidth(2);
      leg->AddEntry(hleg_cv, "CV", "l");

      for (int i = 0; i < nvar; ++i) {
        TH1D *hl = new TH1D(Form("hleg_%d_%d", hleg_uid, i), "", 1, 0, 1);
        int col_i, sty;
        double vval = (i < (int)sh_first.tweakvals.size()) ? sh_first.tweakvals[i] : 0;
        GetVarColor(vval, col_i, sty);
        hl->SetLineColor(col_i); hl->SetLineWidth(2); hl->SetLineStyle(sty);
        std::string lbl = (i < (int)sh_first.tweakvals.size()) ? Form("%.3g", sh_first.tweakvals[i]) : Form("var%d", i);
        leg->AddEntry(hl, lbl.c_str(), "l");
      }
      hleg_uid++;
      leg->Draw();

      // Grid of panels
      c->cd();
      // grid_bot bumped from 0.0 -> 0.04 so the bottom row's x-axis label
      // doesn't sit flush against the canvas edge and get clipped.
      double grid_top = 0.88, grid_bot = 0.04;
      double row_h = (grid_top - grid_bot) / nrows;
      double col_w = 1.0 / ncols;
      double ratio_frac = 0.30;

      for (int iv_abs = var_start; iv_abs < var_end; ++iv_abs) {
        int iv = iv_abs - var_start; // local index within this page
        const VarDef &vd = *page_vars[iv_abs];
        auto vit = var_map.find(vd.branch);
        if (vit == var_map.end()) continue;
        const SystHists &sh = vit->second;

        int col = iv % ncols;
        int row = iv / ncols;
        double x1 = col * col_w, x2 = (col + 1) * col_w;
        double y_top = grid_top - row * row_h;
        double y_bot = y_top - row_h;
        double y_split = y_bot + row_h * ratio_frac;

        // Pad and clone names use a monotonically-increasing global counter
        // so they never collide across pages, channels, or dials. ROOT
        // registers Clone'd histograms in gDirectory; TCanvas::Clear() does
        // not delete them, so reusing a name like "rh0_0" on a later page
        // causes the new Clone to silently reuse the OLD histogram (with
        // stale data) and the new ratio panel ends up empty. Was visible
        // most readily on multi-dial multi-channel runs through --plot-config.
        static int g_panel_uid = 0;
        int uid = ++g_panel_uid;

        // Main pad
        c->cd();
        TPad *pmain = new TPad(Form("m_%d", uid), "", x1, y_split, x2, y_top);
        // TopMargin bumped from 0.04 -> 0.09 to keep the y-axis title's
        // tallest characters (superscripts in d#sigma/dQ^{2}, etc.) inside
        // the pad. LeftMargin bumped 0.16 -> 0.18 for the same reason on
        // the title's outer edge.
        pmain->SetBottomMargin(0.01); pmain->SetTopMargin(0.09);
        pmain->SetLeftMargin(0.18); pmain->SetRightMargin(0.03);
        pmain->Draw(); pmain->cd();

        double ymax = sh.h_cv->GetMaximum();
        for (auto *h : sh.h_var) ymax = std::max(ymax, h->GetMaximum());

        TH1D *frame = (TH1D *)sh.h_cv->Clone(Form("f_%d", uid));
        frame->SetDirectory(nullptr);
        frame->Reset();
        frame->SetMaximum(ymax * 1.35); frame->SetMinimum(0);
        frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTickLength(0.03);
        frame->GetYaxis()->SetTitle(BuildYAxisTitle(vd).c_str());
        frame->GetYaxis()->SetTitleSize(0.06); frame->GetYaxis()->SetTitleOffset(1.25);
        frame->GetYaxis()->SetLabelSize(0.06);
        frame->Draw("AXIS");

        sh.h_cv->SetLineColor(kBlack); sh.h_cv->SetLineWidth(2);
        sh.h_cv->Draw("HIST SAME");

        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          int col_i, sty;
          double vval = (i < (int)sh.tweakvals.size()) ? sh.tweakvals[i] : 0;
          GetVarColor(vval, col_i, sty);
          sh.h_var[i]->SetLineColor(col_i);
          sh.h_var[i]->SetLineWidth(2);
          sh.h_var[i]->SetLineStyle(sty);
          sh.h_var[i]->Draw("HIST SAME");
        }

        // Find the variation with the largest chi^2 against CV; report chi^2
        // and NDF as raw numbers (NOT the chi^2/NDF ratio) so the reader can
        // form the ratio themselves. Also track the max fractional shift
        // across all variations. Display-only: histograms are NOT modified,
        // no bins are filtered out. Bins with h_cv == 0 are skipped from
        // the max fractional shift because the ratio is undefined there.
        double max_chi2 = 0;
        int ndf_at_max = 0;
        double max_frac_shift = 0;
        int nbins_x = sh.h_cv->GetNbinsX();
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          double chi2 = 0;
          int ndf = 0;
          int igood = 0;
          sh.h_cv->Chi2TestX(sh.h_var[i], chi2, ndf, igood, "WW");
          if (std::isfinite(chi2) && chi2 > max_chi2) {
            max_chi2 = chi2;
            ndf_at_max = ndf;
          }
          for (int b = 1; b <= nbins_x; ++b) {
            double cv_val = sh.h_cv->GetBinContent(b);
            if (cv_val == 0) continue;
            double var_val = sh.h_var[i]->GetBinContent(b);
            double frac = std::abs(var_val - cv_val) / std::abs(cv_val);
            if (frac > max_frac_shift) max_frac_shift = frac;
          }
        }

        // Draw stats box in upper-left of main pad. NDC y values stay inside
        // the data area (top margin is 0.09 -> data area top is NDC 0.91).
        // Text size kept small so the lines don't overrun the pad width when
        // there are 3 columns per page.
        TLatex stats;
        stats.SetNDC();
        stats.SetTextSize(0.045);
        stats.SetTextFont(42);
        stats.SetTextAlign(13); // top-left
        stats.DrawLatex(0.21, 0.88,
                         Form("#chi^{2}_{max}/ndf = %.2f / %d",
                              max_chi2, ndf_at_max));
        stats.DrawLatex(0.21, 0.82,
                         Form("max|#Delta|/CV = %.1f%%", max_frac_shift * 100));

        // Ratio pad
        c->cd();
        TPad *pratio = new TPad(Form("r_%d", uid), "", x1, y_bot, x2, y_split);
        pratio->SetTopMargin(0.01); pratio->SetBottomMargin(0.30);
        pratio->SetLeftMargin(0.18); pratio->SetRightMargin(0.03);
        pratio->Draw(); pratio->cd();

        TH1D *rframe = (TH1D *)sh.h_cv->Clone(Form("rf_%d", uid));
        rframe->SetDirectory(nullptr);
        rframe->Reset();
        // Adaptive ratio y-range from bins around the CV peak (where the
        // most events are, hence the most physically meaningful variation):
        // max |ratio| across every variation, restricted to bins with
        // CV >= 30% of peak CV. Drops low-stat tail bins that produce huge
        // ratios from a near-zero denominator. Clamped to [0.05, 5.0].
        double cv_peak = sh.h_cv->GetMaximum();
        double cv_floor = cv_peak * 0.30;
        double zmax_obs = 0.0;
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          for (int b = 1; b <= sh.h_var[i]->GetNbinsX(); ++b) {
            double cv = sh.h_cv->GetBinContent(b);
            if (cv <= 0 || cv < cv_floor) continue;
            double v  = sh.h_var[i]->GetBinContent(b);
            double r  = (v - cv) / cv;
            if (!std::isfinite(r)) continue;
            double ar = std::abs(r);
            if (ar > zmax_obs) zmax_obs = ar;
          }
        }
        double zmax = std::min(std::max(zmax_obs * 1.25, 0.05), 5.0);
        rframe->SetMaximum(+zmax); rframe->SetMinimum(-zmax);
        rframe->GetXaxis()->SetTitle(BuildXAxisTitle(vd).c_str());
        rframe->GetXaxis()->SetTitleSize(0.14); rframe->GetXaxis()->SetTitleOffset(0.9);
        rframe->GetXaxis()->SetLabelSize(0.12);
        rframe->GetYaxis()->SetTitle("#frac{Var-CV}{CV}");
        rframe->GetYaxis()->SetTitleSize(0.12); rframe->GetYaxis()->SetTitleOffset(0.45);
        rframe->GetYaxis()->SetLabelSize(0.10);
        rframe->GetYaxis()->SetNdivisions(505);
        rframe->Draw("AXIS");

        TLine *zero = new TLine(rframe->GetXaxis()->GetXmin(), 0,
                                 rframe->GetXaxis()->GetXmax(), 0);
        zero->SetLineStyle(2); zero->SetLineColor(kGray + 1); zero->Draw();

        // ratio_min_cv_count is interpreted as a FRACTION OF PEAK CV (not a
        // raw event count) so it's scale-invariant. After ScaleToDiffXsec
        // the h_cv bin contents are in differential-xsec units (~1e-38), so
        // comparing against a literal event-count threshold would zero
        // every bin and wipe the ratio lines.
        double cv_peak_for_supp = sh.h_cv->GetMaximum();
        double cv_supp_thresh = cv_peak_for_supp * cliopts::ratio_min_cv_count;
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          TH1D *rh = (TH1D *)sh.h_var[i]->Clone(Form("rh_%d_%d", uid, i));
          rh->SetDirectory(nullptr);
          rh->Add(sh.h_cv, -1.0);
          rh->Divide(sh.h_cv);
          // Suppress ratio bins where CV < ratio_min_cv_count * peak_CV.
          // Default ratio_min_cv_count = 0 means no suppression.
          if (cv_supp_thresh > 0) {
            for (int b = 1; b <= rh->GetNbinsX(); ++b) {
              if (sh.h_cv->GetBinContent(b) < cv_supp_thresh) {
                rh->SetBinContent(b, 0);
                rh->SetBinError(b, 0);
              }
            }
          }
          int col_i, sty;
          double vval = (i < (int)sh.tweakvals.size()) ? sh.tweakvals[i] : 0;
          GetVarColor(vval, col_i, sty);
          rh->SetLineColor(col_i); rh->SetLineWidth(2); rh->SetLineStyle(sty);
          rh->Draw("HIST SAME");
        }
      }
      c->Print(pdfname.c_str());
      } // end ipage loop
    }
  }
  delete c;
}

// ===== ROOT output =========================================================
void WriteROOT(const HistMap &allhists) {
  std::string rootname = cliopts::output_base + ".root";
  TFile *fout = new TFile(rootname.c_str(), "RECREATE");

  for (auto &[ch, param_map] : allhists) {
    TDirectory *chdir = fout->mkdir(ch.c_str());
    for (auto &[par, var_map] : param_map) {
      TDirectory *pdir = chdir->mkdir(par.c_str());
      for (auto &[var, sh] : var_map) {
        TDirectory *vdir = pdir->mkdir(var.c_str());
        vdir->cd();
        sh.h_cv->Write("h_cv");
        for (int i = 0; i < (int)sh.h_var.size(); ++i) {
          std::string hname = (i < (int)sh.tweakvals.size()) ?
            Form("h_var_%.3g", sh.tweakvals[i]) : Form("h_var_%d", i);
          sh.h_var[i]->Write(hname.c_str());
        }
      }
    }
  }

  // Save legend info as a TTree
  fout->cd();
  // (legend colors/styles are deterministic from GetVarColor, no need to store)

  fout->Close();
  std::cout << "ROOT file written to " << rootname << std::endl;
  delete fout;
}

// ===== main ================================================================
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

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  // Raise ROOT's diagnostic threshold before any rendering happens. Without
  // this, ROOT's TH1::Chi2TestX floods stdout with one "less than 10
  // effective events" Info per low-stats page. Was previously only inside
  // FillFromGHEP, so flat-tree mode missed it.
  nusyst::quiet::SetGlobalQuiet();
  // Hide the default ROOT stat box on every histogram. It otherwise sits in
  // the upper-right corner and overlays the y-axis labels in tight panels.
  gStyle->SetOptStat(0);

  if (cliopts::input_file.empty()) {
    std::cout << "[ERROR]: -i <input> is required." << std::endl;
    SayUsage(argv); return 1;
  }

  // Cache fallback: -p without -c means "use the cached kitchen sink and
  // filter to these dials". Same resolution + auto-generation as nusyst
  // inventory / tweaks / response.
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

  // Resolve variables
  std::vector<VarDef> vars;
  if (cliopts::variables.empty()) {
    EnsureRegistryInitialised();
    for (auto &[key, vd] : g_vars_registry) vars.push_back(vd);
  } else {
    for (auto &v : cliopts::variables) vars.push_back(ResolveVar(v));
  }

  // Detect input mode
  TFile *fin = TFile::Open(cliopts::input_file.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cout << "[ERROR]: Cannot open " << cliopts::input_file << std::endl;
    return 2;
  }

  bool flat_tree_mode = false;
  TTree *flat_tree = dynamic_cast<TTree *>(fin->Get(cliopts::treename.c_str()));
  TTree *meta_tree = dynamic_cast<TTree *>(fin->Get("tweak_metadata"));

  if (flat_tree && meta_tree) {
    TObjArray *branches = flat_tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i)
      if (std::string(branches->At(i)->GetName()).find("tweak_responses_") == 0)
        { flat_tree_mode = true; break; }
  }

  HistMap allhists;
  EventCounts counts;

  // Optional per-provider applies_to_channels patterns. Loaded only if a
  // fhicl config is given; empty -> no filtering.
  std::map<std::string, std::vector<std::string>> applies_to_channels;
  if (!cliopts::fclname.empty()) {
    applies_to_channels = LoadAppliesToChannelsMap(cliopts::fclname,
                                                    cliopts::fhicl_key);
  }

  if (flat_tree_mode) {
    std::cout << "Detected flat tree mode." << std::endl;

    // Read param metadata
    std::vector<ParamMeta> params;
    {
      TObjString *name_obj = nullptr;
      int ntw = 0; double tweakvals[200];
      meta_tree->SetBranchAddress("name", &name_obj);
      meta_tree->SetBranchAddress("ntweaks", &ntw);
      meta_tree->SetBranchAddress("tweakvalues", tweakvals);
      for (Long64_t i = 0; i < meta_tree->GetEntries(); ++i) {
        meta_tree->GetEntry(i);
        ParamMeta pm;
        pm.fullname = name_obj->GetString().Data();
        pm.ntweaks = ntw;
        pm.tweakvalues.assign(tweakvals, tweakvals + ntw);
        if (ParamSelected(pm.fullname)) params.push_back(pm);
      }
      meta_tree->ResetBranchAddresses();
    }

    // Filter variables to those with existing branches (skip coslep etc.)
    // Formula-backed vars skip this check -- they compose other branches via
    // TTreeFormula and don't have a single backing branch to test for.
    // 2D mode doesn't render 1D `vars`; skip the filter (and its warnings)
    // entirely.
    std::vector<VarDef> valid_vars;
    if (cliopts::dim == 2) {
      valid_vars = vars;
    } else
    for (auto &vd : vars) {
      if (!vd.formula.empty()) {
        valid_vars.push_back(vd);
        continue;
      }
      bool is_computed = (vd.branch == "leading_proton_KE" || vd.branch == "coslep");
      if (is_computed) {
        if (vd.branch == "leading_proton_KE" && flat_tree->GetBranch("fsprotons_KE"))
          valid_vars.push_back(vd);
        else if (vd.branch == "coslep")
          std::cout << "[INFO]: coslep not available in flat tree mode, skipping." << std::endl;
        else
          std::cout << "[WARN]: " << vd.branch << " not available, skipping." << std::endl;
      } else if (flat_tree->GetBranch(vd.branch.c_str())) {
        valid_vars.push_back(vd);
      } else {
        std::cout << "[WARN]: Branch " << vd.branch << " not found, skipping." << std::endl;
      }
    }
    vars = valid_vars;

    if (cliopts::dim == 2) {
      Hist2DMap allhists2d;
      FillFromFlatTree2D(flat_tree, meta_tree, kPredefined2DPairs, params,
                         applies_to_channels, allhists2d, counts);
      if (allhists2d.empty()) {
        std::cout << "[ERROR]: No 2D histograms were filled." << std::endl;
        return 4;
      }
      Make2DPlots(allhists2d, kPredefined2DPairs, counts.total);
      return 0;  // 2D path bypasses the remaining 1D rendering pipeline.
    }
    FillFromFlatTree(flat_tree, meta_tree, vars, params, applies_to_channels,
                     allhists, counts);
  } else {
    std::cout << "Detected GHEP mode." << std::endl;
    if (cliopts::dim == 2) {
      std::cout << "[ERROR]: --dim 2 requires a flat-tree input (the output\n"
                   "         of `nusyst tweaks`). Run that first:\n"
                   "           nusyst tweaks -p <dial> -i events.ghep.root "
                   "-o tweaks.root\n"
                   "         then pass tweaks.root as -i here." << std::endl;
      return 5;
    }
    // Formula vars need a TTree backing them; warn once and let the per-event
    // loop silently skip them.
    std::vector<std::string> dropped;
    for (auto const &vd : vars)
      if (!vd.formula.empty()) dropped.push_back(vd.branch);
    if (!dropped.empty()) {
      std::cerr << "[WARN]: GHEP mode cannot evaluate TTreeFormula variables; "
                   "skipping:\n        ";
      for (auto const &n : dropped) std::cerr << n << " ";
      std::cerr << "\n        Run `nusyst tweaks` first, then pass the "
                   "resulting events tree as -i to use formula vars."
                << std::endl;
    }
    fin->Close(); delete fin; fin = nullptr;
    if (cliopts::fclname.empty()) {
      std::cout << "[ERROR]: GHEP mode requires -c <config.fcl>" << std::endl;
      return 3;
    }
    FillFromGHEP(cliopts::input_file, cliopts::fclname, vars,
                 applies_to_channels, allhists, counts);
  }

  PruneEmpty(allhists);
  if (allhists.empty()) {
    std::cout << "[ERROR]: No histograms were filled." << std::endl;
    return 4;
  }

  ScaleToDiffXsec(allhists, counts.total);

  // Build inclusive channels: total CC, total NC, total
  // by merging per-channel histograms
  std::set<std::string> param_names;
  for (auto &[ch, pm] : allhists)
    for (auto &[pn, _] : pm) param_names.insert(pn);

  auto AddChannel = [&](const std::string &newkey,
                         std::function<bool(const std::string &)> matcher) {
    for (auto &pn : param_names) {
      for (auto &vd : vars) {
        TH1D *sum_cv = nullptr;
        std::vector<TH1D *> sum_var;
        std::vector<double> tweakvals;

        for (auto &[ch, pm] : allhists) {
          if (!matcher(ch)) continue;
          auto pit = pm.find(pn);
          if (pit == pm.end()) continue;
          auto vit = pit->second.find(vd.branch);
          if (vit == pit->second.end()) continue;
          auto &sh = vit->second;

          if (!sum_cv) {
            std::string tag = newkey + "_" + pn + "_" + vd.branch;
            sum_cv = (TH1D *)sh.h_cv->Clone(("h_cv_" + tag).c_str());
            tweakvals = sh.tweakvals;
            for (int i = 0; i < (int)sh.h_var.size(); ++i)
              sum_var.push_back((TH1D *)sh.h_var[i]->Clone(Form("h_var%d_%s", i, tag.c_str())));
          } else {
            sum_cv->Add(sh.h_cv);
            for (int i = 0; i < (int)sh.h_var.size() && i < (int)sum_var.size(); ++i)
              sum_var[i]->Add(sh.h_var[i]);
          }
        }
        if (sum_cv && sum_cv->GetEntries() > 0) {
          SystHists sh;
          sh.h_cv = sum_cv;
          sh.h_var = sum_var;
          sh.tweakvals = tweakvals;
          allhists[newkey][pn][vd.branch] = sh;
        }
      }
    }
  };

  AddChannel("Total", [](const std::string &) { return true; });
  AddChannel("Total_CC", [](const std::string &ch) { return ch.find("CC_") == 0; });
  AddChannel("Total_NC", [](const std::string &ch) {
    return ch.find("NC_") == 0 && ch.find("Total") == std::string::npos;
  });

  // Content-based prune: drop (channel, dial) pairs where every variation
  // histogram matches CV bin-for-bin. Catches dials whose internal GENIE cuts
  // (e.g. struck-nucleon PDG, FS pion multiplicity, FS-hadron presence for
  // FSI dials) reject every event in the channel even though the channel
  // topology nominally matched the dial's applies_to_channels patterns.
  PruneTrivialParams(allhists);

  if (cliopts::make_root) WriteROOT(allhists);

  if (cliopts::make_pdf) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    std::string pdfname = cliopts::output_base + ".pdf";
    TCanvas *c = new TCanvas("c", "", 1400, 1000);
    c->Print((pdfname + "[").c_str());

    // Summary page
    MakeSummaryPage(c, counts);
    c->Print(pdfname.c_str());

    // Now delegate to MakePlots for the actual histogram pages
    // (MakePlots writes pages to an already-open PDF)
    MakePlots(allhists, vars);

    c->Print((pdfname + "]").c_str());
    delete c;
  }

  if (fin) { fin->Close(); delete fin; }
  return 0;
}
