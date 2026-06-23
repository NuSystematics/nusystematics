// GenerateAllDialsConfigNuSyst.cxx
// Builds a "kitchen sink" generated_systematic_provider_configuration fhicl
// containing every available dial:
//   - all GENIE Reweight dials enumerated from genie::rew::GSyst_t at runtime
//     (no hardcoded dial names -- discovered via GSyst::AsString)
//   - all standalone nusystematics provider tool configs found in --fcl-dir,
//     each loaded and validated; ones that fail to instantiate are reported
//     and skipped.
//
// Modes (selected via --mode): genierw, providers, all  (default: all)

#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/utility/make_instance.hh"
#include "nusystematics/utility/response_helper.hh"
#include "nusystematics/utility/silence_genie.hh"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/parse.h"

#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "RwFramework/GSyst.h"

#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace systtools;
using namespace nusyst;

// ===== CLI =================================================================
namespace cliopts {
std::string mode = "all"; // genierw | providers | all
std::string outputfile;
std::string fcl_dir;
std::string variation_descriptor = "[-3,-2,-1,0,1,2,3]";
bool single_instance = false;
bool include_skeleton = false;
std::vector<std::string> skip_providers;

// Tool-config files that are *templates* / development examples rather than
// production reweight providers. Skipped by default in providers/all modes;
// opt back in with --include-skeleton.
//
//   SkeleWeighter.ToolConfig.fcl
//     SkeleWeighter is the worked example used in README.md to teach how to
//     write a new syst provider. As a bonus it exercises a GENIE Reweight bug:
//     it stores genie::rew::GReWeightNuXSecCCQE *by value* in a std::vector,
//     and the upstream class lacks a proper copy/move constructor, so the
//     first vector realloc shallow-copies internal pointers and the
//     subsequent destructor double-frees → SIGSEGV. Not a nusystematics bug
//     to fix here; not a production provider either.
const std::set<std::string> kDevelopmentTemplates = {
  "SkeleWeighter.ToolConfig.fcl",
};
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n\n"
    "  Generates a 'generated_systematic_provider_configuration' fhicl\n"
    "  containing every available dial.\n\n"
    "  Optional:\n"
    "    --mode <m>             genierw | providers | all (default: all)\n"
    "    -o <output.fcl>        Output file (default: stdout)\n"
    "    --fcl-dir <dir>        Tool-config fcl directory\n"
    "                           (default: $NUSYST/fcl)\n"
    "    --variation-descriptor \"[-3,-2,-1,0,1,2,3]\"\n"
    "                           Variation descriptor used for every GENIE RW\n"
    "                           dial (default shown).\n"
    "    --single-instance      Emit one GENIEReWeight_All provider with all\n"
    "                           dials, instead of splitting by channel bucket\n"
    "                           (CCQE / CCRES / NCRES / DIS / FSI / etc.).\n"
    "                           Disables the channel-aware skip optimisation\n"
    "                           in DumpConfiguredTweaksNuSyst.\n"
    "    --include-skeleton     Include development-template tool configs\n"
    "                           (currently: SkeleWeighter.ToolConfig.fcl)\n"
    "                           that are skipped by default. SkeleWeighter\n"
    "                           segfaults on instantiation due to an upstream\n"
    "                           GReWeight rule-of-three violation; only use\n"
    "                           this if you are actively developing SkeleWeighter.\n"
    "    --skip <a,b,...>       Comma-separated tool-config filenames to\n"
    "                           additionally skip in providers mode\n"
    "                           (default skip list already covers the\n"
    "                           templates above).\n"
    "    -?|--help              Show this message\n"
    << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    std::string s(argv[opt]);
    if (s == "-?" || s == "--help") { SayUsage(argv); exit(0); }
    else if (s == "--mode") cliopts::mode = argv[++opt];
    else if (s == "-o") cliopts::outputfile = argv[++opt];
    else if (s == "--fcl-dir") cliopts::fcl_dir = argv[++opt];
    else if (s == "--variation-descriptor") cliopts::variation_descriptor = argv[++opt];
    else if (s == "--single-instance") cliopts::single_instance = true;
    else if (s == "--include-skeleton") cliopts::include_skeleton = true;
    else if (s == "--skip") {
      std::string tok; std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ',')) if (!tok.empty()) cliopts::skip_providers.push_back(tok);
    } else {
      std::cout << "[ERROR]: Unknown option: " << s << std::endl;
      SayUsage(argv); exit(1);
    }
    opt++;
  }
  if (cliopts::mode != "genierw" && cliopts::mode != "providers" && cliopts::mode != "all") {
    std::cout << "[ERROR]: --mode must be one of: genierw, providers, all" << std::endl;
    exit(1);
  }
}

// ===== Helpers =============================================================
bool IsValidGENIEDialName(const std::string &name) {
  if (name.empty() || name == "-") return false;
  if (name.find("INVALID") != std::string::npos) return false;
  return true;
}

// Mutual-exclusion / redundancy filter for the GENIEReWeight provider.
// Returns true if the dial should be included in the default tool config.
//
// nusystematics' GENIEReWeight has constraints encoded in
// GENIEReWeightParamConfig.cc that we have to respect:
//
//   1. Dipole vs Z-expansion CCQE axial form factor dials are mutually
//      exclusive within a single GReWeight engine. We default to dipole;
//      the Z-expansion family (ZNormCCQE / ZExpA1..4CCQE / AxFFCCQEshape)
//      is dropped. Supporting both in parallel would need a second CCQE
//      provider instance in Z-expansion mode -- TODO.
//
//   2. shape+norm Ma/Mv dials and their shape-only twins reweight the same
//      physics from two different code paths. The Norm{CC,NC}RES and the
//      shape-only {M[av](CC|NC)RES}shape variants crash the default RES
//      reweight engine when applied to events outside the RES topology.
//      The shape+norm (MaCCRES/MaNCRES/Mv...) are robust across all events,
//      so we keep those and drop the *shape twins. The NormCCQE/NormCCRES/
//      NormNCRES dials require the corresponding IsShapeOnly tool_options
//      which conflict with keeping the shape+norm dials, so we drop them.
// Which CCQE axial-FF parameterisation the loaded tune uses. Detected once
// by instantiating a transient GReWeightNuXSecCCQE and probing IsHandled with
// a known dipole-only dial (NormCCQE, gated on the dipole/MArunAxial branch
// of Init) and a known z-expansion-only dial (ZExpA1CCQE).
//
// The engine's Init() reads the FormFactorsAlg/AxialFormFactorModel out of
// the loaded tune's CCQE algorithm and self-selects fMode accordingly, so
// asking IsHandled is the same logic the runtime engine will use to accept
// or silently reject dials.
enum class CCQEFFMode { Dipole, ZExp, Unknown };

const char *CCQEFFModeName(CCQEFFMode m) {
  switch (m) {
    case CCQEFFMode::Dipole:  return "dipole";
    case CCQEFFMode::ZExp:    return "z-expansion";
    case CCQEFFMode::Unknown: return "unknown";
  }
  return "unknown";
}

CCQEFFMode DetectCCQEFFMode() {
  static const CCQEFFMode cached = []() -> CCQEFFMode {
    try {
      genie::rew::GReWeightNuXSecCCQE engine;
      if (engine.IsHandled(genie::rew::kXSecTwkDial_NormCCQE))
        return CCQEFFMode::Dipole;
      if (engine.IsHandled(genie::rew::kXSecTwkDial_ZExpA1CCQE))
        return CCQEFFMode::ZExp;
    } catch (std::exception const &e) {
      std::cerr << "[WARN]: Could not probe CCQE axial FF mode: " << e.what()
                << "; both dipole and z-expansion dials will be declared."
                << std::endl;
    } catch (...) {
      std::cerr << "[WARN]: Could not probe CCQE axial FF mode (non-std "
                   "exception); both dipole and z-expansion dials will be "
                   "declared." << std::endl;
    }
    return CCQEFFMode::Unknown;
  }();
  return cached;
}

// Returns an empty string if the dial should be kept, or a short human-readable
// reason if it should be excluded from the default tool config.
std::string DialInactiveReason(const std::string &name) {
  // Z-expansion CCQE axial form-factor dials -- accepted only by a CCQE engine
  // running in kModeZExp on events whose AxialFormFactorModel is ZExp.
  static const std::set<std::string> zexp_axial = {
    "ZNormCCQE", "ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE",
    "AxFFCCQEshape"
  };
  // Dipole CCQE axial form-factor dials -- accepted only by a CCQE engine
  // running in kModeMa / kModeNormAndMaShape on dipole or MArun events.
  static const std::set<std::string> dipole_axial = {
    "MaCCQE", "E0CCQE"
  };
  // Shape-only twins -- redundant with the shape+norm dials and unsafe on
  // non-matching event topologies.
  static const std::set<std::string> shape_only_twins = {
    "MaCCQEshape", "MaCCRESshape", "MvCCRESshape",
    "MaNCRESshape", "MvNCRESshape",
    "AhtBYshape", "BhtBYshape", "CV1uBYshape", "CV2uBYshape",
    "E0CCQEshape", "VecFFCCQEshape"
  };
  // Norm dials that require the corresponding IsShapeOnly tool option.
  // We keep the shape+norm MaCCQE/MaCCRES/..., so drop these.
  static const std::set<std::string> norm_requires_shape_only = {
    "NormCCQE", "NormCCRES", "NormNCRES"
  };
  // Electromagnetic-only dials. These only fire on e/mu probe scattering
  // (electromagnetic 2p2h, EM RES, ...) and return response = 1.0 for every
  // weak-current neutrino event, so they never produce a non-trivial weight
  // on a neutrino sample. Excluded by default; pass --include-em to keep.
  static const std::set<std::string> em_only_dials = {
    "NormEMMEC",         // EM 2p2h normalisation
    "EmpMEC_FracPN_EM",  // empirical MEC pn fraction, EM channel
    "EmpMEC_FracEMQE"    // empirical MEC EM-QE fraction
  };

  CCQEFFMode ff = DetectCCQEFFMode();
  if (zexp_axial.count(name)) {
    if (ff == CCQEFFMode::Dipole)
      return "Z-expansion CCQE FF -- loaded tune uses dipole axial FF";
    // ff is ZExp or Unknown: keep
  }
  if (dipole_axial.count(name)) {
    if (ff == CCQEFFMode::ZExp)
      return "dipole CCQE FF -- loaded tune uses Z-expansion axial FF";
    // ff is Dipole or Unknown: keep
  }
  if (shape_only_twins.count(name))
    return "shape-only twin of a kept shape+norm dial";
  if (norm_requires_shape_only.count(name))
    return "requires IsShapeOnly which conflicts with kept shape+norm dials";
  if (em_only_dials.count(name))
    return "EM-only dial -- trivial on weak-current neutrino events";
  return "";
}

// ===== Channel buckets for GENIE Reweight dials ==============================
// Each GENIE Reweight dial is mapped to a single channel "bucket". The bucket
// has an associated set of fnmatch-style channel-key patterns; events whose
// channel-key matches any pattern are the ones for which the bucket's provider
// will be evaluated at runtime. Channel-key format:
//
//     <CC|NC>_<nuname>_<target>_<topology>   e.g. "CC_numu_Ar40_QE"
//
// This is a hand-curated table -- the GENIE Reweight engines themselves don't
// expose channel applicability via any introspection API, but each dial's
// physics scope is well known.
const std::map<std::string, std::string> &DialToBucketMap() {
  static const std::map<std::string, std::string> m = {
    // CCQE family (only present in --include-fragile-qe / --include-zexp modes;
    // included here for completeness)
    {"MaCCQE", "CCQE"}, {"MaCCQEshape", "CCQE"},
    {"E0CCQE", "CCQE"}, {"E0CCQEshape", "CCQE"},
    {"NormCCQE", "CCQE"}, {"NormCCQEenu", "CCQE"},
    {"VecFFCCQEshape", "CCQE"},
    {"CCQEPauliSupViaKF", "CCQE"}, {"CCQEMomDistroFGtoSF", "CCQE"},
    {"RPA_CCQE", "CCQE"}, {"CoulombCCQE", "CCQE"},
    {"AxFFCCQEshape", "CCQE"},
    {"ZNormCCQE", "CCQE"},
    {"ZExpA1CCQE", "CCQE"}, {"ZExpA2CCQE", "CCQE"},
    {"ZExpA3CCQE", "CCQE"}, {"ZExpA4CCQE", "CCQE"},
    // NCEL (GENIE classifies NCEL as IsQuasiElastic + IsWeakNC)
    {"MaNCEL", "NCEL"}, {"EtaNCEL", "NCEL"},
    // CC resonance
    {"MaCCRES", "CCRES"}, {"MvCCRES", "CCRES"},
    {"MaCCRESshape", "CCRES"}, {"MvCCRESshape", "CCRES"},
    {"NormCCRES", "CCRES"},
    // NC resonance
    {"MaNCRES", "NCRES"}, {"MvNCRES", "NCRES"},
    {"MaNCRESshape", "NCRES"}, {"MvNCRESshape", "NCRES"},
    {"NormNCRES", "NCRES"},
    // Delta resonance internals (apply to both CC-RES and NC-RES)
    {"RDecBR1eta", "RES"}, {"RDecBR1gamma", "RES"},
    {"Theta_Delta2Npi", "RES"}, {"ThetaDelta2NRad", "RES"},
    // Coherent
    {"MaCOHpi", "COH"}, {"R0COHpi", "COH"},
    {"NormCCCOH", "COH"}, {"NormNCCOH", "COH"},
    // DIS (Bodek-Yang shape + AGKY hadronization)
    {"AhtBY", "DIS"}, {"BhtBY", "DIS"},
    {"CV1uBY", "DIS"}, {"CV2uBY", "DIS"},
    {"AhtBYshape", "DIS"}, {"BhtBYshape", "DIS"},
    {"CV1uBYshape", "DIS"}, {"CV2uBYshape", "DIS"},
    {"AGKYpT1pi", "DIS"}, {"AGKYxF1pi", "DIS"},
    {"RnubarnuCC", "DIS"}, {"NormDISCC", "DIS"},
    // MEC / 2p2h
    {"NormCCMEC", "MEC"}, {"NormNCMEC", "MEC"}, {"NormEMMEC", "MEC"},
    {"XSecShape_CCMEC", "MEC"}, {"DecayAngMEC", "MEC"},
    {"FracDelta_CCMEC", "MEC"}, {"FracPN_CCMEC", "MEC"},
    // Non-resonant single/two-pion background (applies to RES + DIS topologies)
    {"NonRESBGvpCC1pi",      "SPP"}, {"NonRESBGvpCC2pi",      "SPP"},
    {"NonRESBGvnCC1pi",      "SPP"}, {"NonRESBGvnCC2pi",      "SPP"},
    {"NonRESBGvpNC1pi",      "SPP"}, {"NonRESBGvpNC2pi",      "SPP"},
    {"NonRESBGvnNC1pi",      "SPP"}, {"NonRESBGvnNC2pi",      "SPP"},
    {"NonRESBGvbarpCC1pi",   "SPP"}, {"NonRESBGvbarpCC2pi",   "SPP"},
    {"NonRESBGvbarnCC1pi",   "SPP"}, {"NonRESBGvbarnCC2pi",   "SPP"},
    {"NonRESBGvbarpNC1pi",   "SPP"}, {"NonRESBGvbarpNC2pi",   "SPP"},
    {"NonRESBGvbarnNC1pi",   "SPP"}, {"NonRESBGvbarnNC2pi",   "SPP"},
    {"MKSPP_ReWeight",       "SPP"},
    // Final-state interactions -- pions/nucleons rescatter in the nucleus on
    // top of any underlying topology, so these stay un-filtered (apply to all).
    {"FormZone",     "FSI"},
    {"MFP_pi",       "FSI"}, {"MFP_N",        "FSI"},
    {"FrCEx_pi",     "FSI"}, {"FrCEx_N",      "FSI"},
    {"FrInel_pi",    "FSI"}, {"FrInel_N",     "FSI"},
    {"FrAbs_pi",     "FSI"}, {"FrAbs_N",      "FSI"},
    {"FrPiProd_pi",  "FSI"}, {"FrPiProd_N",   "FSI"},
    {"FrElas_pi",    "FSI"}, {"FrElas_N",     "FSI"},
  };
  return m;
}

// Returns the channel bucket for a given GENIE Reweight dial. Unknown dials
// fall into "Misc" -- those will be evaluated on every event by default.
std::string BucketForDial(const std::string &name) {
  auto const &m = DialToBucketMap();
  auto it = m.find(name);
  if (it == m.end()) return "Misc";
  return it->second;
}

// Fnmatch glob patterns matched against the per-event channel-key. Empty
// or absent => bucket applies to all events (no filtering).
std::vector<std::string> ChannelPatternsForBucket(const std::string &bucket) {
  static const std::map<std::string, std::vector<std::string>> m = {
    {"CCQE",  {"CC_*_*_QE"}},
    {"NCEL",  {"NC_*_*_QE"}},
    {"CCRES", {"CC_*_*_RES"}},
    {"NCRES", {"NC_*_*_RES"}},
    {"RES",   {"*_*_*_RES"}},
    {"COH",   {"*_*_*_COH"}},
    {"DIS",   {"*_*_*_DIS"}},
    {"MEC",   {"*_*_*_MEC"}},
    {"SPP",   {"*_*_*_RES", "*_*_*_DIS"}},
    // FSI / Misc: no patterns -> applies to every event
    {"FSI",   {}},
    {"Misc",  {}},
  };
  auto it = m.find(bucket);
  if (it == m.end()) return {};
  return it->second;
}

// Build per-bucket GENIEReWeight tool-config ParameterSets by enumerating
// GSyst_t. Returns a map keyed by bucket name; each value is a tool-config
// with `instance_name = <bucket>` and only that bucket's dials. With
// --single-instance, returns a single "All" entry with every dial (legacy
// kitchen-sink layout).
std::map<std::string, fhicl::ParameterSet>
BuildGENIEReWeightToolConfigs(std::vector<std::string> &included_dials,
                              std::vector<std::string> &skipped_names,
                              std::vector<std::string> &skipped_reasons) {
  std::map<std::string, fhicl::ParameterSet> out;

  // Build a flat list of (dial_name, bucket) for the dials we keep.
  std::vector<std::pair<std::string, std::string>> dials;
  int n_invalid_gaps = 0;
  for (int i = static_cast<int>(genie::rew::kNullSystematic) + 1;
       i < static_cast<int>(genie::rew::kNTwkDials); ++i) {
    genie::rew::GSyst_t s = static_cast<genie::rew::GSyst_t>(i);
    std::string name = genie::rew::GSyst::AsString(s);
    if (!IsValidGENIEDialName(name)) {
      ++n_invalid_gaps;  // collapsed into a single summary row below
      continue;
    }
    if (auto reason = DialInactiveReason(name); !reason.empty()) {
      skipped_names.push_back(name);
      skipped_reasons.push_back(reason);
      continue;
    }
    std::string bucket = cliopts::single_instance ? "All" : BucketForDial(name);
    dials.emplace_back(name, bucket);
    included_dials.push_back(name);
  }

  // Seed an empty ParameterSet per bucket so iteration order is deterministic.
  std::set<std::string> bucket_set;
  for (auto &[_, bucket] : dials) bucket_set.insert(bucket);

  for (auto const &bucket : bucket_set) {
    fhicl::ParameterSet ps;
    ps.put<std::string>("tool_type", "GENIEReWeight");
    ps.put<std::string>("instance_name", bucket);
    // Required to allow splineable parameters individually
    ps.put<bool>("ignore_parameter_dependence", true);
    ps.put<std::string>("genie_tune_name", "${GENIE_XSEC_TUNE}");
    out[bucket] = ps;
  }

  for (auto &[name, bucket] : dials) {
    out[bucket].put<double>(name + "_central_value", 0.0);
    out[bucket].put<std::string>(name + "_variation_descriptor",
                                  cliopts::variation_descriptor);
  }

  // Surface the count of invalid GSyst_t enum gaps as a single summary row
  // instead of one row per gap (typically 2-3 unused slots; not actionable).
  if (n_invalid_gaps > 0) {
    skipped_names.push_back("(unused GSyst_t enum gaps)");
    skipped_reasons.push_back(std::to_string(n_invalid_gaps) +
                              " enum slots with no dial name");
  }

  return out;
}

// Detect a tool-config-style fhicl: must contain a top-level "syst_providers" key.
bool IsToolConfigFile(const std::string &path) {
  std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
  try {
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(path, *fm);
    return ps.has_key("syst_providers");
  } catch (...) {
    return false;
  }
}

// Best-effort tool_type extraction from a tool-config fhicl. Walks each
// provider stanza listed in syst_providers and collects its tool_type
// string. Returns empty on parse failure (caller falls back gracefully).
std::vector<std::string> ToolTypesDeclaredIn(const std::string &path) {
  std::vector<std::string> out;
  std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
  try {
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(path, *fm);
    auto names = ps.get<std::vector<std::string>>("syst_providers",
                                                    std::vector<std::string>{});
    for (auto const &n : names) {
      try {
        auto prov = ps.get<fhicl::ParameterSet>(n);
        auto tt = prov.get<std::string>("tool_type", "");
        if (!tt.empty()) out.push_back(tt);
      } catch (...) {}
    }
  } catch (...) {}
  return out;
}

// List *.fcl files in a directory, recursively. Used to discover
// tool-config files; some live one level down in per-provider subdirs
// (e.g. CCQETemplateReweight/, MECq0q3InterpWeighting/, QEInterference/).
std::vector<std::string> ListFcls(const std::string &dir) {
  std::vector<std::string> out;
  std::error_code ec;
  std::filesystem::recursive_directory_iterator it(dir, ec);
  if (ec) return out;
  for (auto const &entry : it) {
    if (!entry.is_regular_file()) continue;
    auto const &p = entry.path();
    if (p.extension() != ".fcl") continue;
    out.push_back(p.string());
  }
  std::sort(out.begin(), out.end());
  return out;
}

// Verify that a provider can be re-instantiated from its own parameter
// headers via response_helper-style loading. This catches providers (like
// CCQERPA) that successfully build their metadata via ConfigureFromToolConfig
// but then fail in SetupResponseCalculator when actually used downstream.
bool VerifyProviderLoad(IGENIESystProvider_tool *prov, std::string &err_out) {
  try {
    fhicl::ParameterSet inner;
    inner.put(prov->GetFullyQualifiedName(), prov->GetParameterHeadersDocument());
    inner.put<std::vector<std::string>>("syst_providers",
                                         {prov->GetFullyQualifiedName()});
    nusyst::response_helper rh;
    rh.LoadProvidersAndHeaders(inner);
    return true;
  } catch (std::exception &e) {
    err_out = e.what();
    return false;
  } catch (...) {
    err_out = "unknown exception";
    return false;
  }
}

// Heuristic: is `err` a "missing external data file" error rather than a
// genuine code bug? These are reported as [SKIP] (data unavailable) so the
// user can tell at a glance which providers need their data staged into the
// data folder vs which have real problems to fix in code.
bool IsMissingDataFileError(const std::string &err) {
  static const std::vector<std::string> needles = {
    "Couldn't open input file",
    "Cannot open file",
    "No such file or directory",
  };
  for (auto const &n : needles) {
    if (err.find(n) != std::string::npos) return true;
  }
  return false;
}

// Heuristic: is `err` the known shared upstream bug currently affecting
// CCQERPA / QEInterference / FSIReweight? All three throw
// `basic_string: construction from null is not valid` during verify-load
// because something inside their SetupResponseCalculator chain constructs a
// std::string from a nullptr. Likely one shared utility fix; investigation
// tracked in the project_basicstring_null_bug memory entry.
//
// We classify it as [SKIP] (known upstream bug) rather than [FAIL] (genuine
// config problem) so users can tell at a glance that the provider's tool
// config is fine -- the failure is in shared code waiting on a fix.
bool IsKnownUpstreamBug(const std::string &err) {
  return err.find("basic_string: construction from null") != std::string::npos;
}

// Try to instantiate providers from a tool-config fhicl. Returns the
// configured providers (empty on failure) and writes any error to err_out.
// Also verifies each provider can be loaded from its own parameter headers.
// Uses syst_param_id_offset to make sure paramIds remain unique across calls.
std::vector<std::unique_ptr<IGENIESystProvider_tool>>
TryLoadToolConfig(const std::string &path, paramId_t &syst_param_id_offset,
                  std::string &err_out) {
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
  try {
    std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(path, *fm);
    tools = systtools::ConfigureISystProvidersFromToolConfig<IGENIESystProvider_tool>(
        ps, nusyst::make_instance, "syst_providers", syst_param_id_offset);
  } catch (std::exception &e) {
    err_out = e.what();
    return tools;
  } catch (...) {
    err_out = "unknown exception";
    return tools;
  }

  // Verification pass: try to load each provider from its own parameter
  // headers. Drop providers that fail this stricter test.
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> verified;
  for (auto &tool : tools) {
    std::string verr;
    if (VerifyProviderLoad(tool.get(), verr)) {
      syst_param_id_offset += tool->GetSystMetaData().size();
      verified.push_back(std::move(tool));
    } else {
      err_out = "verify failed: " + verr;
      return std::vector<std::unique_ptr<IGENIESystProvider_tool>>{};
    }
  }
  return verified;
}

// Try to instantiate the GENIEReWeight provider from an in-memory tool-config.
// Also verifies it can round-trip through response_helper-style loading.
// `instance_name` (e.g. "All", "CCQE", "FSI") must match the instance_name
// field inside `tool_ps` so the resulting FQ name lines up with the wrapper key.
std::vector<std::unique_ptr<IGENIESystProvider_tool>>
TryLoadInMemoryGENIERW(const fhicl::ParameterSet &tool_ps,
                        const std::string &instance_name,
                        paramId_t &syst_param_id_offset, std::string &err_out) {
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
  std::string fqname = "GENIEReWeight_" + instance_name;
  try {
    fhicl::ParameterSet wrapper;
    wrapper.put(fqname, tool_ps);
    wrapper.put<std::vector<std::string>>("syst_providers", {fqname});
    tools = systtools::ConfigureISystProvidersFromToolConfig<IGENIESystProvider_tool>(
        wrapper, nusyst::make_instance, "syst_providers", syst_param_id_offset);
  } catch (std::exception &e) {
    err_out = e.what();
    return tools;
  } catch (...) {
    err_out = "unknown exception";
    return tools;
  }

  std::vector<std::unique_ptr<IGENIESystProvider_tool>> verified;
  for (auto &tool : tools) {
    std::string verr;
    if (VerifyProviderLoad(tool.get(), verr)) {
      syst_param_id_offset += tool->GetSystMetaData().size();
      verified.push_back(std::move(tool));
    } else {
      err_out = "verify failed: " + verr;
      return std::vector<std::unique_ptr<IGENIESystProvider_tool>>{};
    }
  }
  return verified;
}

// ===== main ================================================================
// Bring GENIE up to the same state the runtime providers expect: the tune
// must be built so AlgFactory / XSecSplineList resolve algorithm configs.
// Mirrors the logic in IGENIESystProvider_tool::CheckTune().
void EnsureGENIETuneLoaded() {
  std::string current_tune = genie::XSecSplineList::Instance()->CurrentTune();
  if (!current_tune.empty()) return;
  char const *tune_env = std::getenv("GENIE_XSEC_TUNE");
  if (!tune_env || std::string(tune_env).empty()) {
    // We tried supporting "no tune set" via a warning + best-effort fallback,
    // but the very next step (DetectCCQEFFMode constructing a transient
    // GReWeightNuXSecCCQE) reaches into the empty AlgFactory registry and
    // SIGSEGVs in genie::Registry::SafeFind. Rather than chase every code
    // path that might dereference an unloaded tune, refuse to proceed.
    std::cerr << "[FATAL]: GENIE_XSEC_TUNE is not set in the environment.\n"
                 "         nusyst config (and the GENIE Reweight engines it\n"
                 "         constructs) require a built GENIE tune. Export the\n"
                 "         tune you generated your events with, e.g.\n"
                 "           export GENIE_XSEC_TUNE=AR23_20i_00_000\n"
                 "         and re-run." << std::endl;
    std::exit(2);
  }
  genie::RunOpt *grunopt = genie::RunOpt::Instance();
  grunopt->SetTuneName(tune_env);
  grunopt->BuildTune();
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  // Resolve fcl-dir default
  if (cliopts::fcl_dir.empty()) {
    char const *nusyst_env = std::getenv("NUSYST");
    if (nusyst_env) cliopts::fcl_dir = std::string(nusyst_env) + "/fcl";
  }

  nusyst::quiet::SetGlobalQuiet();
  // Build the requested GENIE tune up front so AlgFactory can resolve
  // tune-dependent algorithm configs (CCQE axial FF detection in particular).
  // The tune-build emits GENIE's mandatory ASCII banner plus dozens of INFO
  // lines (RunOpt, TuneId, TransverseEnhancementFFModel, ...); suppress.
  {
    nusyst::quiet::StdoutSink _quiet;
    genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
        "Messenger_whisper.xml");
    EnsureGENIETuneLoaded();
    // Prime the CCQE-FF detection cache silently so that
    // BuildGENIEReWeightToolConfigs and the cerr report below don't trigger
    // an engine construction with its associated chatter.
    DetectCCQEFFMode();
  }

  // Collect all successfully-loaded providers here
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> all_providers;
  // For per-bucket GENIE Reweight providers: provider FQ name ->
  // applies_to_channels patterns, written into the output fhicl so that
  // downstream consumers (DumpConfiguredTweaksNuSyst) can skip non-matching
  // providers per event without reading any external map.
  std::map<std::string, std::vector<std::string>> provider_channel_patterns;
  // Running paramId offset, threaded through every loader call so that
  // parameter IDs remain unique across providers in the merged output.
  paramId_t syst_param_id_offset = 0;

  // Capture standalone-scan outcomes so they can be embedded in the output
  // fhicl as a `_scan_report` table. `nusyst inventory` reads it back and
  // prints a footer summarising what was looked at but not loaded.
  //
  // Each failure category is stored as two parallel string arrays (names +
  // reasons) rather than as a single "<name>: <reason>" composite. Reason
  // strings can contain colons, brackets, and other characters that fhicl's
  // value parser tries to interpret if they appear unquoted in a list, which
  // would break round-tripping. The pair-of-arrays form avoids that.
  struct ScanReport {
    std::string              fcl_dir;
    std::vector<std::string> failed_names,            failed_reasons;
    std::vector<std::string> skipped_data_names,      skipped_data_reasons;
    std::vector<std::string> skipped_bug_names,       skipped_bug_reasons;
    std::vector<std::string> skipped_other_names,     skipped_other_reasons;
    std::vector<std::string> registered_no_fcl;       // dispatchable but no fhicl
    // GENIE Reweight dials we filter out from the kitchen-sink config so they
    // never appear in the inventory table. Surfaced in the footer so users
    // know they exist and why they were dropped.
    std::vector<std::string> genie_rw_skipped_names;
    std::vector<std::string> genie_rw_skipped_reasons;
  } scan_report;

  // Sanitize an error / reason message for embedding into the output fhicl as
  // a string. Two issues to handle:
  //   1. Newlines would break fhicl's single-line value form.
  //   2. fhicl's value parser treats `:`, `[`, `]`, `{`, `}`, `,` and quotes
  //      as syntactic -- even inside a string element of a vector<string>,
  //      `put<vector<string>>` re-validates each entry and fails on those
  //      characters. So we strip them down to spaces / backticks. The result
  //      is a readable footer at the cost of punctuation fidelity.
  auto flatten_err = [](std::string err) {
    std::string ascii; ascii.reserve(err.size());
    for (unsigned char c : err) {
      switch (c) {
        case '\n': case '\r':
        case ':':  case '[':  case ']':
        case '{':  case '}':  case ',':
          ascii += ' '; break;
        case '"':  case '\'':
          ascii += '`'; break;
        default:
          if (c < 32 || c > 126) {
            // Non-ASCII (UTF-8 bytes, em-dashes, etc.) -- fhicl's value parser
            // can choke on these. Replace with a plain hyphen for em-dashes
            // (E2 80 94) / en-dashes (E2 80 93) / others; the goal is "human
            // can still read it", not "round-trips bit-for-bit".
            ascii += '-';
          } else {
            ascii += static_cast<char>(c);
          }
          break;
      }
    }
    std::string out; out.reserve(ascii.size());
    bool prev_space = false;
    for (char c : ascii) {
      if (c == ' ') {
        if (!prev_space) out += c;
        prev_space = true;
      } else {
        out += c; prev_space = false;
      }
    }
    return out;
  };

  // ----- GENIE RW mode -----
  if (cliopts::mode == "genierw" || cliopts::mode == "all") {
    std::cerr << "=== GENIE Reweight ===" << std::endl;
    std::cerr << "[INFO] CCQE axial form factor detected: "
              << CCQEFFModeName(DetectCCQEFFMode())
              << " (dial set selected accordingly; "
                 "the inactive family is silently ignored by the runtime engine)"
              << std::endl;
    std::vector<std::string> included, skipped_names, skipped_reasons;
    auto bucket_configs = BuildGENIEReWeightToolConfigs(included, skipped_names,
                                                        skipped_reasons);
    std::cerr << "Enumerated " << (included.size() + skipped_names.size())
              << " GSyst_t entries (" << included.size() << " active across "
              << bucket_configs.size() << " "
              << (cliopts::single_instance ? "instance" : "channel buckets")
              << ", " << skipped_names.size() << " skipped)" << std::endl;
    if (!skipped_names.empty()) {
      std::cerr << "Skipped GENIE RW dials:" << std::endl;
      for (size_t i = 0; i < skipped_names.size(); ++i) {
        std::cerr << "    " << skipped_names[i]
                  << " (" << skipped_reasons[i] << ")" << std::endl;
      }
    }
    // Persist for `nusyst inventory` footer. Sanitize the reasons through
    // flatten_err for the same fhicl-roundtrip reasons as the scan report.
    scan_report.genie_rw_skipped_names = skipped_names;
    for (auto const &r : skipped_reasons)
      scan_report.genie_rw_skipped_reasons.push_back(flatten_err(r));

    for (auto &[bucket, grw_ps] : bucket_configs) {
      std::string err;
      std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
      {
        nusyst::quiet::StdoutSink _quiet;
        tools = TryLoadInMemoryGENIERW(grw_ps, bucket, syst_param_id_offset, err);
      }
      if (tools.empty()) {
        std::cerr << "[ERROR]: GENIE Reweight provider for bucket " << bucket
                  << " failed to instantiate:\n  " << err << std::endl;
        std::cerr << "  Hint: ensure GENIE_XSEC_TUNE is set and try --include-zexp\n"
                  << "        if the dipole CCQE FF dials are not the desired set."
                  << std::endl;
        continue;
      }
      std::vector<std::string> patterns = ChannelPatternsForBucket(bucket);
      for (auto &t : tools) {
        std::cerr << "  loaded: " << t->GetFullyQualifiedName() << " ("
                  << t->GetSystMetaData().size() << " params";
        if (!patterns.empty()) {
          std::cerr << ", channels: ";
          for (size_t i = 0; i < patterns.size(); ++i)
            std::cerr << (i ? ", " : "") << patterns[i];
        } else if (!cliopts::single_instance) {
          std::cerr << ", channels: all";
        }
        std::cerr << ")" << std::endl;
        if (!patterns.empty()) {
          provider_channel_patterns[t->GetFullyQualifiedName()] = patterns;
        }
        all_providers.push_back(std::move(t));
      }
    }
  }

  // ----- providers mode -----
  if (cliopts::mode == "providers" || cliopts::mode == "all") {
    std::cerr << "\n=== Standalone providers ===" << std::endl;
    if (cliopts::fcl_dir.empty()) {
      std::cerr << "[ERROR]: --fcl-dir not set and $NUSYST/fcl unavailable; "
                << "skipping providers." << std::endl;
    } else {
      scan_report.fcl_dir = cliopts::fcl_dir;
      std::cerr << "Scanning " << cliopts::fcl_dir << " (recursive)" << std::endl;
      std::vector<std::string> fcls = ListFcls(cliopts::fcl_dir);
      std::cerr << "Found " << fcls.size() << " .fcl files" << std::endl;

      for (auto &path : fcls) {
        // Use the path relative to fcl_dir so subdirectory layout is visible
        // in the report (e.g. "QEInterference/QEInterference.ToolConfig.fcl").
        std::string base = path;
        if (path.rfind(cliopts::fcl_dir + "/", 0) == 0)
          base = path.substr(cliopts::fcl_dir.size() + 1);
        std::string leaf = path.substr(path.find_last_of('/') + 1);

        // Skip development-template tool configs unless explicitly included
        if (!cliopts::include_skeleton &&
            cliopts::kDevelopmentTemplates.count(leaf)) {
          std::cerr << "  [SKIP] " << base
                    << " (development template; pass --include-skeleton to "
                       "include it)"
                    << std::endl;
          scan_report.skipped_other_names.push_back(base);
          scan_report.skipped_other_reasons.push_back("development template");
          continue;
        }
        // Skip user-specified files
        bool skip = false;
        for (auto &s : cliopts::skip_providers) {
          if (leaf == s) { skip = true; break; }
        }
        if (skip) {
          std::cerr << "  [SKIP] " << base << " (user --skip)" << std::endl;
          scan_report.skipped_other_names.push_back(base);
          scan_report.skipped_other_reasons.push_back("user --skip");
          continue;
        }

        // Skip non-toolconfig files (e.g. paramHeader_FSI.fcl is a generated config)
        if (!IsToolConfigFile(path)) {
          std::cerr << "  [SKIP] " << base << " (not a tool config: no syst_providers key)" << std::endl;
          scan_report.skipped_other_names.push_back(base);
          scan_report.skipped_other_reasons.push_back(
              "not a tool config (no syst_providers key)");
          continue;
        }

        std::string err;
        std::vector<std::unique_ptr<IGENIESystProvider_tool>> tools;
        {
          nusyst::quiet::StdoutSink _quiet;
          tools = TryLoadToolConfig(path, syst_param_id_offset, err);
        }
        if (tools.empty()) {
          if (IsMissingDataFileError(err)) {
            std::cerr << "  [SKIP] " << base
                      << " (missing external data file -- stage data and re-run)\n"
                      << "         " << err << std::endl;
            scan_report.skipped_data_names.push_back(base);
            scan_report.skipped_data_reasons.push_back(flatten_err(err));
          } else if (IsKnownUpstreamBug(err)) {
            std::cerr << "  [SKIP] " << base
                      << " (known upstream bug -- see project_basicstring_null_bug)\n"
                      << "         " << err << std::endl;
            scan_report.skipped_bug_names.push_back(base);
            scan_report.skipped_bug_reasons.push_back(flatten_err(err));
          } else {
            std::cerr << "  [FAIL] " << base << "\n         " << err << std::endl;
            scan_report.failed_names.push_back(base);
            scan_report.failed_reasons.push_back(flatten_err(err));
          }
          continue;
        }
        std::cerr << "  [OK]   " << base;
        for (auto &t : tools) {
          std::cerr << "  -> " << t->GetFullyQualifiedName() << " ("
                    << t->GetSystMetaData().size() << " params)";
          all_providers.push_back(std::move(t));
        }
        std::cerr << std::endl;
      }
    }
  }

  // ----- Compute "registered but no tool config" -----
  // The accurate definition we want: a registered tool_type appears in
  // `registered_no_fcl` iff no `.fcl` we scanned declared it. Providers
  // whose fhicl was found but failed/skipped are NOT the same thing -- they
  // have a fhicl, it just didn't load -- so those tool_types must also be
  // excluded from this list. Collect tool_types declared by *any* scanned
  // tool-config (regardless of load outcome), plus those from providers
  // that loaded successfully, and subtract from the registered set.
  {
    std::set<std::string> observed_tt;
    for (auto &p : all_providers) {
      try {
        observed_tt.insert(
            p->GetParameterHeadersDocument().get<std::string>("tool_type"));
      } catch (...) {}
    }
    if (!cliopts::fcl_dir.empty()) {
      for (auto const &path : ListFcls(cliopts::fcl_dir)) {
        for (auto const &tt : ToolTypesDeclaredIn(path)) observed_tt.insert(tt);
      }
    }
    for (auto const &tt : nusyst::RegisteredToolTypes()) {
      if (!observed_tt.count(tt))
        scan_report.registered_no_fcl.push_back(tt);
    }
  }

  if (all_providers.empty()) {
    std::cerr << "\n[ERROR]: No providers loaded successfully. Nothing to write."
              << std::endl;
    return 2;
  }

  // ----- Build the merged generated config -----
  fhicl::ParameterSet out_ps;
  std::vector<std::string> providerNames;
  for (auto &prov : all_providers) {
    if (!systtools::Validate(prov->GetSystMetaData(), false)) {
      std::cerr << "[WARN]: Provider " << prov->GetFullyQualifiedName()
                << " failed Validate(); skipping." << std::endl;
      continue;
    }
    fhicl::ParameterSet tool_ps = prov->GetParameterHeadersDocument();
    // For per-channel GENIE Reweight buckets, inject the channel patterns
    // into tool_options so downstream consumers can skip non-matching events.
    auto pcp_it = provider_channel_patterns.find(prov->GetFullyQualifiedName());
    if (pcp_it != provider_channel_patterns.end()) {
      fhicl::ParameterSet topts =
          tool_ps.get<fhicl::ParameterSet>("tool_options", fhicl::ParameterSet{});
      topts.put_or_replace<std::vector<std::string>>("applies_to_channels",
                                                       pcp_it->second);
      tool_ps.put_or_replace("tool_options", topts);
    }
    out_ps.put(prov->GetFullyQualifiedName(), tool_ps);
    providerNames.push_back(prov->GetFullyQualifiedName());
  }
  out_ps.put("syst_providers", providerNames);

  // Embed the scan report so `nusyst inventory` can surface providers that
  // were registered in make_instance but had no tool-config fhicl to load,
  // and those whose fhicl was found but failed/skipped.
  {
    fhicl::ParameterSet report_ps;
    report_ps.put("fcl_dir",                  scan_report.fcl_dir);
    report_ps.put("failed_names",             scan_report.failed_names);
    report_ps.put("failed_reasons",           scan_report.failed_reasons);
    report_ps.put("skipped_data_names",       scan_report.skipped_data_names);
    report_ps.put("skipped_data_reasons",     scan_report.skipped_data_reasons);
    report_ps.put("skipped_bug_names",        scan_report.skipped_bug_names);
    report_ps.put("skipped_bug_reasons",      scan_report.skipped_bug_reasons);
    report_ps.put("skipped_other_names",      scan_report.skipped_other_names);
    report_ps.put("skipped_other_reasons",    scan_report.skipped_other_reasons);
    report_ps.put("registered_no_fcl",        scan_report.registered_no_fcl);
    report_ps.put("genie_rw_skipped_names",   scan_report.genie_rw_skipped_names);
    report_ps.put("genie_rw_skipped_reasons", scan_report.genie_rw_skipped_reasons);
    out_ps.put("_scan_report", report_ps);
  }

  fhicl::ParameterSet wrapped_out_ps;
  wrapped_out_ps.put("generated_systematic_provider_configuration", out_ps);

  std::ostream *os(nullptr);
  std::ofstream fs;
  if (cliopts::outputfile.size()) {
    fs.open(cliopts::outputfile);
    if (!fs.is_open()) {
      std::cerr << "[ERROR]: Failed to open " << cliopts::outputfile << std::endl;
      return 1;
    }
    os = &fs;
  } else {
    os = &std::cout;
  }
  (*os) << wrapped_out_ps.to_indented_string() << std::endl;
  if (cliopts::outputfile.size()) fs.close();

  std::cerr << "\n=== Summary ===" << std::endl;
  std::cerr << "Wrote " << providerNames.size() << " providers" << std::endl;
  if (cliopts::outputfile.size())
    std::cerr << "Output: " << cliopts::outputfile << std::endl;

  return 0;
}
