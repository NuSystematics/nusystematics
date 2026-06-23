#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"

#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/channel_classification.hh"
#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/genie_metadata.hh"
#include "nusystematics/utility/KinVarUtils.hh"

#include "nusystematics/utility/response_helper.hh"
#include "nusystematics/utility/silence_genie.hh"

#include "fhiclcpp/ParameterSet.h"

#include <map>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Conventions/Units.h"

#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
using namespace systtools;
using namespace nusyst;
using namespace genie;
using namespace genie::rew;

NEW_SYSTTOOLS_EXCEPT(unexpected_number_of_responses);

namespace cliopts {
std::string fclname = "";
std::string genie_input = "";
std::string genie_branch_name = "gmcrec";
std::string outputfile = "";
std::string envvar = "FHICL_FILE_PATH";
std::string fhicl_key = "generated_systematic_provider_configuration";
size_t NMax = std::numeric_limits<size_t>::max();
size_t NSkip = 0;
int n_threads = 1;     // -j / --threads N: spawn N worker processes via fork()
int worker_id = -1;    // -1 == orchestrator / serial; >=0 == fork child index
std::vector<std::string> parameters; // -p substr filter on dial prettyName
#ifndef NO_ART
int lookup_policy = 1;
#endif
} // namespace cliopts

// Cache-resolution constants mirror those in DeclaredDialTestNuSyst so the
// two tools share the same auto-generation behaviour: if no -c is passed,
// fall back to $NUSYST_INVENTORY_FCL and then $nusystematics_ROOT/fcl/
// nusyst_inventory.fcl (then $NUSYST/fcl/..., then /tmp/...), auto-
// generating via GenerateAllDialsConfigNuSyst on first use.
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

// Returns true if `name` matches any substring in cliopts::parameters, or
// if the filter is empty (no filter = keep everything).
bool DialMatchesFilter(const std::string &name) {
  if (cliopts::parameters.empty()) return true;
  for (auto const &sub : cliopts::parameters)
    if (name.find(sub) != std::string::npos) return true;
  return false;
}

struct TweakSummaryTree {
  TFile *f = NULL;
  TTree *t = NULL;
  TTree *m = NULL;

  TweakSummaryTree(std::string const &fname, bool AddWeights) {
    f = new TFile(fname.c_str(), "RECREATE");
    t = new TTree("events", "");
    if(AddWeights){
      m = new TTree("tweak_metadata", "");
    }
    t->SetDirectory(f);
  }
  ~TweakSummaryTree() {
    f->Write();
    f->Close();
    delete f;
  }

  // TH: Add variables for for output weight tree
  int Mode, nucleon_pdg, target_pdg;
  std::string Mode_str;
  double xsec;
  double Emiss, Emiss_preFSI, pmiss, pmiss_preFSI;
  double Emiss_GENIE;
  double q0, Enu_true, q3, Q2, w;
  double plep;
  int nu_pdg;
  double e_nu_GeV;
  int tgt_A;
  int tgt_Z;
  bool is_cc;
  bool is_qe;
  bool is_mec;
  int mec_topology;
  bool is_res;
  int res_channel;
  bool is_dis;
  double EAvail_GeV;
  // Calorimetric reconstructed neutrino energy: EAvail + E_lep. Same
  // definition as PlotSystVariationsNuSyst's EventVars.Ereco_cal.
  double Ereco_cal;
  // FS particle multiplicities (no KE threshold) -- derivable from
  // fsparticles_pdg via Sum$(fsparticles_pdg==2212) etc., but exposed as
  // dedicated scalar branches because they're commonly plotted.
  int nproton, npip, npim, npi0, nneutron;
  std::vector<int> fsi_pdgs;
  std::vector<int> fsi_codes;
  // TKI
  std::vector<double> fsprotons_KE;
  // Highest-|p| proton momentum magnitude -- different from
  // fsprotons_KE[0] in general (sorted by KE, not |p|) though they
  // coincide for non-pathological samples.
  double leading_proton_p;
  // Bjorken scaling variable x = Q^2 / (2 M_N q0). Useful for DIS plots.
  double Bjorken_x;
  // Final-state lepton momentum projected along (parallel) and perpendicular
  // to the incoming-neutrino direction. Both in GeV.
  double plep_L;
  double plep_T;
  double deltaPT, deltaalphaT;
  // Experimenting signal selection here..
  bool IsSignal_ICARUS_1muNp0pi;

  // Full per-particle info, captured at three stages of the interaction so
  // downstream tools can build arbitrary TTreeFormula expressions without
  // recompiling. Each (pdg, px, py, pz, E) tuple is index-aligned across
  // the five vectors of its stage.
  //
  //   fsparticles_*        kIStStableFinalState   final state after FSI
  //   prefsi_particles_*   kIStHadronInTheNucleus pre-FSI hadrons (the
  //                                                primary interaction
  //                                                products before INTRANUKE
  //                                                propagates them)
  //   isparticles_*        kIStInitialState       initial state (incoming
  //                                                neutrino, target nucleus,
  //                                                struck nucleon)
  std::vector<int>    fsparticles_pdg;
  std::vector<double> fsparticles_px;
  std::vector<double> fsparticles_py;
  std::vector<double> fsparticles_pz;
  std::vector<double> fsparticles_E;

  std::vector<int>    prefsi_particles_pdg;
  std::vector<double> prefsi_particles_px;
  std::vector<double> prefsi_particles_py;
  std::vector<double> prefsi_particles_pz;
  std::vector<double> prefsi_particles_E;

  std::vector<int>    isparticles_pdg;
  std::vector<double> isparticles_px;
  std::vector<double> isparticles_py;
  std::vector<double> isparticles_pz;
  std::vector<double> isparticles_E;

  std::vector<int> ntweaks;
  std::vector<std::vector<double>> tweak_branches;
  std::vector<double> paramCVResponses;
  std::map<paramId_t, size_t> tweak_indices;

  TObjString *meta_name;
  int meta_n;
  std::vector<double> meta_tweak_values;

  void AddBranches(response_helper &phh) {
    
    // TH: Add branches for output weights tree
    t->Branch("Mode", &Mode, "Mode/I");
    t->Branch("Mode_str", &Mode_str);
    t->Branch("xsec", &xsec, "xsec/D");
    t->Branch("Emiss", &Emiss, "Emiss/D");
    t->Branch("Emiss_preFSI", &Emiss_preFSI, "Emiss_preFSI/D");
    t->Branch("Emiss_GENIE", &Emiss_GENIE, "Emiss_GENIE/D");
    t->Branch("pmiss", &pmiss, "pmiss/D");
    t->Branch("pmiss_preFSI", &pmiss_preFSI, "pmiss_preFSI/D");
    t->Branch("q0", &q0, "q0/D");
    t->Branch("Q2", &Q2, "Q2/D");
    t->Branch("q3", &q3, "q3/D");
    t->Branch("w", &w, "w/D");
    t->Branch("Enu_true", &Enu_true, "Enu_true/D");
    t->Branch("nu_pdg", &nu_pdg, "nu_pdg/I");
    t->Branch("plep", &plep, "plep/D");
    t->Branch("nucleon_pdg", &nucleon_pdg, "nucleon_pdg/I");
    t->Branch("target_pdg", &target_pdg, "target_pdg/I");
      
    size_t vector_idx = 0;
    t->Branch("tgt_A", &tgt_A, "tgt_A/I");
    t->Branch("tgt_Z", &tgt_Z, "tgt_Z/I");
    t->Branch("is_cc", &is_cc, "is_cc/O");
    t->Branch("is_qe", &is_qe, "is_qe/O");
    t->Branch("is_mec", &is_mec, "is_mec/O");
    t->Branch("mec_topology", &mec_topology, "mec_topology/I");
    t->Branch("is_res", &is_res, "is_res/O");
    t->Branch("res_channel", &res_channel, "res_channel/I");
    t->Branch("is_dis", &is_dis, "is_dis/O");
    t->Branch("EAvail_GeV", &EAvail_GeV, "EAvail_GeV/D");
    t->Branch("Ereco_cal",  &Ereco_cal,  "Ereco_cal/D");
    t->Branch("nproton",    &nproton,    "nproton/I");
    t->Branch("npip",       &npip,       "npip/I");
    t->Branch("npim",       &npim,       "npim/I");
    t->Branch("npi0",       &npi0,       "npi0/I");
    t->Branch("nneutron",   &nneutron,   "nneutron/I");
    t->Branch("fsi_pdgs", "vector<int>", &fsi_pdgs);
    t->Branch("fsi_codes", "vector<int>", &fsi_codes);
    t->Branch("fsprotons_KE", "vector<double>", &fsprotons_KE);
    t->Branch("leading_proton_p", &leading_proton_p, "leading_proton_p/D");
    t->Branch("Bjorken_x",        &Bjorken_x,        "Bjorken_x/D");
    t->Branch("plep_L",           &plep_L,           "plep_L/D");
    t->Branch("plep_T",           &plep_T,           "plep_T/D");

    // Per-FS-particle 4-vectors + PDG (index-aligned across the five branches).
    t->Branch("fsparticles_pdg", "vector<int>",    &fsparticles_pdg);
    t->Branch("fsparticles_px",  "vector<double>", &fsparticles_px);
    t->Branch("fsparticles_py",  "vector<double>", &fsparticles_py);
    t->Branch("fsparticles_pz",  "vector<double>", &fsparticles_pz);
    t->Branch("fsparticles_E",   "vector<double>", &fsparticles_E);

    // Pre-FSI hadrons (kIStHadronInTheNucleus): the primary interaction
    // products before INTRANUKE propagates them.
    t->Branch("prefsi_particles_pdg", "vector<int>",    &prefsi_particles_pdg);
    t->Branch("prefsi_particles_px",  "vector<double>", &prefsi_particles_px);
    t->Branch("prefsi_particles_py",  "vector<double>", &prefsi_particles_py);
    t->Branch("prefsi_particles_pz",  "vector<double>", &prefsi_particles_pz);
    t->Branch("prefsi_particles_E",   "vector<double>", &prefsi_particles_E);

    // Initial state (kIStInitialState): incoming neutrino, target nucleus,
    // struck nucleon.
    t->Branch("isparticles_pdg", "vector<int>",    &isparticles_pdg);
    t->Branch("isparticles_px",  "vector<double>", &isparticles_px);
    t->Branch("isparticles_py",  "vector<double>", &isparticles_py);
    t->Branch("isparticles_pz",  "vector<double>", &isparticles_pz);
    t->Branch("isparticles_E",   "vector<double>", &isparticles_E);
    t->Branch("deltaPT", &deltaPT, "deltaPT/D");
    t->Branch("deltaalphaT", &deltaalphaT, "deltaalphaT/D");
    t->Branch("IsSignal_ICARUS_1muNp0pi", &IsSignal_ICARUS_1muNp0pi, "IsSignal_ICARUS_1muNp0pi/O");


    if(m){

      for (paramId_t pid : phh.GetParameters()) { // Need to size vectors first so
                                                  // that realloc doesn't upset
                                                  // the TBranches
        SystParamHeader const &hdr = phh.GetHeader(pid);
        if (hdr.isResponselessParam) {
          continue;
        }

        // -p filter: skip dials whose prettyName doesn't match any of the
        // user-supplied substrings. The full kitchen-sink is still computed
        // per event (response_helper iterates all providers); we just don't
        // emit tree branches for the filtered-out dials, so the output file
        // stays small. For a compute-speed win, hand-roll a tool config.
        if (!DialMatchesFilter(hdr.prettyName)) {
          continue;
        }

        // Backward-compatibility; isCorrection header used to have empty paramVariations
        if (hdr.isCorrection) {
          ntweaks.emplace_back(1);
        } else {
          ntweaks.emplace_back(hdr.paramVariations.size());
        }
        tweak_branches.emplace_back();
        std::fill_n(std::back_inserter(tweak_branches.back()), ntweaks.back(), 1);
        tweak_indices[pid] = vector_idx;

        if (ntweaks.back() > int(meta_tweak_values.size())) {
          meta_tweak_values.resize(ntweaks.back());
        }
        vector_idx++;
      }
      std::fill_n(std::back_inserter(paramCVResponses), ntweaks.size(), 1);

      meta_name = nullptr;
      m->Branch("name", &meta_name);
      m->Branch("ntweaks", &meta_n, "ntweaks/I");
      m->Branch("tweakvalues", meta_tweak_values.data(),
                "tweakvalues[ntweaks]/D");

      for (paramId_t pid : phh.GetParameters()) {
        SystParamHeader const &hdr = phh.GetHeader(pid);
        if (hdr.isResponselessParam) {
          continue;
        }
        // -p filter: same logic as the first pass above; skip dials that
        // never got into tweak_indices.
        if (!DialMatchesFilter(hdr.prettyName)) {
          continue;
        }
        size_t idx = tweak_indices[pid];

        // Find the IGENIESystProvider_tool(ISystProviderTool) for this pid
        int matched_idx_sp = -1;
        for(int idx_sp=0; idx_sp<phh.GetSystProvider().size(); idx_sp++){
          if(phh.GetSystProvider()[idx_sp]->ParamIsHandled(pid)){
            matched_idx_sp = idx_sp;
          }
        }
        if(matched_idx_sp<0){
          printf("[WeightUpdater::CreateGlobalTree] IGENIESystProvider_tool not found from pid = %d\n", int(pid));
          abort();
        }

        std::string this_full_name = phh.GetSystProvider()[matched_idx_sp]->GetFullyQualifiedName()+"_"+hdr.prettyName;

        std::stringstream ss_ntwk("");
        ss_ntwk << "ntweaks_" << this_full_name;
        t->Branch(ss_ntwk.str().c_str(), &ntweaks[idx],
                  (ss_ntwk.str() + "/I").c_str());

        std::stringstream ss_twkr("");
        ss_twkr << "tweak_responses_" << this_full_name;
        t->Branch(ss_twkr.str().c_str(), tweak_branches[idx].data(),
                  (ss_twkr.str() + "[" + ss_ntwk.str() + "]/D").c_str());

        std::stringstream ss_twkcv("");
        ss_twkcv << "paramCVWeight_" << this_full_name;
        t->Branch(ss_twkcv.str().c_str(), &paramCVResponses[idx],
                  (ss_twkcv.str() + "/D").c_str());

        *meta_name = this_full_name.c_str();
        meta_n = ntweaks[idx];
        // Backward-compatibility; isCorrection header used to have empty paramVariations
        if (hdr.isCorrection) {
          meta_tweak_values[0] = hdr.centralParamValue;
        } else {
          std::copy_n(hdr.paramVariations.begin(), meta_n,
                      meta_tweak_values.begin());
        }

        // In multi-process mode, only worker 0 writes the metadata tree so the
        // hadd-merged output has a single un-duplicated copy.
        if (cliopts::worker_id <= 0) {
          m->Fill();
        }
      }

    } // IF m set

  }

  // Clear weight vectors
  void Clear() {
    std::fill_n(ntweaks.begin(), ntweaks.size(), 0);
    std::fill_n(paramCVResponses.begin(), ntweaks.size(), 1);
  }
  void Add(event_unit_response_t const &eu) {
    for (std::pair<paramId_t, size_t> idx_id : tweak_indices) {
      size_t resp_idx = GetParamContainerIndex(eu, idx_id.first);
      if (resp_idx != systtools::kParamUnhandled<size_t>) {
        ParamResponses const &resp = eu[resp_idx];
        if (tweak_branches[idx_id.second].size() != resp.responses.size()) {
          throw unexpected_number_of_responses()
              << "[ERROR]: Expected " << ntweaks[idx_id.second]
              << " responses from parameter " << resp.pid << ", but found "
              << resp.responses.size();
        }
        ntweaks[idx_id.second] = resp.responses.size();
        std::copy_n(resp.responses.begin(), ntweaks[idx_id.second],
                    tweak_branches[idx_id.second].begin());
      } else {
        ntweaks[idx_id.second] = 7;
        std::fill_n(tweak_branches[idx_id.second].begin(),
                    ntweaks[idx_id.second], 1);
      }
    }
  }
  void Add(event_unit_response_w_cv_t const &eu) {
    for (std::pair<paramId_t, size_t> idx_id : tweak_indices) {
      size_t resp_idx = GetParamContainerIndex(eu, idx_id.first);
      if (resp_idx != systtools::kParamUnhandled<size_t>) {
        VarAndCVResponse const &prcw = eu[resp_idx];
        if (tweak_branches[idx_id.second].size() != prcw.responses.size()) {
          throw unexpected_number_of_responses()
              << "[ERROR]: Expected " << ntweaks[idx_id.second]
              << " responses from parameter " << prcw.pid << ", but found "
              << prcw.responses.size();
        }
        ntweaks[idx_id.second] = prcw.responses.size();
        std::copy_n(prcw.responses.begin(), ntweaks[idx_id.second],
                    tweak_branches[idx_id.second].begin());
        paramCVResponses[idx_id.second] = prcw.CV_response;

      } else {
        ntweaks[idx_id.second] = 7;
        std::fill_n(tweak_branches[idx_id.second].begin(),
                    ntweaks[idx_id.second], 1);
        paramCVResponses[idx_id.second] = 1;
      }
    }
  }

  void Fill() { t->Fill(); }
};

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-c <config.fcl>  : fhicl file to read (parameter-headers\n"
               "\t                   format, i.e. the output of\n"
               "\t                   `nusyst config`). Optional: if omitted,\n"
               "\t                   resolves to $NUSYST_INVENTORY_FCL\n"
               "\t                   then $NUSYST/fcl/nusyst_inventory.fcl,\n"
               "\t                   auto-generating via `nusyst config`\n"
               "\t                   if absent.\n"
               "\t-p <par1,par2,...>: Filter to dials whose prettyName\n"
               "\t                   contains any of the comma-separated\n"
               "\t                   substrings. Providers whose dials all\n"
               "\t                   miss the filter are not constructed and\n"
               "\t                   not evaluated (compute-side skip);\n"
               "\t                   non-matching dials within a kept\n"
               "\t                   provider are still computed but not\n"
               "\t                   written to the output tree.\n"
               "\t-k <list key>    : fhicl key to look for parameter headers,\n"
               "\t                   "
               "\"generated_systematic_provider_configuration\"\n"
               "\t                   by default.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-b <NtpMCEventRecord branch name>   : Name of the NtpMCEventRecord branch (default:gmcrec)\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
               "\t-o <out.root>    : File to write validation canvases to.\n"
               "\t-j|--threads N   : Spawn N worker subprocesses via fork()\n"
               "\t                   to process disjoint event ranges in\n"
               "\t                   parallel, then hadd the parts into the\n"
               "\t                   final output. Default 1 = single process.\n"
               "\t                   (Implemented as multi-process, not\n"
               "\t                   std::thread, because GENIE's reweight\n"
               "\t                   stack carries non-thread-safe globals.)\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-p") {
      std::string tok;
      std::istringstream ss(argv[++opt]);
      while (std::getline(ss, tok, ','))
        if (!tok.empty()) cliopts::parameters.push_back(tok);
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::genie_input = argv[++opt];
    } else if (std::string(argv[opt]) == "-b") {
      cliopts::genie_branch_name = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-j" ||
               std::string(argv[opt]) == "--threads") {
      cliopts::n_threads = str2T<int>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

typedef IGENIESystProvider_tool SystProv;

fhicl::ParameterSet ReadParameterSet(char const *[]) {
  // TODO
  std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
  return fhicl::ParameterSet::make(cliopts::fclname, *fm);
}

// Re-parse the parameter-headers fhicl to extract per-provider
// `applies_to_channels` patterns (written by GenerateAllDialsConfigNuSyst's
// per-bucket emission). Returns an empty map if no provider declares the key,
// in which case no channel-aware skipping happens and behaviour matches the
// pre-optimisation baseline.
std::map<std::string, std::vector<std::string>>
LoadAppliesToChannelsMap(const std::string &fclname,
                         const std::string &top_key) {
  std::map<std::string, std::vector<std::string>> out;
  try {
    std::unique_ptr<cet::filepath_maker> fm =
        std::make_unique<cet::filepath_lookup_nonabsolute>("FHICL_FILE_PATH");
    fhicl::ParameterSet raw = fhicl::ParameterSet::make(fclname, *fm);
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
  } catch (std::exception &e) {
    std::cerr << "[WARN]: Failed to load applies_to_channels map from "
              << fclname << ": " << e.what()
              << ". Falling back to evaluating every provider on every event."
              << std::endl;
    return {};
  }
  return out;
}

// Build a "trivial" response for a provider -- one VarAndCVResponse per
// non-responseless dial, with CV=1 and every variation weight =1. Used as a
// drop-in for `GetEventVariationAndCVResponse` when the provider is skipped
// for the current event's channel.
event_unit_response_w_cv_t
TrivialResponseFor(IGENIESystProvider_tool const &provider) {
  event_unit_response_w_cv_t out;
  for (auto const &hdr : provider.GetSystMetaData()) {
    if (hdr.isResponselessParam) continue;
    VarAndCVResponse r;
    r.pid = hdr.systParamId;
    r.CV_response = 1.0;
    size_t nvar = hdr.isCorrection ? 1 : hdr.paramVariations.size();
    r.responses.assign(nvar, 1.0);
    out.push_back(std::move(r));
  }
  return out;
}

// Forward declarations for the multi-process dispatcher (see end of file).
int RunSerial();
int DispatchWorkers();

int RunSerial() {
  // Cache fallback: a `-p` filter implicitly says "I want some reweights,
  // configure them from the cached kitchen sink". Resolve fclname against
  // $NUSYST_INVENTORY_FCL / /tmp/nusyst_inventory.fcl and auto-generate via
  // `nusyst config` if missing, mirroring DeclaredDialTestNuSyst's logic.
  // Plain `nusyst tweaks -i ghep.root -o out.root` (no -c, no -p) keeps the
  // pre-existing "no reweights" behaviour for backwards compatibility.
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

  bool RunNuSyst = true;
  if (!cliopts::fclname.size()) {
    RunNuSyst = false;
    std::cout << "-c is not given, running without evaluating reweights" << std::endl;
  }
  if (!cliopts::genie_input.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option. "
              << "(Run with -? for usage.)" << std::endl;
    return 1;
  }

  nusyst::quiet::SetGlobalQuiet();

  // Silence GENIE's banner + NOTICE chatter unconditionally. The event loop
  // touches GENIE per event (NeutReactionCode, ProcInfo, ...) regardless of
  // whether we're running reweights, and the very first GENIE call triggers
  // Messenger::Instance() which prints the multi-screen banner. Wrap that
  // first touch in a StdoutSink and raise NOTICE+ to WARN so the per-event
  // GHepUtils notices stop printing.
  {
    nusyst::quiet::StdoutSink _quiet;
    genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
        "Messenger_whisper.xml");
  }

  response_helper phh;
  std::map<std::string, std::vector<std::string>> applies_to_channels;
  size_t n_filtered_providers = 0;
  if(RunNuSyst){
    // Load the parameter-headers fhicl as a ParameterSet so a -p filter can
    // drop providers before they're constructed. Without the filter, the
    // effect is identical to the old `response_helper(fclname)` one-shot.
    fhicl::ParameterSet raw_ps;
    {
      std::unique_ptr<cet::filepath_maker> fm =
          std::make_unique<cet::filepath_lookup_nonabsolute>("FHICL_FILE_PATH");
      raw_ps = fhicl::ParameterSet::make(cliopts::fclname, *fm);
    }
    fhicl::ParameterSet gen_ps =
        raw_ps.get<fhicl::ParameterSet>(cliopts::fhicl_key);

    // Provider-level -p filter: drop any provider whose dials all miss the
    // filter. These are not instantiated and not evaluated per event.
    // Providers with a partial match are kept whole (still compute their
    // full dial set); the per-dial output filter below trims their branches
    // from the tree.
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
            if (DialMatchesFilter(pretty)) { any_match = true; break; }
          } catch (...) {}
        }
        if (any_match) {
          kept.push_back(pname);
        } else {
          gen_ps.erase(pname);
        }
      }
      gen_ps.put_or_replace<std::vector<std::string>>("syst_providers", kept);
      std::cerr << "[INFO]: -p filter kept " << kept.size() << " of "
                << provider_names.size() << " providers ("
                << (provider_names.size() - kept.size())
                << " skipped -- neither constructed nor evaluated).\n";
    }

    {
      nusyst::quiet::StdoutSink _quiet;
      genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
          "Messenger_whisper.xml");
      phh.LoadProvidersAndHeaders(gen_ps);
    }

    // Re-derive applies_to_channels from the (possibly filtered) PS instead
    // of re-reading the file -- keeps the two views in sync.
    auto provider_names =
        gen_ps.get<std::vector<std::string>>("syst_providers",
                                              std::vector<std::string>{});
    for (auto const &pname : provider_names) {
      try {
        fhicl::ParameterSet prov = gen_ps.get<fhicl::ParameterSet>(pname);
        fhicl::ParameterSet topts =
            prov.get<fhicl::ParameterSet>("tool_options",
                                           fhicl::ParameterSet{});
        auto patterns = topts.get<std::vector<std::string>>(
            "applies_to_channels", std::vector<std::string>{});
        if (!patterns.empty()) {
          applies_to_channels[pname] = std::move(patterns);
        }
      } catch (...) {}
    }
    for (auto &sp : phh.GetSystProvider()) {
      if (applies_to_channels.count(sp->GetFullyQualifiedName())) {
        ++n_filtered_providers;
      }
    }
    std::cout << "[INFO]: " << n_filtered_providers << " of "
              << phh.GetSystProvider().size()
              << " providers declare applies_to_channels; "
              << (phh.GetSystProvider().size() - n_filtered_providers)
              << " will be evaluated on every event." << std::endl;
  }

  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(cliopts::genie_input.c_str())) {
    std::cout << "[ERROR]: Failed to find any TTrees named "
              << std::quoted("gtree") << ", from TChain::Add descriptor: "
              << std::quoted(cliopts::genie_input) << "." << std::endl;
    return 3;
  }

  // Surface the GENIE version/tune that generated this sample so the user can
  // diagnose version-skew issues (e.g. samples from GENIE 3.04 carry
  // QE-event phase-space labels that GENIE 3.06's reweight can't transform).
  // Only the orchestrator (worker_id < 0) reports -- workers stay quiet.
  if (cliopts::worker_id < 0) {
    nusyst::metadata::ReportSampleInfo(cliopts::genie_input);
  }

  size_t NEvs = gevs->GetEntries();

  if (!NEvs) {
    std::cout << "[ERROR]: Input TChain contained no entries." << std::endl;
    return 4;
  }

  if( cliopts::NSkip >= NEvs ){
    printf("[ERROR]: NSkip is larger than NEvs; (NSkip, NEvs) = (%ld, %ld)\n", cliopts::NSkip, NEvs);
    return 5;
  }

  genie::NtpMCEventRecord *GenieNtpl = nullptr;

  if (gevs->SetBranchAddress(cliopts::genie_branch_name.c_str(), &GenieNtpl) != TTree::kMatch) {
    std::cout << "[ERROR]: Failed to set branch address on ghep tree."
              << std::endl;
    return 6;
  }

  TweakSummaryTree tst(cliopts::outputfile.c_str(), RunNuSyst);
  tst.AddBranches(phh);

  size_t NToRead = std::min(NEvs, cliopts::NMax);
  // Progress every 5% of events, capped at 1000 so big runs don't go
  // silent for a minute+ between prints (FSI dials are ~20 ms/event and
  // a NToRead/20 stride at 100k = 5000 events = ~90 s between updates).
  size_t NToShout = std::min<size_t>(NToRead / 20, 1000);
  NToShout = NToShout ? NToShout : 1;

  // Silence GENIE Reweight engines' direct-cout chatter (FSI dials are the
  // worst offenders) for the whole event loop in one shot. /dev/null is
  // opened once here rather than re-opened per event -- the previous
  // per-event StdoutSink scope made FrAbs_pi-style runs ~100x slower from
  // the syscall storm. Progress prints already go to cerr, which the
  // StdoutSink doesn't touch.
  nusyst::quiet::StdoutSink _quiet_event_loop;

  for (size_t ev_it = cliopts::NSkip; ev_it < NToRead; ++ev_it) {
    // Start-of-event: reset tweak containers so later Add() fills cleanly.
    // (This does NOT touch physics scalars like w, q0, Q2.)
    tst.Clear();

    gevs->GetEntry(ev_it);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;

    tst.xsec = GenieGHep.XSec()/genie::units::cm2;

    genie::GHepParticle *FSLep = GenieGHep.FinalStatePrimaryLepton();
    genie::GHepParticle *ISLep = GenieGHep.Probe();
    genie::GHepParticle *nucleon = GenieGHep.HitNucleon();
    
    TLorentzVector FSLepP4 = *FSLep->P4();
    TLorentzVector ISLepP4 = *ISLep->P4();
    TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

    tst.Mode = genie::utils::ghep::NeutReactionCode(&GenieGHep);
    tst.Mode_str = GenieGHep.Summary()->AsString();
    tst.Emiss = GetEmiss(GenieGHep, false);
    tst.Emiss_preFSI = GetEmiss(GenieGHep, true);
    tst.pmiss = GetPmiss(GenieGHep, false);
    tst.pmiss_preFSI = GetPmiss(GenieGHep, true);

    if (GenieGHep.HitNucleon() == NULL){
      tst.Emiss_GENIE = -999;
    }
    else {
      tst.Emiss_GENIE = GenieGHep.HitNucleon()->RemovalEnergy();
    }

    tst.q0 = emTransfer.E();
    tst.Q2 = -emTransfer.Mag2();
    tst.q3 = emTransfer.Vect().Mag();
    
    // Compute W from (q0, Q2) and nucleon mass
    // W^2 = M_N^2 + 2 M_N q0 - Q^2  (q0 = energy transfer)
    // Use consistent nucleon mass value from nuisance
    double MN = 0.93827203; // GeV
    double W2 = MN * MN + 2.0 * MN * tst.q0 - tst.Q2;
    tst.w = (W2 > 0.0) ? std::sqrt(W2) : -1.0;
    
    tst.Enu_true = ISLepP4.E();
    tst.nu_pdg = ISLep->Pdg();

    tst.plep = FSLepP4.Vect().Mag();
    if (nucleon == NULL) {tst.nucleon_pdg = -999;}
    else{tst.nucleon_pdg = nucleon->Pdg();}

    if(GenieGHep.TargetNucleus()){
      tst.target_pdg = GenieGHep.TargetNucleus()->Pdg();
      tst.tgt_A = GenieGHep.TargetNucleus()->A();
      tst.tgt_Z = GenieGHep.TargetNucleus()->Z();
    }
    else{
      tst.target_pdg = 1000010010;
      tst.tgt_A = 1;
      tst.tgt_Z = 1;
    }

    tst.is_cc = GenieGHep.Summary()->ProcInfo().IsWeakCC();
    tst.is_qe = GenieGHep.Summary()->ProcInfo().IsQuasiElastic();
    tst.is_mec = GenieGHep.Summary()->ProcInfo().IsMEC();
    tst.mec_topology = -1;
    if (tst.is_mec) {
      tst.mec_topology = e2i(GetQELikeTarget(GenieGHep));
    }
    tst.is_res = GenieGHep.Summary()->ProcInfo().IsResonant();
    tst.res_channel = 0;
    if (tst.is_res) {
      tst.res_channel = SPPChannelFromGHep(GenieGHep);
    }
    tst.is_dis = GenieGHep.Summary()->ProcInfo().IsDeepInelastic();

    tst.EAvail_GeV = GetErecoil_MINERvA_LowRecoil(GenieGHep);
    // Calorimetric reconstructed neutrino energy = EAvail + E_lep. Matches
    // PlotSystVariationsNuSyst's EventVars.Ereco_cal definition.
    tst.Ereco_cal = tst.EAvail_GeV + FSLepP4.E();

    // Bjorken x = Q^2 / (2 M_N q0). Guard q0 -> 0 with a tiny floor.
    {
      double safe_q0 = (tst.q0 > 1e-9) ? tst.q0 : 1e-9;
      tst.Bjorken_x = tst.Q2 / (2.0 * MN * safe_q0);
    }

    // pL = lepton momentum along the IS-neutrino direction,
    // pT = lepton momentum perpendicular to it. Falls back to (plep, 0) if
    // the neutrino is along +z (no preferred direction in some test files).
    {
      TVector3 nuDir = ISLepP4.Vect();
      double nu_mag = nuDir.Mag();
      if (nu_mag > 1e-9) {
        nuDir *= 1.0 / nu_mag;
        TVector3 lp = FSLepP4.Vect();
        double pL = lp.Dot(nuDir);
        double pmag = lp.Mag();
        double pT2 = std::max(0.0, pmag * pmag - pL * pL);
        tst.plep_L = pL;
        tst.plep_T = std::sqrt(pT2);
      } else {
        tst.plep_L = tst.plep;
        tst.plep_T = 0.0;
      }
    }

    // loop over particles
    int ip=-1;
    GHepParticle * p = 0;
    TIter event_iter(&GenieGHep);

    std::vector<int> fsi_pdgs;
    std::vector<int> fsi_codes;

    // Reset per-FS-particle vectors + multiplicity counters for the new event.
    tst.fsparticles_pdg.clear();
    tst.fsparticles_px.clear();
    tst.fsparticles_py.clear();
    tst.fsparticles_pz.clear();
    tst.fsparticles_E.clear();
    tst.prefsi_particles_pdg.clear();
    tst.prefsi_particles_px.clear();
    tst.prefsi_particles_py.clear();
    tst.prefsi_particles_pz.clear();
    tst.prefsi_particles_E.clear();
    tst.isparticles_pdg.clear();
    tst.isparticles_px.clear();
    tst.isparticles_py.clear();
    tst.isparticles_pz.clear();
    tst.isparticles_E.clear();
    tst.nproton = 0;
    tst.npip = 0;
    tst.npim = 0;
    tst.npi0 = 0;
    tst.nneutron = 0;

    // Particle loop
    std::vector<GHepParticle *> protons;

    // ICARUS 1muNp0pi signal definition
    unsigned int nMu_1muNp0pi(0), nP_1muNp0pi(0), nPi_1muNp0pi(0);
    unsigned int nPhoton_1muNp0pi(0), nMesons_1muNp0pi(0), nBaryonsAndPi0_1muNp0pi(0);
    double maxMomentumP_1muNp0pi = -999.;
    bool passProtonPCut_1muNp0pi = false;

    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      ip++;

      // Skip particles not rescattered by the actual hadron transport code
      int  pdgc       = p->Pdg();
      bool is_pion    = pdg::IsPion   (pdgc);
      bool is_nucleon = pdg::IsNucleon(pdgc);
      bool is_proton = pdg::IsProton(pdgc);
      bool is_kaon = pdg::IsKaon( pdgc );

      GHepStatus_t ist  = p->Status();

      // For FSI study
      if(ist==kIStHadronInTheNucleus){
        // Kaon FSIs can't currently be reweighted. Just update (A, Z) based on
        // the particle's daughters and move on.
        if ( is_pion || is_nucleon ){
          int fsi_code = p->RescatterCode();
          fsi_pdgs.push_back(pdgc);
          fsi_codes.push_back(fsi_code);
        }
        // Record the pre-FSI hadron 4-vector regardless of FSI-code support.
        // These are the primary interaction products before INTRANUKE
        // propagates them.
        TLorentzVector const *p4 = p->P4();
        tst.prefsi_particles_pdg.push_back(pdgc);
        tst.prefsi_particles_px.push_back(p4->Px());
        tst.prefsi_particles_py.push_back(p4->Py());
        tst.prefsi_particles_pz.push_back(p4->Pz());
        tst.prefsi_particles_E.push_back(p4->E());
      }

      // Initial state: probe + target nucleus + struck nucleon (whichever
      // entries GENIE flags as IS for this interaction).
      if (ist == genie::kIStInitialState) {
        TLorentzVector const *p4 = p->P4();
        tst.isparticles_pdg.push_back(pdgc);
        tst.isparticles_px.push_back(p4->Px());
        tst.isparticles_py.push_back(p4->Py());
        tst.isparticles_pz.push_back(p4->Pz());
        tst.isparticles_E.push_back(p4->E());
      }

      // Stable final state particle
      if(ist==genie::kIStStableFinalState){
        // Record this particle's full 4-vector so downstream TTreeFormula
        // expressions in `nusyst plots` can compute arbitrary quantities
        // (e.g. leading-momentum proton energy) without recompiling.
        TLorentzVector const *p4 = p->P4();
        tst.fsparticles_pdg.push_back(pdgc);
        tst.fsparticles_px.push_back(p4->Px());
        tst.fsparticles_py.push_back(p4->Py());
        tst.fsparticles_pz.push_back(p4->Pz());
        tst.fsparticles_E.push_back(p4->E());

        // Multiplicity counters (no KE threshold; mirrors EventVars in
        // PlotSystVariationsNuSyst's GHEP-mode fill).
        if      (pdgc ==  2212) tst.nproton++;
        else if (pdgc ==   211) tst.npip++;
        else if (pdgc ==  -211) tst.npim++;
        else if (pdgc ==   111) tst.npi0++;
        else if (pdgc ==  2112) tst.nneutron++;

        // All FS protons for generic purpose
        if(is_proton){
          protons.push_back(p);
        }

        // ICARUS 1muNp0pi
        double momentum = p->P4()->Vect().Mag();

        bool PassMuonPCut = (momentum > 0.226);
        if ( abs(pdgc) == 13 ) {
          if (PassMuonPCut) nMu_1muNp0pi+=1;
        }

        if ( abs(pdgc) == 2212 ) {
          nP_1muNp0pi+=1;
          if ( momentum > maxMomentumP_1muNp0pi ) {
            maxMomentumP_1muNp0pi = momentum;
            passProtonPCut_1muNp0pi = (momentum > 0.4 && momentum < 1.);
          }
        }

        if ( abs(pdgc) == 111 || abs(pdgc) == 211 ) nPi_1muNp0pi+=1;
        // CHECK A SIMILAR DEFINITION AS MINERVA FOR EXTRA REJECTION OF UNWANTED THINGS IN SIGNAL DEFN.
        if ( abs(pdgc) == 22 && p->E() > 0.01 ) nPhoton_1muNp0pi+=1;
        else if ( abs(pdgc) == 211 || abs(pdgc) == 321 || abs(pdgc) == 323 ||
                  pdgc == 111 || pdgc == 130 || pdgc == 310 || pdgc == 311 ||
                  pdgc == 313 || abs(pdgc) == 221 || abs(pdgc) == 331 ) nMesons_1muNp0pi+=1;
        else if ( pdgc == 3112 || pdgc == 3122 || pdgc == 3212 || pdgc == 3222 ||
                  pdgc == 4112 || pdgc == 4122 || pdgc == 4212 || pdgc == 4222 ||
                  pdgc == 411 || pdgc == 421 || pdgc == 111 ) nBaryonsAndPi0_1muNp0pi+=1;

      }


    } // END particle loop
    tst.fsi_pdgs = fsi_pdgs;
    tst.fsi_codes = fsi_codes;

    // Final state protons
    // - Sort protons in descending order of KE
    std::sort(protons.begin(), protons.end(), 
              [](GHepParticle* a, GHepParticle* b) {
                  return a->KinE() > b->KinE();
              });

    tst.fsprotons_KE.clear();
    tst.leading_proton_p = -999;
    double max_p_proton = -1;
    for(const auto& proton: protons){
      double this_KE = proton->KinE();
      tst.fsprotons_KE.push_back(this_KE);
      double pmag = proton->P4()->Vect().Mag();
      if (pmag > max_p_proton) { max_p_proton = pmag; tst.leading_proton_p = pmag; }
    }
    // Calculate TKI
    double deltaPT = -999.;
    double deltaalphaT = -999.;
    if(protons.size()>0){
      deltaPT = CalcTKI_deltaPT(
        FSLepP4.Vect(),
        protons[0]->P4()->Vect(),
        ISLepP4.Vect()
      );
      deltaalphaT = CalcTKI_deltaalphaT(
        FSLepP4.Vect(),
        protons[0]->P4()->Vect(),
        ISLepP4.Vect()
      );
    }

    tst.deltaPT = deltaPT;
    tst.deltaalphaT = deltaalphaT;

    // ICARUS 1muNp0pi
    tst.IsSignal_ICARUS_1muNp0pi = nMu_1muNp0pi==1 && 
                                   nP_1muNp0pi>0 && passProtonPCut_1muNp0pi &&
                                   nPi_1muNp0pi==0 &&
                                   nPhoton_1muNp0pi==0 &&
                                   nMesons_1muNp0pi==0 &&
                                   nBaryonsAndPi0_1muNp0pi==0;

    // In multi-process mode only worker 0 prints progress (others would
    // interleave illegibly); single-process behaves as before.
    // Progress goes to cerr so it survives the per-event StdoutSink that
    // silences direct cout writes from the GENIE Reweight engines (FSI
    // dials in particular dump per-event diagnostics through plain cout,
    // bypassing the GENIE Messenger).
    if (cliopts::worker_id <= 0 && !(ev_it % NToShout)) {
      std::cerr << (ev_it ? "\r" : "") << "Event #" << ev_it << "/" << NToRead
                << ", Interaction: " << GenieGHep.Summary()->AsString()
                << std::flush;
    }

    

    // Calcuate weights -- per-provider loop so we can skip providers whose
    // `applies_to_channels` patterns don't match this event's topology.
    // Skipped providers contribute a trivial response (CV=1, all variations=1)
    // for each of their dials, preserving the event's flat-tree shape.
    if(RunNuSyst){
      event_unit_response_w_cv_t resp;
      std::string ch = (n_filtered_providers > 0)
                           ? nusyst::channel::MakeChannelKey(GenieGHep)
                           : std::string{};
      auto &providers = phh.GetSystProvider();
      // cout is already redirected to /dev/null by the loop-scope
      // StdoutSink above, so direct cout writes from inside FSI / other
      // reweight engines are silenced without per-event syscall cost.
      for (auto &sp : providers) {
        auto it = applies_to_channels.find(sp->GetFullyQualifiedName());
        bool applies = (it == applies_to_channels.end())
                           ? true
                           : nusyst::channel::MatchesAny(ch, it->second);
        event_unit_response_w_cv_t prov_resp =
            applies ? sp->GetEventVariationAndCVResponse(GenieGHep)
                    : TrivialResponseFor(*sp);
        for (auto &er : prov_resp) resp.push_back(std::move(er));
      }
      tst.Add(resp);
    }

    tst.Fill();
    
    // TH: Very important to clear this object to avoid memory issues!
    GenieNtpl->Clear();

  }
  std::cout << std::endl;
  return 0;
}

// Spawn cliopts::n_threads worker subprocesses via fork(), each handling a
// disjoint slice of the input event range. Workers write to <out>.partXX
// files; the orchestrator hadd's them into <out> on success and deletes the
// parts. Returns 0 on success, non-zero if any worker failed or hadd failed.
//
// fork() (rather than std::thread / OpenMP) is mandatory: GENIE's reweight
// machinery keeps non-thread-safe global state (genie::Messenger, the tune
// singletons, RandomGen, GReWeight engine caches). Each worker gets its own
// address space and initialises GENIE independently.
int DispatchWorkers() {
  if (cliopts::outputfile.empty()) {
    std::cerr << "[ERROR]: --threads N>1 requires -o <out.root>." << std::endl;
    return 1;
  }

  // Peek at the input to compute event count without touching GENIE. The
  // TChain is destroyed (closing any opened files) before fork().
  size_t NEvs = 0;
  {
    TChain peek("gtree");
    if (!peek.Add(cliopts::genie_input.c_str())) {
      std::cerr << "[ERROR]: Could not open " << cliopts::genie_input << std::endl;
      return 1;
    }
    NEvs = peek.GetEntries();
  }
  if (!NEvs) {
    std::cerr << "[ERROR]: Input chain is empty." << std::endl;
    return 1;
  }

  size_t user_start = cliopts::NSkip;
  size_t user_end = std::min(cliopts::NMax, NEvs);
  if (user_end <= user_start) {
    std::cerr << "[ERROR]: -s/-N leaves zero events to process." << std::endl;
    return 1;
  }
  size_t total = user_end - user_start;
  int N = std::min<int>(cliopts::n_threads, static_cast<int>(total));
  size_t per = (total + N - 1) / N;

  std::vector<pid_t> pids;
  std::vector<std::string> partfiles;
  std::string base = cliopts::outputfile;

  std::cout << "[INFO]: Dispatching " << N << " workers across "
            << total << " events (~" << per << " per worker)." << std::endl;

  for (int k = 0; k < N; ++k) {
    size_t start = user_start + static_cast<size_t>(k) * per;
    size_t end = std::min(start + per, user_end);
    if (start >= end) break;
    std::string part = base + ".part" + std::to_string(k);

    pid_t pid = fork();
    if (pid == 0) {
      // Child: configure as worker k and run the serial path.
      cliopts::worker_id = k;
      cliopts::n_threads = 1;
      cliopts::NSkip = start;
      cliopts::NMax = end;
      cliopts::outputfile = part;
      int rc = RunSerial();
      // _exit so we don't run global destructors twice.
      _exit(rc == 0 ? 0 : 1);
    } else if (pid < 0) {
      std::perror("fork");
      return 1;
    }
    pids.push_back(pid);
    partfiles.push_back(part);
    std::cout << "  worker " << k << " (pid " << pid << "): events ["
              << start << ", " << end << ") -> " << part << std::endl;
  }

  int orchestrator_rc = 0;
  for (size_t i = 0; i < pids.size(); ++i) {
    int status = 0;
    waitpid(pids[i], &status, 0);
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
      std::cerr << "[ERROR]: worker " << i << " (pid " << pids[i]
                << ") exited abnormally (status " << status << ")." << std::endl;
      orchestrator_rc = 2;
    }
  }
  if (orchestrator_rc != 0) {
    std::cerr << "[ERROR]: At least one worker failed; not running hadd. "
              << "Part files left in place for inspection:" << std::endl;
    for (auto &p : partfiles) std::cerr << "    " << p << std::endl;
    return orchestrator_rc;
  }

  std::string cmd = "hadd -f " + base;
  for (auto &p : partfiles) cmd += " " + p;
  std::cout << "[INFO]: Merging part files: " << cmd << std::endl;
  int hadd_rc = std::system(cmd.c_str());
  if (hadd_rc != 0) {
    std::cerr << "[ERROR]: hadd failed (exit code " << hadd_rc
              << "); leaving part files in place." << std::endl;
    return 3;
  }

  for (auto &p : partfiles) std::remove(p.c_str());
  std::cout << "[INFO]: Done. Output: " << base << std::endl;
  return 0;
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);
  if (cliopts::n_threads > 1 && cliopts::worker_id < 0) {
    return DispatchWorkers();
  }
  return RunSerial();
}
