#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/KinVarUtils.hh"

#ifdef NO_ART
#include "nusystematics/artless/response_helper.hh"
#endif

#include "fhiclcpp/ParameterSet.h"

#ifndef NO_ART
#include "cetlib/filepath_maker.h"
#else
#include "fhiclcpp/make_ParameterSet.h"
#endif

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace systtools;
using namespace nusyst;

NEW_SYSTTOOLS_EXCEPT(unexpected_number_of_responses);

struct TweakSummaryTree {
  TFile *f;
  TTree *t;
  TTree *m;

  TweakSummaryTree(std::string const &fname) {
    f = new TFile(fname.c_str(), "RECREATE");
    t = new TTree("events", "");
    m = new TTree("tweak_metadata", "");
    t->SetDirectory(f);
  }
  ~TweakSummaryTree() {
    f->Write();
    f->Close();
    delete f;
  }

  // TH: Add variables for for output weight tree
  int Mode;
  double Emiss, Emiss_preFSI, pmiss, pmiss_preFSI;
  std::vector<int> ntweaks;
  std::vector<std::vector<double>> tweak_branches;
  std::vector<double> paramCVResponses;
  std::map<paramId_t, size_t> tweak_indices;

  TObjString *meta_name;
  int meta_n;
  std::vector<double> meta_tweak_values;

  void AddBranches(ParamHeaderHelper const &phh) {
    
	// TH: Add branches for output weights tree
  t->Branch("Mode", &Mode, "Mode/I");
  t->Branch("Emiss", &Emiss, "Emiss/D");
  t->Branch("Emiss_preFSI", &Emiss_preFSI, "Emiss_preFSI/D");
  t->Branch("pmiss", &pmiss, "pmiss/D");
  t->Branch("pmiss_preFSI", &pmiss_preFSI, "pmiss_preFSI/D");
    
	size_t vector_idx = 0;
    for (paramId_t pid : phh.GetParameters()) { // Need to size vectors first so
                                                // that realloc doesn't upset
                                                // the TBranches
      SystParamHeader const &hdr = phh.GetHeader(pid);
      if (hdr.isResponselessParam) {
        continue;
      }

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
      size_t idx = tweak_indices[pid];

      std::stringstream ss_ntwk("");
      ss_ntwk << "ntweaks_" << hdr.prettyName;
      t->Branch(ss_ntwk.str().c_str(), &ntweaks[idx],
                (ss_ntwk.str() + "/I").c_str());

      std::stringstream ss_twkr("");
      ss_twkr << "tweak_responses_" << hdr.prettyName;
      t->Branch(ss_twkr.str().c_str(), tweak_branches[idx].data(),
                (ss_twkr.str() + "[" + ss_ntwk.str() + "]/D").c_str());

      std::stringstream ss_twkcv("");
      ss_twkcv << "paramCVWeight_" << hdr.prettyName;
      t->Branch(ss_twkcv.str().c_str(), &paramCVResponses[idx],
                (ss_twkcv.str() + "/D").c_str());

      *meta_name = hdr.prettyName.c_str();
      meta_n = ntweaks[idx];
      std::copy_n(hdr.paramVariations.begin(), meta_n,
                  meta_tweak_values.begin());

      m->Fill();
    }
  }

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

namespace cliopts {
std::string fclname = "";
std::string genie_input = "";
std::string outputfile = "";
std::string envvar = "FHICL_FILE_PATH";
std::string fhicl_key = "generated_systematic_provider_configuration";
size_t NMax = std::numeric_limits<size_t>::max();
#ifndef NO_ART
int lookup_policy = 1;
#endif
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
#ifndef NO_ART
               "\t-l <policy_id>   : FHICL_FILE_PATH lookup policy:\n"
               "\t                    0 : cet::filepath_maker\n"
               "\t                   {1}: cet::filepath_lookup\n"
               "\t                    2 : cet::filepath_lookup_nonabsolute\n"
               "\t                    3 : cet::filepath_lookup_after1\n"
               "\t-p <envvar name> : Environment variable to use when searching"
               " for fhicl. \n"
               "\t                   FHICL_FILE_PATH by default.\n"
#endif
               "\t-c <config.fcl>  : fhicl file to read.\n"
               "\t-k <list key>    : fhicl key to look for parameter headers,\n"
               "\t                   "
               "\"generated_systematic_provider_configuration\"\n"
               "\t                   by default.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-o <out.root>    : File to write validation canvases to.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
#ifndef NO_ART
    } else if (std::string(argv[opt]) == "-l") {
      cliopts::lookup_policy = systtools::str2T<int>(argv[++opt]);
      if (cliopts::lookup_policy > 3 || cliopts::lookup_policy < 0) {
        std::cout << "[ERROR]: -l expected to be passed an integer between 0 "
                     "and 3."
                  << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-p") {
      cliopts::envvar = argv[++opt];
      char const *ev = getenv(cliopts::envvar.c_str());
      if (!ev) {
        std::cout << "[ERROR]: Could not read environment variable:"
                  << std::quoted(cliopts::envvar)
                  << ". Please supply a variable containing a valid path list."
                  << std::endl;
        SayUsage(argv);
        exit(1);
      }
#endif
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::genie_input = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::outputfile = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

#ifdef NO_ART
typedef IGENIESystProvider_tool SystProv;
#else
typedef ISystProviderTool SystProv;
#endif

fhicl::ParameterSet ReadParameterSet(char const *[]) {

#ifndef NO_ART
  char const *ev = nullptr;
  if (cliopts::lookup_policy != 0) {
    ev = getenv(cliopts::envvar.c_str());
    if (!ev) {
      std::cout << "[ERROR]: Could not read environment variable:\""
                << cliopts::envvar
                << "\". Please supply a variable containing a valid path list "
                   "via the -p command line option."
                << std::endl;
      SayUsage(argv);
      exit(1);
    }
  }
#endif

  fhicl::ParameterSet ps;
#ifndef NO_ART
  std::unique_ptr<cet::filepath_maker> fm(nullptr);

  switch (cliopts::lookup_policy) {
  case 0: {
    fm = std::make_unique<cet::filepath_maker>();
    break;
  }
  case 1: {
    fm = std::make_unique<cet::filepath_lookup>(ev);
    break;
  }
  case 2: {
    fm = std::make_unique<cet::filepath_lookup_nonabsolute>(ev);

    break;
  }
  case 3: {
    fm = std::make_unique<cet::filepath_lookup_after1>(ev);
    break;
  }
  default: {}
  }
  ps = fhicl::ParameterSet::make(cliopts::fclname, *fm);
#else
  ps = fhicl::make_ParameterSet(cliopts::fclname);
#endif
  return ps;
}

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);
  if (!cliopts::fclname.size()) {
    std::cout << "[ERROR]: Expected to be passed a -c option." << std::endl;
    SayUsage(argv);
    return 1;
  }
  if (!cliopts::genie_input.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

#ifndef NO_ART
  fhicl::ParameterSet ps = ReadParameterSet(argv);
  std::vector<std::unique_ptr<SystProv>> syst_providers;
  syst_providers =
      systtools::ConfigureISystProvidersFromParameterHeaders<SystProv>(
          ps.get<fhicl::ParameterSet>(cliopts::fhicl_key));

  if (!syst_providers.size()) {
    throw response_helper_found_no_parameters()
        << "[ERROR]: Expected to load some systematic providers from input: "
        << std::quoted(cliopts::fclname);
  }

  systtools::param_header_map_t configuredParameterHeaders =
      systtools::BuildParameterHeaders(syst_providers);

  if (!configuredParameterHeaders.size()) {
    throw response_helper_found_no_parameters()
        << "[ERROR]: Expected systematric providers loaded from input: "
        << std::quoted(cliopts::fclname)
        << " to provide some parameter headers.";
  }

  ParamHeaderHelper phh(configuredParameterHeaders);
#else
  response_helper phh(cliopts::fclname);
#endif

  TChain *gevs = new TChain("gtree");
  if (!gevs->Add(cliopts::genie_input.c_str())) {
    std::cout << "[ERROR]: Failed to find any TTrees named "
              << std::quoted("gtree") << ", from TChain::Add descriptor: "
              << std::quoted(cliopts::genie_input) << "." << std::endl;
    return 3;
  }

  size_t NEvs = gevs->GetEntries();

  if (!NEvs) {
    std::cout << "[ERROR]: Input TChain contained no entries." << std::endl;
    return 4;
  }

  genie::NtpMCEventRecord *GenieNtpl = nullptr;

  if (gevs->SetBranchAddress("gmcrec", &GenieNtpl) != TTree::kMatch) {
    std::cout << "[ERROR]: Failed to set branch address on ghep tree."
              << std::endl;
    return 5;
  }

  TweakSummaryTree tst(cliopts::outputfile.c_str());
  tst.AddBranches(phh);

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  size_t NToRead = std::min(NEvs, cliopts::NMax);
  size_t NToShout = NToRead / 20;
  NToShout = NToShout ? NToShout : 1;
  for (size_t ev_it = 0; ev_it < NToRead; ++ev_it) {
    gevs->GetEntry(ev_it);
    genie::EventRecord const &GenieGHep = *GenieNtpl->event;

    tst.Mode = genie::utils::ghep::NeutReactionCode(&GenieGHep);
    tst.Emiss = GetEmiss(GenieGHep, false);
    tst.Emiss_preFSI = GetEmiss(GenieGHep, true);
    tst.pmiss = GetPmiss(GenieGHep, false);
    tst.pmiss_preFSI = GetPmiss(GenieGHep, true);
  
	// TH: add in select number of variable calculations here
    //genie::Target const &tgt = GenieGHep.Summary()->InitState().Tgt();
    //genie::GHepParticle *FSLep = GenieGHep.FinalStatePrimaryLepton();
    //genie::GHepParticle *ISLep = GenieGHep.Probe();
    //TLorentzVector FSLepP4 = *FSLep->P4();
    //TLorentzVector ISLepP4 = *ISLep->P4();
    //TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
    //tst.Mode = GenieGHep.EventGenerationMode();

    if (!(ev_it % NToShout)) {
      std::cout << (ev_it ? "\r" : "") << "Event #" << ev_it << "/" << NToRead
                << ", Interaction: " << GenieGHep.Summary()->AsString()
                << std::flush;
    }

    tst.Clear();


#ifndef NO_ART
    event_unit_response_w_cv_t resp;
    for (auto &sp : syst_providers) {
      systtools::ExtendEventUnitResponse(
          resp, sp->GetEventVariationAndCVResponse(GenieGHep));
    }
#else
    event_unit_response_w_cv_t resp =
        phh.GetEventVariationAndCVResponse(GenieGHep);
#endif

    tst.Add(resp);
    tst.Fill();
  }
  std::cout << std::endl;
}
