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

#include "nusystematics/utility/response_helper.hh"

#include "fhiclcpp/ParameterSet.h"

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
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
using namespace systtools;
using namespace nusyst;
using namespace genie;
using namespace genie::rew;

NEW_SYSTTOOLS_EXCEPT(unexpected_number_of_responses);

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
  std::vector<int> fsi_pdgs;
  std::vector<int> fsi_codes;
  // TKI
  std::vector<double> fsprotons_KE;
  double deltaPT, deltaalphaT;
  // Experimenting signal selection here..
  bool IsSignal_ICARUS_1muNp0pi;

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
    t->Branch("fsi_pdgs", "vector<int>", &fsi_pdgs);
    t->Branch("fsi_codes", "vector<int>", &fsi_codes);
    t->Branch("fsprotons_KE", "vector<double>", &fsprotons_KE);
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
        // For a correction dial, hdr.paramVariations is empty, so manually fill the vector
        if (hdr.isCorrection) {
          meta_tweak_values[0] = hdr.centralParamValue;
        } else {
          std::copy_n(hdr.paramVariations.begin(), meta_n,
                      meta_tweak_values.begin());
        }

        m->Fill();
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

namespace cliopts {
std::string fclname = "";
std::string genie_input = "";
std::string genie_branch_name = "gmcrec";
std::string outputfile = "";
std::string envvar = "FHICL_FILE_PATH";
std::string fhicl_key = "generated_systematic_provider_configuration";
size_t NMax = std::numeric_limits<size_t>::max();
size_t NSkip = 0;
#ifndef NO_ART
int lookup_policy = 1;
#endif
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-c <config.fcl>  : fhicl file to read.\n"
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

int main(int argc, char const *argv[]) {
  HandleOpts(argc, argv);

  bool RunNuSyst = true;
  if (!cliopts::fclname.size()) {
    RunNuSyst = false;
    std::cout << "-c is not given, running without evaluating reweights" << std::endl;
  }
  if (!cliopts::genie_input.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

  response_helper phh;
  if(RunNuSyst){
    phh = response_helper(cliopts::fclname);
  }

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

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile(
      "Messenger_whisper.xml");

  size_t NToRead = std::min(NEvs, cliopts::NMax);
  size_t NToShout = NToRead / 20;
  NToShout = NToShout ? NToShout : 1;
  for (size_t ev_it = cliopts::NSkip; ev_it < NToRead; ++ev_it) {
    // Start-of-event: reset tweak containers so later Add() fills cleanly.
    // (This does NOT touch physics scalars like w, q0, Q2.)
    tst.Clear();
    // Optional: set physics sentinels to catch logic mistakes in quick scans.
    tst.w   = -999.0;

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

// ------------------------------------------------------------
    // Compute W from (q0, Q2) and a nucleon mass.
    // W^2 = M_N^2 + 2 M_N q0 - Q^2, with q0 = energy transfer (nu)
    // Use GENIE's nucleon mass if available; otherwise default to 0.939 GeV.
    // ------------------------------------------------------------
    {
      double MN = 0.939;               // default GeV
  if (nucleon) MN = nucleon->Mass();  // use GENIE's struck-nucleon mass if available

  const double W2 = MN*MN + 2.0*MN*tst.q0 - tst.Q2;
  tst.w = (W2 > 0.0) ? std::sqrt(W2) : -1.0; // guard small negatives
    }

    
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

    // loop over particles
    int ip=-1;
    GHepParticle * p = 0;
    TIter event_iter(&GenieGHep);

    std::vector<int> fsi_pdgs;
    std::vector<int> fsi_codes;

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
      }

      // Stable final state particle
      if(ist==genie::kIStStableFinalState){
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
    for(const auto& proton: protons){
      const double this_KE = proton->KinE();
      tst.fsprotons_KE.push_back(this_KE);
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

    if (!(ev_it % NToShout)) {
      std::cout << (ev_it ? "\r" : "") << "Event #" << ev_it << "/" << NToRead
                << ", Interaction: " << GenieGHep.Summary()->AsString()
                << std::flush;
    }

    

    // Calcuate weights
    if(RunNuSyst){
      event_unit_response_w_cv_t resp = phh.GetEventVariationAndCVResponse(GenieGHep);
      tst.Add(resp);
    }
    // --- Compute W from (q0, Q2) and nucleon mass *right before Fill* ---
    //    W^2 = M_N^2 + 2 M_N q0 - Q^2  (q0 = energy transfer)
    {
      // Force single-nucleon mass for consistent comparison (no cluster mass)
      const double MN = 0.939; // GeV
      const double W2 = MN * MN + 2.0 * MN * tst.q0 - tst.Q2;
      tst.w = (W2 > 0.0) ? std::sqrt(W2) : -1.0;
    }

    // Optional quick debug every ~50k events to prove values are changing:
    // if ((ev_it % 50000) == 0) {
    //   std::cout << "\n[ev " << ev_it << "] q0=" << tst.q0
    //             << " Q2=" << tst.Q2 << " W=" << tst.w << std::endl;
    // }
    tst.Fill();
    
    // TH: Very important to clear this object to avoid memory issues!
    GenieNtpl->Clear();

  }
  std::cout << std::endl;
}
