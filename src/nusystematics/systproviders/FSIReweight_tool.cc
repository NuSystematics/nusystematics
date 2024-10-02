#include "nusystematics/systproviders/FSIReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

FSIReweight::FSIReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      fsiReweightCalculator(nullptr),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData FSIReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[FSIReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "FSIReweight",
                                                 phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet templateManifest =
      cfg.get<fhicl::ParameterSet>("FSIReweight_input_manifest");
  tool_options.put("FSIReweight_input_manifest", templateManifest);

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  save_map = cfg.get<bool>("save_map", false);
  tool_options.put("save_map", save_map);

  return smd;
}

bool FSIReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[FSIReweight::SetupResponseCalculator] called" << std::endl;

  fhicl::ParameterSet templateManifest =
      tool_options.get<fhicl::ParameterSet>("FSIReweight_input_manifest");

  fsiReweightCalculator = std::make_unique<FSIReweightCalculator>( templateManifest );

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "FSIReweight");

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  save_map = tool_options.get<bool>("save_map", false);
  if (save_map) {
    std::cout << "[[FSIReweight::SetupResponseCalculator] save_map is True" << std::endl;
    outfile_map = new TFile("FSI_2Dmap.root", "recreate");
    h_KEini_Ebias = new TH2D("h_KEini_Ebias", "Ebias vs KEini; KEini [GeV]; Ebias [GeV]", 40, 0, 2, 40, -0.9, 1.1);
  }


  return true;
}

event_unit_response_t
FSIReweight::GetEventResponse(genie::EventRecord const &ev) {

  // when the event is not applicable for this type of reweighting,
  // use GetDefaultEventResponse() to return an auto-1.-filled vector

  //if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
  //    !ev.Summary()->ProcInfo().IsWeakCC()) {
  //  return this->GetDefaultEventResponse();
  //}

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  genie::GHepParticle *ISNuc = ev.TargetNucleus();
  genie::GHepParticle *ISHitNucleon = ev.HitNucleon();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  //std::cout<<ev.Summary()->AsString()<<std::flush;
  //ev.Print();
  //std::cout<<std::endl;

  vector<genie::GHepParticle*> preFSIhadron_list;
  /*for (int ip=ISHitNucleon->FirstDaughter(); ip<=ISHitNucleon->LastDaughter(); ip++) {
    genie::GHepParticle *preFSIhadron = ev.Particle(ip);
    cout<<preFSIhadron->Pdg()<<"\t"<<preFSIhadron->Status()<<endl;
    if (preFSIhadron->Status() == 14)
      preFSIhadron_list.push_back(preFSIhadron);
  }*/
  TObjArrayIter piter(&ev);
  genie::GHepParticle * par = 0;
  while( (par = (genie::GHepParticle *) piter.Next()) ) {
    if (par->Status() == 14) {
      if (ev.Particle(par->FirstMother())->Status() == 14) continue;
      preFSIhadron_list.push_back(par);
    }
  } // get all pre-FSI hadrons
  
  TLorentzVector FSLepP4 = *FSLep->P4(); // l
  TLorentzVector ISLepP4 = *ISLep->P4(); // nu
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double AngleLeps = FSLepP4.Vect().Angle( ISLepP4.Vect() );
  double CAngleLeps = TMath::Cos(AngleLeps);
  double SAngleLeps = TMath::Sin(AngleLeps);
  if(AngleLeps>=M_PI/2.) SAngleLeps *= -1.;

  std::array<double, 2> bin_kin;
  bin_kin = {emTransfer.Vect().Mag(), emTransfer.E()};

  // now make the output
  systtools::event_unit_response_t resp;

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];
  resp.push_back( {hdr.systParamId, {}} );
  vector<int> npar_sel = {1, 0, 0, 0, 0}; // edit the number of {p, n, pip, pi0, pim} (temporarily added inline for test, should be added as parameters in fcl file)
  bool allevents = true; // if allevents==true, npar_sel will be disabled
  int np = 0;
  for (int num : npar_sel) { np += num; }
  map<int,int> npar_idx;
  npar_idx[2212] = 0; npar_idx[2112] = 1; npar_idx[211] = 2; npar_idx[111] = 3; npar_idx[-211] = 4;
  vector<int> npar = {0, 0, 0, 0, 0};
  for (double var : hdr.paramVariations) {
    double weight = 1;
    if ((!allevents)&&!(preFSIhadron_list.size() == np)) { // selection on the number of pre-FSI particles
      weight = -199;
      //weight = 1; // do not reweight events with no pre-FSI particles at all
    }
    else {
      double KEini_0, KEini_1;
      double Ehad_tot = 0;
      int idx = 0;
      bool rewei = true;
      for (auto IShad : preFSIhadron_list) {
        //if (!(IShad->Pdg()==2212||IShad->Pdg()==211)) { // selection on proton-only for test
        //  weight = -99;
        //  rewei = false;
        //  break;
        //}
        int had_pdg = npar_idx[ IShad->Pdg() ];
        ++npar[had_pdg];
        if ((!allevents)&&(npar[had_pdg] > npar_sel[had_pdg])) {
          weight = -99;
          rewei = false;
          break;
        }
        
        if (IShad->FirstDaughter()==-1) {
          //cout<<"No daughter particle for this pre-FSI hadron..."<<endl; // due to no daughter info in INCL and G4BC
          continue;
        }
        vector<genie::GHepParticle*> FSdaughter_list;
        get_FS_daughters(IShad, FSdaughter_list, ev);
        double Eini = IShad->E();
        double KEini = IShad->KinE();
        double Ehad = 0;
        for (auto FSpar : FSdaughter_list) {
          int pdg = FSpar->Pdg();
          double totE = FSpar->E();
          double kinE = FSpar->KinE();
          if ( pdg==211 || pdg==-211 || pdg==111) // pion
            Ehad += totE;
          else if ( pdg==2212 ) // proton
            Ehad += kinE;
          else if ( pdg==2112 ) // neutron
            ;
          else if ( pdg==22 ) // gamma
            Ehad += totE;
          //else cout<<"$#@ "<<FSpar->Name()<<endl;
        }
        double Ebias = 0;
        if (idx==0) {
          KEini_0 = KEini;
          idx++;
        }
        else {
          KEini_1 = KEini;
        }
        
        if (IShad->RescatterCode()==1) { // non-FSI
          //weight = 1.6;
          continue;
        }
        Ehad_tot += Ehad;
        switch (IShad->Pdg()) {
          case 2212:
            Ebias = (KEini - Ehad) / KEini;
            break;
          case 2112:
            Ebias = (KEini - Ehad) / KEini;
            break;
          case 211:
            Ebias = (Eini - Ehad) / Eini;
            break;
          case 111:
            Ebias = (Eini - Ehad) / Eini;
            break;
          case -211:
            Ebias = (Eini - Ehad) / Eini;
            break;
          default:
            //weight = -99;
            //break;
            weight = 1;
            continue;
        }
        if (save_map) h_KEini_Ebias->Fill( KEini, Ebias );
        //cout<<IShad->Name()<<": KEini "<<KEini<<"; Ehad "<<Ehad<<"; Ebias "<<Ebias<<endl;
        double this_reweight = fsiReweightCalculator->GetFSIReweight(KEini, Ebias, var, IShad->Pdg());
        weight *= this_reweight;
        //cout<<IShad->Name()<<": KEini "<<KEini<<"; Ebias "<<Ebias<<"; weight "<<this_reweight<<endl;
      }
      //if (rewei) {
      //  double tmp = KEini_0;
      //  if (tmp<KEini_1) {
      //    KEini_0 = KEini_1;
      //    KEini_1 = tmp;
      //  }
      //  double Ebias = (KEini_0+KEini_1-Ehad_tot)/(KEini_0+KEini_1);
      //  weight = fsiReweightCalculator->GetFSIReweight_2par(KEini_0, KEini_1, Ebias, var, 2);
      //}
    }
    resp.back().responses.push_back( weight );
  }

  if (fill_valid_tree) {

    pdgfslep = ev.FinalStatePrimaryLepton()->Pdg();
    momfslep = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();

    Pdgnu = ISLep->Pdg();
    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

    QELikeTarget_t qel_targ = GetQELikeTarget(ev);
    QELTarget = e2i(qel_targ);

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    valid_tree->Fill();
  }

  return resp;

}

std::string FSIReweight::AsString() { return ""; }

void FSIReweight::InitValidTree() {

  valid_file = new TFile("MINERvAq3q0WeightValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("QELTarget", &QELTarget);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
}

FSIReweight::~FSIReweight() {
  if (save_map) {
    h_KEini_Ebias->Write("hA2018_proton_KEini_vs_Ebias");
    outfile_map->Write();
  }
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
