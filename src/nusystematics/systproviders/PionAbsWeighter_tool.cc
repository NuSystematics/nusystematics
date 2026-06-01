#include "nusystematics/systproviders/PionAbsWeighter_tool.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "nusystematics/utility/exceptions.hh"
///include "RwFramework/GSyst.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/HadronTransport/INukeHadroFates2018.h"

using namespace nusyst;
using namespace systtools;
using namespace fhicl;

//std::array<std::string, 4> ParamPrettyNames = {"p1", "p2", "p3", "p4"};

// constructor passes up configuration object to base class for generic tool
// initialization and initialises our local copies of paramIds to unconfigured
//  flag values
PionAbsWeighter::PionAbsWeighter(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params),
    PionAbsCalculator(nullptr)
{
    for( std::vector<std::string>::iterator it_desc = descriptors.begin();
       it_desc != descriptors.end(); ++it_desc ) {
    ResponseParameterIndices.emplace_back( kParamUnhandled<size_t> );
  }

}     
// pidx_Params{kParamUnhandled<size_t>, kParamUnhandled<size_t>,
     //             kParamUnhandled<size_t>, kParamUnhandled<size_t>} {}

PionAbsWeighter::~PionAbsWeighter() {}

SystMetaData PionAbsWeighter::BuildSystMetaData(fhicl::ParameterSet const &cfg,
                                              paramId_t id) {
  SystMetaData smd;
  //SystParamHeader responseParam;

  /*
  // loop through the four named parameters that this tool provides
  for (std::string const &pname : ParamPrettyNames) {
    SystParamHeader phdr;

    // Set up parameter definition with a standard tool configuration form
    // using helper function.
    if (ParseFhiclToolConfigurationParameter(ps, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;

      // set any parameter-specific ParamHeader metadata here
      phdr.isSplineable = true;

      // add it to the metadata list to pass back.
      smd.push_back(phdr);
    }
  }//loop over parameters
*/
 //Use descriptors instead of PrettyNames
  // Obtain each named parameter from the descriptors in the fcl
  for( std::vector<std::string>::iterator it_desc = descriptors.begin();
       it_desc != descriptors.end(); ++it_desc ) {
    SystParamHeader phdr;
    std::string pname = *it_desc;
    if( ParseFhiclToolConfigurationParameter(cfg, pname, phdr, id) ) {
      if( verbosity_level > 4 ) {
	std::cout << "[DEBUG]: Found parameter " << pname << " with id = " << id << std::endl;
      } // verbose output
      phdr.systParamId = id++; // increment id here
      smd.push_back(phdr);
    }
  } // loop over named parameters	

    ParameterSet templateManifest =
      cfg.get<ParameterSet>("PionAbs_input_manifest");

  if (!cfg.has_key("PionAbs_input_manifest") ||
      !cfg.is_key_to_table("PionAbs_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for "
           "PionAbs, expected to find a FHiCL table keyed by "
           "PionAbs_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }

  tool_options.put("PionAbs_input_manifest", templateManifest);

  // Put any options that you want to propagate to the ParamHeaders options
  tool_options.put("verbosity_level", cfg.get<int>("verbosity_level", 0));

  return smd;
}

bool PionAbsWeighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  
  // Silence GENIE
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  // loop through the named parameters that this tool provides, check that they
  // are configured, and grab their id in the current systmetadata and set up
  // and pre-calculations/configurations required.
  //for (size_t i = 0; i < ParamPrettyNames.size(); ++i) {

    /*if (!HasParam(md, ParamPrettyNames[i])) {
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter " << ParamPrettyNames[i]
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }
    */

    //Use Descriptors instead of PrettyNames:
    for( std::vector<std::string>::iterator it_desc = descriptors.begin();
       it_desc != descriptors.end(); ++it_desc ) {
    if( !HasParam( md, *it_desc ) ) {
      continue;
    }

    //pidx_Params[i] = GetParamIndex(md, ParamPrettyNames[i]);
    
       // Get the parameter index
    ResponseParameterIndices[it_desc - descriptors.begin()] = GetParamIndex(md, *it_desc);
    if( verbosity_level > 1 ) {
      std::ostringstream asts;
      std::vector<double> param_variations = md[GetParamIndex(md, *it_desc)].paramVariations;
      asts << "[INFO]: Configured parameter " << *it_desc << " with variations:\n\t[";
      for( double & var: param_variations ) asts << " " << var << ",";
      asts << " ]";
      std::cout << asts.str() << std::endl;
    } // verbose output about configuration 

    //auto param_md = md[pidx_Params[i]];

    //CVs[i] = param_md.centralParamValue;
    //Variations[i] = param_md.paramVariations;

      // Get manifests for options
  if (!tool_options.has_key("PionAbs_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: PionAbs parameter exists in the "
           "SystMetaData, but "
           "no PionAbs_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "PionAbsWeighter::BuildSystMetaData, but wasn't, this is a "
           "bug, "
           "please report to the maintainer.";
  }

  fhicl::ParameterSet const &templateManifest =
      tool_options.get<fhicl::ParameterSet>(
          "PionAbs_input_manifest");

  // Initialise the calculator
  PionAbsCalculator = std::make_unique<PionAbsResponseCalculator>(templateManifest);

      } //loop over descriptors
  // returning cleanly
  return true;
}

event_unit_response_t
PionAbsWeighter::GetEventResponse(genie::EventRecord const &ev) {
 
  //if (verbosity_level > 4) std::cout<<"In GetEventResponse"<<std::endl;

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  //want to explicitly fill the weight vector with 1 if our event isn't resonant and absorption
  bool resFlag = false;
  bool absFlag = false;
  
  if (ev.Summary()->ProcInfo().IsResonant() ) {
    resFlag = true;
  }
  
  if (verbosity_level > 5 && resFlag) std::cout<<"Event is RES"<<std::endl;

  //demand that the event contain no pion in the final state:
  int NfPip      = 0; // number of \pi^+'s         in final state
  int NfPim      = 0; // number of \pi^-'s         in final state
  int NfPi0      = 0; // number of \pi^0's         in final state
  
  TObjArrayIter piter(&ev);
  genie::GHepParticle * p = 0;
  int ip=-1;

  while( (p = (genie::GHepParticle *) piter.Next())) {
    ip++;
    int pdgc = p->Pdg();
    int ist  = p->Status();
    // only final state particles
    if(ist != genie::kIStStableFinalState) continue;
    // don't count final state lepton as part of the hadronic system
    if(ev.Particle(ip)->FirstMother()==0) continue;

    if (pdgc == genie::kPdgPiP        ) NfPip++;
    if (pdgc == genie::kPdgPiM        ) NfPim++;
    if (pdgc == genie::kPdgPi0        ) NfPi0++;

  }

  //Now do a loop to check for the existence of intermediate state pions
  int niPip = 0, niPim = 0, niPi0 = 0;

  int pionAbsIndex = -1;
  //does a pion undergo an absorption interaction?
  TObjArrayIter mpiter(&ev);
      genie::GHepParticle * part = 0;
      int imp=-1;
      while( (part = (genie::GHepParticle *) mpiter.Next())) {
        imp++;
        int pdgc = part->Pdg();
        int ist  = part->Status();
        // only pions
	if(std::abs(pdgc) != genie::kPdgPiP && pdgc != genie::kPdgPi0) continue; //require an intermediate pion
        if(ist != genie::kIStHadronInTheNucleus /* && ist != genie::kIstDISPreFragmHadronicState */ ) continue;
                switch( pdgc ) {
                        case genie::kPdgPiP : niPip++; break;
                        case genie::kPdgPiM : niPim++; break;
                        case genie::kPdgPi0 : niPi0++; break;
                }
        // look at the rescattering code
        // It's super annoying this returns an int...
        // Define pion abs by the presence of an intermediate pion which absorbs
        genie::INukeFateHA_t rescat = static_cast<genie::INukeFateHA_t>( part->RescatterCode() );
        if( rescat == genie::kIHAFtAbs ) {absFlag = true; pionAbsIndex = imp;} //technically takes the last pion abs in an event. Almost never 2 in 1 event.
	}
   //Define pion absorption based on the presence of intermediate state pion and absence of final state pion
   //if ( (niPip + niPim + niPi0 > 0) && (NfPip + NfPim + NfPi0 == 0) ) isPionAbs = true;

  if (verbosity_level > 5 && absFlag) std::cout<<"We have a pion abs event."<<std::endl;
  //Now that we know we have a pion abs, assign the reweights:
  // loop through and calculate weights
   for( std::vector<std::string>::iterator it_desc = descriptors.begin();
       it_desc != descriptors.end(); ++it_desc ) {
 
   // this parameter wasn't configured, nothing to do
    if (ResponseParameterIndices[it_desc - descriptors.begin()] == kParamUnhandled<size_t>) {
      continue;
    }    
    // initialize the response array with this paramId
    resp.push_back({md[ResponseParameterIndices[it_desc - descriptors.begin()]].systParamId, {}});

    if (verbosity_level>5) std::cout<<"Initialized response array"<<std::endl;
    // loop through variations for this parameter
    std::vector<double> param_variations = md[ResponseParameterIndices[it_desc - descriptors.begin()]].paramVariations;
    for (size_t v_it = 0; v_it < param_variations.size(); ++v_it) {
	if(verbosity_level>5) std::cout<<"In variation "<<v_it<<std::endl;

      // put the response weight for this variation of this parameter into the
      // response object
	double this_weight;
	double this_variation = param_variations[v_it];

    	//calculate the weight here
    	if (resFlag && absFlag) {
    		if ( it_desc - descriptors.begin() == 0) this_weight = PionAbsWeighter::ReweighthADiff(ev,pionAbsIndex,this_variation,false);
    		if ( it_desc - descriptors.begin() == 1) this_weight = PionAbsWeighter::ReweighthADiff(ev,pionAbsIndex,this_variation,true);
    		if ( it_desc - descriptors.begin() == 2) this_weight = PionAbsWeighter::ReweightQDFraction(ev,pionAbsIndex,this_variation);	
	}
	else this_weight = 1;
      resp.back().responses.push_back( this_weight );


      if (verbosity_level > 3) {
        std::cout << "[DEBG]: For parameter " << descriptors.at(it_desc - descriptors.begin())
                  << " at variation[" << v_it << "] = " << param_variations[v_it]
                  << " calculated weight: " << resp.back().responses.back()
                  << std::endl;
      }//if verbose
    }//loop over variations
  }//loop over parameters

  return resp;
}//GetEventResponse


//helper functions:

//vector of stable final state nucleon indices descended from the absorbed pion
std::vector<int> PionAbsWeighter::IndicesOfAbsNucleons(genie::EventRecord const &ev, int pionAbsIndex) {
   
    std::vector<int>this_IndicesOfAbsNucleons;
    for (int i=0; i<ev.GetEntries(); ++i) {
        const auto* p = ev.Particle(i);

	//particle must be a nucleon
        if (!(p->Pdg() == 2212 || p->Pdg() == 2112)) continue;
	if (!(p->Status() == genie::kIStStableFinalState)) continue;

        bool ancestorIsAbsorbedPion = false;
        if (!ancestorIsAbsorbedPion) {
            int mom = p->FirstMother();
            while (mom >= 0) {
                const auto* m = ev.Particle(mom);
                if (!m) break;
                if (*m == *(ev.Particle(pionAbsIndex))) {ancestorIsAbsorbedPion = true; break;}
                mom = m->FirstMother();
            }
        }
        if (ancestorIsAbsorbedPion) this_IndicesOfAbsNucleons.push_back(i);
    }
return this_IndicesOfAbsNucleons;
}

//difference of multiplicity of absorption nucleons above threshold
int PionAbsWeighter::DiffAbsNucAboveThresh(genie::EventRecord const &ev, int pionAbsIndex, double threshold) {
    std::vector<int> this_IndicesOfAbsNucleons = PionAbsWeighter::IndicesOfAbsNucleons(ev, pionAbsIndex);

    int num_p = 0;
    int num_n = 0;

    //loop over absorption nucleons to see if there is a nuclear cluster parent
    for (int nuc_index : this_IndicesOfAbsNucleons) {
        const auto* p = ev.Particle(nuc_index);
	double this_ke = p->Energy() - p->Mass();
	bool above_threshold = (this_ke >= threshold);
	if (p->Pdg() == 2212 && above_threshold) num_p++;
	if (p->Pdg() == 2112 && above_threshold) num_n++;
    }

    return num_p - num_n; 
}

std::string PionAbsWeighter::LabelBase(genie::EventRecord const &ev,int pionAbsIndex) {
    int probe = ev.Particle(pionAbsIndex)->Pdg();
    std::string probeStr = (probe==211  ? "piPlus"                                                   
 		          : probe==-211 ? "piMinus"
                          : "pi0");
    std::string keStr;
    float ke = ev.Particle(pionAbsIndex)->Energy() - ev.Particle(pionAbsIndex)->Mass();

    //round to the nearest 0.1 GeV

    float roundedKE = std::round(ke * 10.0f) / 10.0f;

    // Clamp to [0.1, 1.0]
    if (roundedKE < 0.1f) roundedKE = 0.1f;
    if (roundedKE > 1.0f) roundedKE = 1.0f;

    if (std::fabs(roundedKE - 1.0) < 1e-6) {
       keStr = "1.0";
       } else {
       keStr = Form("%g", roundedKE);
       }
    
    int target = 1000180400;
    std::string labelBase = Form("%s_%sGeV_%i", probeStr.c_str(), keStr.c_str(), target);
    return labelBase;
}

double PionAbsWeighter::ReweighthADiff(genie::EventRecord const &ev,int pionAbsIndex,double this_variation,bool is_INCL) {
    double weight = 1.0;
    
    const double threshold = 2e-3; //20 MeV, assumes energy in GeV
    //calculate the number of nucleons with KE above threshold
    int diff = PionAbsWeighter::DiffAbsNucAboveThresh(ev, pionAbsIndex, threshold);

    int pionPDG = 0;
    float KEpi_GeV = 0;

    TObjArrayIter mpiter(&ev);
      genie::GHepParticle * part = 0;
      int imp=-1;
      while( (part = (genie::GHepParticle *) mpiter.Next())) {
        imp++;
        if (imp == pionAbsIndex){
        pionPDG = part->Pdg();
        KEpi_GeV = part->Energy() - part->Mass();
	break;
        }
	else continue;
      }
    
    if (verbosity_level > 5) std::cout<<"KEpi_GeV = "<<KEpi_GeV<<"; diff = "<<diff<<std::endl;
    //load the ROOT histogram from the file, read the bin content corresponding to diff
    double weight_central_value = PionAbsCalculator->GetWeight(is_INCL, false, pionPDG, KEpi_GeV, diff); //doSum = false

    if (verbosity_level > 5) std::cout<<"Back in ReweighthADiff: weight_central_value = "<<weight_central_value<<std::endl;
    if (verbosity_level > 5) std::cout<<"This variation = "<<this_variation<<std::endl;

    weight = 1 + (weight_central_value - 1)*this_variation;

    //return the reweight for that value of nucleon difference above threshold
    return weight;
}

double PionAbsWeighter::ScalableErf(double this_variation, double asymptotic_unc=0.5, double scale = 1.0) {
    double ev = std::erf(this_variation * scale);
    if (verbosity_level > 5) std::cout<<"this_variation: "<<this_variation<<"; scale: "<<scale<<"; asymptotic_unc: "<<asymptotic_unc<<"; ev: "<<ev<<"; return value: "<<(1-asymptotic_unc) + asymptotic_unc*(1+ev)<<std::endl;
    return 1 + asymptotic_unc*ev;
}

double PionAbsWeighter::ReweightQDFraction(genie::EventRecord const &ev,int pionAbsIndex,double this_variation) {
    bool isMNAbs = false;

    //Determine the value of MN reweight given QD reweight
    float ke = ev.Particle(pionAbsIndex)->Energy() - ev.Particle(pionAbsIndex)->Mass(); //ke is the pion ke in GeV
    if (verbosity_level > 5) std::cout<<"In ReweightQDFraction: pion ke = "<<ke<<std::endl;
    float fRemnA = 40.0;
    float f_QD = 1.14*(.903-0.00189*fRemnA)*(1.35-0.00467*(ke*1000)); //empirical QD fraction from hAcode, ke in GeV -> *1000 in MeV
    if (f_QD < 0 || f_QD > 1) return 1.0;
    if (verbosity_level > 5) std::cout<<"QD fraction: "<<f_QD<<std::endl; 
    
    double asymptotic_unc = 1.0;
    //If QD_Reweight > 1 + (1-f_QD)/f_QD, then MN_Reweight is negative. That's bad.
    if (asymptotic_unc > (1-f_QD)/f_QD) asymptotic_unc = (1-f_QD)/f_QD;

    double QD_Reweight = PionAbsWeighter::ScalableErf(this_variation, asymptotic_unc, 0.55);
    //This value makes sure to preserve the total number of pion abs events
    double MN_Reweight = (1 + (1 - QD_Reweight)*f_QD/(1-f_QD));

    //now determine if the event is QD or MN by checking if an absorption nucleon has a cluster ancestor
    //first grab the nucleons associated with the pion absorption:
    std::vector<int> this_IndicesOfAbsNucleons = PionAbsWeighter::IndicesOfAbsNucleons(ev, pionAbsIndex);

    //loop over absorption nucleons to see if there is a nuclear cluster parent
    for (int nuc_index : this_IndicesOfAbsNucleons) {
	const auto* p = ev.Particle(nuc_index);
	// Look for a cluster ancestor (if true -- MN; if false -- QD)
        bool sawCluster = false;
        if (!sawCluster) {
            int mom = p->FirstMother();
            while (mom >= 0) {
                const auto* m = ev.Particle(mom);
                if (!m) break;
                if (m->Pdg() == genie::kPdgCompNuclCluster) {sawCluster = true; break;}
                mom = m->FirstMother();
            }
	}
    if (sawCluster) isMNAbs = true;
    }//loop over particles

    if (isMNAbs) return MN_Reweight;
    else return QD_Reweight;
}

