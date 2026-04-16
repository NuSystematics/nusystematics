#include "nusystematics/systproviders/QEInterference_tool.hh"
#include "nusystematics/utility/exceptions.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"
#include "Framework/GHEP/GHepParticle.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

QEInterference::QEInterference(ParameterSet const & ps) :
  IGENIESystProvider_tool(ps),
  QEIntfResponseCalculator(nullptr)/*,
				     ResponseParameterIdx(kParamUnhandled<size_t>) {}*/
{
}

QEInterference::~QEInterference() {}

bool QEInterference::SetupResponseCalculator(ParameterSet const & tool_options)
{
  // Silence GENIE
  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_whisper.xml");

  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  q0BinEdges = tool_options.get<std::vector<double>>("q0_bin_edges");
  if( verbosity_level > 3 ) {
    std::cout << "[INFO]: Using bin edges: ";
    for(const auto q0BinEdge: q0BinEdges) std::cout << q0BinEdge << ", ";
    std::cout << std::endl;
  }

  Nq0Bins = q0BinEdges.size()-1;
  ResponseParameterIndices.resize(Nq0Bins);

  // Check the metadata makes sense
  SystMetaData const & md = GetSystMetaData();
  for(unsigned int i_q0bin=0; i_q0bin<Nq0Bins; i_q0bin++){
  
    std::string desc = "QEIntf_dial_"+std::to_string(i_q0bin);

    if( !HasParam( md, desc ) ) {
      throw incorrectly_configured()
	<< "[ERROR]: Expected to find parameter named "
	<< std::quoted(desc);
    } // Problem: Parameter not found.
    
    // Get the parameter index
    ResponseParameterIndices[i_q0bin] = GetParamIndex(md, desc);
    if( verbosity_level > 1 ) {
      std::ostringstream asts;
      std::vector<double> param_variations = md[GetParamIndex(md, desc)].paramVariations;
      asts << "[INFO]: Configured parameter " << desc << " with variations:\n\t[";
      for( double & var: param_variations ) asts << " " << var << ",";
      asts << " ]";
      std::cout << asts.str() << std::endl;
    } // verbose output about configuration
    
  } // loop over descriptors

  // Get manifests for options
  if (!tool_options.has_key("QEInterference_input_manifest")) {
    throw systtools::invalid_ToolOptions()
        << "[ERROR]: QEInterference parameter exists in the "
           "SystMetaData, but "
           "no QEInterference_input_manifest key can be found on the "
           "tool_options table. This reweighting requires input histograms "
           "that must be specified. This should have been caught by  "
           "QEInterference::BuildSystMetaData, but wasn't, this is a "
           "bug, "
           "please report to the maintainer.";
  }

  fhicl::ParameterSet const &templateManifest =
      tool_options.get<fhicl::ParameterSet>(
          "QEInterference_input_manifest");

  // Initialise the calculator
  QEIntfResponseCalculator = std::make_unique<QEInterferenceResponseCalculator>(templateManifest);

  return true;
}

SystMetaData QEInterference::BuildSystMetaData(ParameterSet const & cfg,
						     paramId_t id) 
{

  SystMetaData smd;

  std::vector<double> this_q0BinEdges = cfg.get<std::vector<double>>("q0_bin_edges");
  unsigned int this_Nq0Bins = this_q0BinEdges.size()-1;

  int this_verbosity_level = cfg.get<int>("verbosity_level", 0);
  if( this_verbosity_level > 3 ) {
    std::cout << "[INFO]: Using bin edges: ";
    for(const auto q0BinEdge: this_q0BinEdges) std::cout << q0BinEdge << ", ";
    std::cout << std::endl;
  }

  for(unsigned int i_q0bin=0; i_q0bin<this_Nq0Bins; i_q0bin++){

    std::string desc = "QEIntf_dial_"+std::to_string(i_q0bin);

    SystParamHeader phdr;
    std::string pname = desc;
    if( ParseFhiclToolConfigurationParameter(cfg, pname, phdr, id) ) {
      if( this_verbosity_level > 4 ) {
        std::cout << "[DEBUG]: Found parameter " << pname << " with id = " << id << std::endl;
      } // verbose output
      phdr.systParamId = id++; // increment id here
      smd.push_back(phdr);
    }
  } // loop over named parameters

  tool_options.put("q0_bin_edges", this_q0BinEdges);

  if( !cfg.has_key("QEInterference_input_manifest") ) {
    throw invalid_ToolConfigurationFHiCL() << "No q0 bin edges!!!";
  }
  ParameterSet templateManifest =
      cfg.get<ParameterSet>("QEInterference_input_manifest");

  if (!cfg.has_key("QEInterference_input_manifest") ||
      !cfg.is_key_to_table("QEInterference_input_manifest")) {
    throw invalid_ToolConfigurationFHiCL()
        << "[ERROR]: When configuring calculated variations for "
           "QEInterference, expected to find a FHiCL table keyed by "
           "QEInterference_input_manifest describing the location of "
           "the histogram inputs. See "
           "nusystematics/responsecalculators/"
           "TemplateResponseCalculatorBase.hh for the layout.";
  }

  tool_options.put("QEInterference_input_manifest", templateManifest);

  tool_options.put("verbosity_level", this_verbosity_level);

  return smd;
}

std::string QEInterference::AsString() { return "QEInterference"; }

event_unit_response_t QEInterference::GetEventResponse(genie::EventRecord const & ev) 
{ 
  // Anything but QE is not handled and should be default
  // Use GetDefaultEventResponse() to return an auto-1.-filled vector

  if ( !ev.Summary()->ProcInfo().IsQuasiElastic() ||
       (!ev.Summary()->ProcInfo().IsWeakCC() &&
	!ev.Summary()->ProcInfo().IsWeakNC()) ) {
    return this->GetDefaultEventResponse();
  }

  // Load up the Enu, Q0, and Q3. Also the PDG code!
  genie::GHepParticle * neutrino = ev.Probe();
  genie::GHepParticle * FSLep    = ev.FinalStatePrimaryLepton();

  if( ! neutrino || ! FSLep ) {
    throw incorrectly_generated()
      << "[ERROR]: Failed to find neutrino and final-state primary lepton in event: "
      << ev.Summary()->AsString();
  }

  TLorentzVector p4_nu  = *( neutrino->GetP4() );
  double Enu = p4_nu.E();
  int    pdg = neutrino->Pdg();
  int current = static_cast<int>( ev.Summary()->ProcInfo().IsWeakNC() ); // 0 if CC, 1 if NC
  
  TLorentzVector p4_lep = *( FSLep->GetP4() );
  
  // Q0 == Enu - Elep
  // Q3 == ||\vec{p3nu} - \vec{p3lep}||
  TLorentzVector delta_p4 = p4_nu - p4_lep;
  double Q0 = delta_p4.E();
  double Q3 = delta_p4.P();

  // Make the output
  event_unit_response_t resp;

  SystMetaData const & md = GetSystMetaData();

  // From the tweak dial and the strength, construct the weight.
  // First, get the raw response from the histogram
  double raw_response = QEIntfResponseCalculator->GetWeight( pdg, current, Enu, Q0, Q3 );

  // Initialise response arrays
  // We'll be careful here. Use resize to construct the vectors

  for(unsigned int i_q0bin=0; i_q0bin<Nq0Bins; i_q0bin++){

    SystParamHeader const & sph = md[ResponseParameterIndices[i_q0bin]];

    resp.push_back({sph.systParamId, {}});

    if(q0BinEdges[i_q0bin+1] < Q0 || q0BinEdges[i_q0bin] >= Q0) { // Default response
      std::vector<double> vec(sph.paramVariations.size(), 1.0);
      resp.back().responses = vec;
    } else {
      std::vector<double> vec;
      for( const double & twk : sph.paramVariations ) {
        double this_reweight = std::max( 0.0, (1.0 - twk) + twk * raw_response );
        vec.push_back(this_reweight);
      }
      resp.back().responses = vec;
    }
  }

  return resp;
}
