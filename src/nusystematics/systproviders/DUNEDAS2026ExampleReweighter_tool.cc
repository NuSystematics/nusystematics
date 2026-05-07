#include "nusystematics/systproviders/DUNEDAS2026ExampleReweighter_tool.hh"

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/responsecalculators/DUNEDAS2026ExampleReweighter_calculator.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

DUNEDAS2026ExampleReweighter::DUNEDAS2026ExampleReweighter(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_DialA(systtools::kParamUnhandled<size_t>),
      pidx_DialB(systtools::kParamUnhandled<size_t>){

}

SystMetaData DUNEDAS2026ExampleReweighter::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[DUNEDAS2026ExampleReweighter::BuildSystMetaData] Called" << std::endl;

  SystMetaData smd;

  // Name of the dials that are supported by this module
  std::vector<std::string> AvailPNames = {"DialA", "DialB"};

  // Loop over available names and check if they are specified in ToolConfig
  for(std::string const &pname: AvailPNames){
    systtools::SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(cfg, pname, phdr, firstId)) {
      printf("[DUNEDAS2026ExampleReweighter::BuildSystMetaData] %s is found from ToolConfig\n", pname.c_str());
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }
  if(smd.size()==0){
    std::cout << "[DUNEDAS2026ExampleReweighter::BuildSystMetaData] No dial is set" << std::endl;
  }

  // You can extra parameters for the module;
  // Use get<T> function to retrieve the value from ToolConfig,
  // and then run "put" on "tool_options" (defined in the header of this class as a  member variable)
  // T can be string, bool, int, unsigned, float, double, std::string or even a new fhicl::ParameterSet

  std::string OPT_STRING = cfg.get<std::string>("OPT_STRING", ""); // second argument is the default when OPT_STRING does not exist
  tool_options.put("OPT_STRING", OPT_STRING);

  bool OPT_BOOL = cfg.get<bool>("OPT_BOOL", false);
  tool_options.put("OPT_BOOL", OPT_BOOL);

  fhicl::ParameterSet OPT_PSET = cfg.get<fhicl::ParameterSet>("OPT_PSET");
  tool_options.put("OPT_PSET", OPT_PSET);

  return smd;
}

bool DUNEDAS2026ExampleReweighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[DUNEDAS2026ExampleReweighter::SetupResponseCalculator] Called" << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();

  if(HasParam(md, "DialA")){
    pidx_DialA = GetParamIndex(md, "DialA");
  }
  if(HasParam(md, "DialB")){
    pidx_DialB = GetParamIndex(md, "DialB");
  }

  // Parameters in tool_options
  std::string OPT_STRING = tool_options.get<std::string>("OPT_STRING", "");
  bool OPT_BOOL = tool_options.get<bool>("OPT_BOOL", false);
  fhicl::ParameterSet OPT_PSET = tool_options.get<fhicl::ParameterSet>("OPT_PSET");

  return true;
}

event_unit_response_t
DUNEDAS2026ExampleReweighter::GetEventResponse(genie::EventRecord const &ev) {

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  //------ DDAS Exercise 1-2 START
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector FSLepP4 = *FSLep->P4();
  double Q2 = 0;
  //------ DDAS Exercise 1-2 END

  // now make the output
  // 1) Make an empty object
  systtools::event_unit_response_t resp;
  systtools::SystMetaData const &md = GetSystMetaData();

  // If pidx_DialA is found and set from SetupResponseCalculator,
  // it must be different from systtools::kParamUnhandled<size_t>.
  // Then we evaluate the reweight for DialA
  if (pidx_DialA != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_DialA].systParamId, {}} );
    for (double var : md[pidx_DialA].paramVariations) {
      // var is pariations (e.g., -1, 0, 1...)
      resp.back().responses.push_back( GetReweight_DialA(Q2, var) );
    } 
  }
  // Same for DialB
  if (pidx_DialB != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_DialB].systParamId, {}} );
    for (double var : md[pidx_DialB].paramVariations) {
      resp.back().responses.push_back( GetReweight_DialB(Q2, var) );
    }
  }

  return resp;
}

DUNEDAS2026ExampleReweighter::~DUNEDAS2026ExampleReweighter() {
} 
