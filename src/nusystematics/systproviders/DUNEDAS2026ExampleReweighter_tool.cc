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
      pidx_DDAS(systtools::kParamUnhandled<size_t>){

}

SystMetaData DUNEDAS2026ExampleReweighter::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[DUNEDAS2026ExampleReweighter::BuildSystMetaData] Called" << std::endl;

  SystMetaData smd;

  std::string pname = "DIALNAME";
  systtools::SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, pname, phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  std::string OPT_STRING = cfg.get<std::string>("OPT_STRING", ""); // second argument is the default when OPT_STRING does not exist
  tool_options.put("OPT_STRING", OPT_STRING);

  bool OPT_BOOL = cfg.get<bool>("OPT_BOOL", false);
  tool_options.put("OPT_BOOL", OPT_BOOL);

  return smd;
}

bool DUNEDAS2026ExampleReweighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[DUNEDAS2026ExampleReweighter::SetupResponseCalculator] Called" << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();

  //------ DDAS TASK START
  // Check if the Dial 
  if(HasParam(md, "DIALNAME")){
    pidx_DDAS = GetParamIndex(md, "DIALNAME");
  }
  //------ DDAS TASK END

  return true;
}

event_unit_response_t
DUNEDAS2026ExampleReweighter::GetEventResponse(genie::EventRecord const &ev) {

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();

  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  // now make the output
  // 1) Make an empty object
  systtools::event_unit_response_t resp;
  systtools::SystMetaData const &md = GetSystMetaData();

  if (pidx_DDAS != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_DDAS].systParamId, {}} );
    for (double var : md[pidx_DDAS].paramVariations) {
      resp.back().responses.push_back( GetMyWeight(FSLepP4.E()) );
    } 
  }

  return resp;
}

DUNEDAS2026ExampleReweighter::~DUNEDAS2026ExampleReweighter() {
} 
