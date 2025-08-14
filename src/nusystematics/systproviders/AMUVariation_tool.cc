#include "nusystematics/systproviders/AMUVariation_tool.hh"
#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/string_parsers.hh"
#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

AMUVariation::AMUVariation(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      AMUCalculator(nullptr),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>) {}

SystMetaData AMUVariation::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[AMUVariation::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "AMUVariation",
                                                 phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet sf_Manifest =
      cfg.get<fhicl::ParameterSet>("AMUVariation_input_manifest");
  tool_options.put("AMUVariation_input_manifest", sf_Manifest);

  return smd;

}

bool AMUVariation::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[AMUVariation::SetupResponseCalculator] called" << std::endl;

  fhicl::ParameterSet sf_Manifest =
      tool_options.get<fhicl::ParameterSet>("AMUVariation_input_manifest");

  std::string F1FilePath = sf_Manifest.get<std::string>("F1FilePath");
  std::string F2FilePath = sf_Manifest.get<std::string>("F2FilePath");
  std::string F3FilePath = sf_Manifest.get<std::string>("F3FilePath");

  // Convert all environmental varaible
  F1FilePath = expand_env_vars(F1FilePath);
  F2FilePath = expand_env_vars(F2FilePath);
  F3FilePath = expand_env_vars(F3FilePath);

  AMUCalculator = std::make_unique<weightAMUvariations>(F1FilePath.c_str(), F2FilePath.c_str(), F3FilePath.c_str());

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "AMUVariation");

  return true;

}

event_unit_response_t
AMUVariation::GetEventResponse(genie::EventRecord const &ev) {

  // Check if applicable
  if (!ev.Summary()->ProcInfo().IsDeepInelastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC()) {
    return this->GetDefaultEventResponse();
  }

  // now make the output
  systtools::event_unit_response_t resp;

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  // Get one-sigma reweight first
  const genie::Kinematics &kin = ev.Summary()->Kine();
  double Bjorken_x = kin.x(true);
  double Bjorken_y = kin.y(true);
  // nu
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector ISLepP4 = *ISLep->P4();
  double Enu = ISLepP4.E();
  int is_antinu = ISLep->Pdg() > 0 ? 0 : 1;

  // charged lepton
  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  double Mlep = FSLep->Mass();
/*
  std::cout << "[AMUVariation::GetEventResponse] Bjorken_x = " << Bjorken_x << std::endl;
  std::cout << "[AMUVariation::GetEventResponse] Bjorken_y = " << Bjorken_y << std::endl;
  std::cout << "[AMUVariation::GetEventResponse] Enu = " << Enu << std::endl;
  std::cout << "[AMUVariation::GetEventResponse] Mlep = " << Mlep << std::endl;
  std::cout << "[AMUVariation::GetEventResponse] is_antinu = " << is_antinu << std::endl;
*/
  double weights[10];
  AMUCalculator->getWeight(Bjorken_x, Bjorken_y, Enu, Mlep, is_antinu, weights);
  double rw_RatioToNLO = weights[1];
  double onesig_RatioToNLO = rw_RatioToNLO - 1.;
/*
  std::cout << "[AMUVariation::GetEventResponse] -> rw_RatioToNLO = " << rw_RatioToNLO << std::endl;
  std::cout << "[AMUVariation::GetEventResponse] -> onesig_RatioToNLO = " << onesig_RatioToNLO << std::endl;
*/
  resp.push_back( {hdr.systParamId, {}} );
  for (double var : hdr.paramVariations) {
    double this_reweight = 1.0 + var * onesig_RatioToNLO;
    resp.back().responses.push_back( this_reweight );
  }

  return resp;
}

std::string AMUVariation::AsString() { return "AMUVariation"; }

AMUVariation::~AMUVariation() {
} 
