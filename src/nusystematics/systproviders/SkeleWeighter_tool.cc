#include "nusystematics/systproviders/SkeleWeighter_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "RwFramework/GSyst.h"

using namespace nusyst;
using namespace systtools;

std::array<std::string, 4> ParamPrettyNames = {"p1", "p2", "p3", "p4"};

// constructor passes up configuration object to base class for generic tool
// initialization and initialises our local copies of paramIds to unconfigured
//  flag values
SkeleWeighter::SkeleWeighter(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_Params{kParamUnhandled<size_t>, kParamUnhandled<size_t>,
                  kParamUnhandled<size_t>, kParamUnhandled<size_t>} {}

SystMetaData SkeleWeighter::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                              paramId_t firstId) {
  SystMetaData smd;
  SystParamHeader responseParam;

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
  }

  // Put any options that you want to propagate to the ParamHeaders options
  tool_options.put("verbosity_level", ps.get<int>("verbosity_level", 0));

  return smd;
}

bool SkeleWeighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  // loop through the named parameters that this tool provides, check that they
  // are configured, and grab their id in the current systmetadata and set up
  // and pre-calculations/configurations required.
  for (size_t i = 0; i < ParamPrettyNames.size(); ++i) {

    if (!HasParam(md, ParamPrettyNames[i])) {
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter " << ParamPrettyNames[i]
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }
    pidx_Params[i] = GetParamIndex(md, ParamPrettyNames[i]);

    if (verbosity_level > 1) {
      std::cout << "[INFO]: Have parameter " << ParamPrettyNames[i]
                << " in SystMetaData with ParamId: " << pidx_Params[i]
                << ". Configuring." << std::endl;
    }

    auto param_md = md[pidx_Params[i]];

    CVs[i] = param_md.centralParamValue;
    Variations[i] = param_md.paramVariations;

    // in the concrete version of this example we're going to configure a
    // GENIE GReWeightNuXSecCCQE instance for each configured variation.
    for (auto v : Variations[i]) {
      static std::array<genie::rew::GSyst_t, 4> const gdials{
          genie::rew::kXSecTwkDial_ZExpA1CCQE,
          genie::rew::kXSecTwkDial_ZExpA2CCQE,
          genie::rew::kXSecTwkDial_ZExpA3CCQE,
          genie::rew::kXSecTwkDial_ZExpA4CCQE};

      // instantiate a new weight engine instance for this variation
      ReWeightEngines[i].emplace_back();
      ReWeightEngines[i].back().SetMode(
          genie::rew::GReWeightNuXSecCCQE::kModeZExp);
      // set its systematic
      ReWeightEngines[i].back().SetSystematic(gdials[i], v);
      // configure it to weight events
      ReWeightEngines[i].back().Reconfigure();
      if (verbosity_level > 2) {
        std::cout << "[LOUD]: Configured GReWeightNuXSecCCQE instance for "
                     "GENIE dial: "
                  << genie::rew::GSyst::AsString(gdials[i])
                  << " at variation: " << v << std::endl;
      }
    }
  }
  // returning cleanly
  return true;
}

event_unit_response_t
SkeleWeighter::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  // return early if this event isn't one we provide responses for
  if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC() ||
      ev.Summary()->ExclTag().IsCharmEvent()) {
    return resp;
  }

  // loop through and calculate weights
  for (size_t i = 0; i < ParamPrettyNames.size(); ++i) {
    // this parameter wasn't configured, nothing to do
    if (pidx_Params[i] == kParamUnhandled<size_t>) {
      continue;
    }

    // initialize the response array with this paramId
    resp.push_back({md[pidx_Params[i]].systParamId, {}});

    // loop through variations for this parameter
    for (size_t v_it = 0; v_it < Variations[i].size(); ++v_it) {

      // put the response weight for this variation of this parameter into the
      // response object
      resp.back().responses.push_back(ReWeightEngines[i][v_it].CalcWeight(ev));
      if (verbosity_level > 3) {
        std::cout << "[DEBG]: For parameter " << md[pidx_Params[i]].prettyName
                  << " at variation[" << v_it << "] = " << Variations[i][v_it]
                  << " calculated weight: " << resp.back().responses.back()
                  << std::endl;
      }
    }
  }

  return resp;
}
