# `nusystematics`

Implements neutrino interaction systematics for GENIE3 events (including an 
interface to GReWeight) in the `systematicstools` framework.

## To Build

`nusystematics` will by default build `systematicstools` for you, but requires 
ROOT v6+ and GENIE v3+ to be set up in the environment at the time of 
configuration. Requires cmake 3.21+ (see [here](#fetch-cmake-321) for updating CMake hints).

Having checked out the repository in `/path/to/repo`:

```bash
cd /path/to/repo
mkdir build; cd build
cmake ..
make install -j 8
```

### What if you want to use your own copy of systematicstools

If you want to use a source distribution of `systematicstools` that you are 
modifying, then configure with:

```bash
cmake .. -DCPM_systematicstools_SOURCE=/path/to/systematicstools
```

if you want to use a binary distribution that you have built elsewhere, you need
to make sure that `systematicstools_ROOT` is defined in the environment, and 
then it should be picked up automatically by `CPMFindPackage`. Check the 
configuration output for a line like:

```bash
-- CPM: using local package systematicstools@(...)
```

## Adding a new systprovider

First off, read documentation describing systematicstools design choices 
[here](https://github.com/jedori0228/systematicstools#introduction).

Configuration is 'two step', first from a [_tool configuration_](https://github.com/jedori0228/systematicstools/blob/develop/doc/ToolConfiguration.md#tool-configuration) 
FHiCL file, which produces the [_parameter headers_](https://github.com/jedori0228/systematicstools/blob/develop/doc/ParameterHeaders.md) 
FHiCL file that is used for the second step. The parameter headers file is then 
used to configure an instance of your systprovider into a state where it can 
produce event responses. For more high-level information on writing a provider, 
see [here](https://github.com/jedori0228/systematicstools/blob/develop/doc/WritingAProvider.md).

The rest of this section will give short, documented examples of a skeleton 
systprovider implementation. Full implementation can be found in [SkeleWeighter_tool.hh](src/nusystematics/systproviders/SkeleWeighter_tool.hh) and [SkeleWeighter_tool.cc](src/nusystematics/systproviders/SkeleWeighter_tool.cc)


### Header

```c++
#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "RwCalculators/GReWeightNuXSecCCQE.h"

#include "TFile.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <string>

class SkeleWeighter : public nusyst::IGENIESystProvider_tool {
public:
  explicit SkeleWeighter(fhicl::ParameterSet const &);

  //'First' configuration step: tool configuration
  // - takes 'arbitrary' FHiCL configuration and
  //   produces SystMetaData object which can later be used to configure a
  //   specific set of parameter values to be calculated
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  // - Reads the preconstructed SystMetaData produced by BuildSystMetaData
  //   to configure an instance of this class to calculate weights
  // - Recieves a copy of the tool_options instance constructed by 
  //   BuildSystMetaData as an argument
  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  // Used to pass arbitrary FHiCL options from the tool configuration to the
  //   parameter headers.
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  // Can add as much or as little stateful information here for use when
  // representing this instance as a string.
  std::string AsString() { return "SkeleWeighter"; }

  ~SkeleWeighter(){}

private:
  // arbitrary additional configuration from the tool configuration/parameter
  // headers can be storeds here
  fhicl::ParameterSet tool_options;

  // The ParamHeaders id of the four free parameters provided by this
  // systprovider
  std::array<size_t, 4> pidx_Params;

  // The configured variations to precalculate for the four parameters
  std::array<double, 4> CVs;
  std::array<std::vector<double>, 4> Variations;

  // Concrete example is more clear with some actual implementation, we will use
  // some GENIE reweight engines.
  std::array<std::vector<genie::rew::GReWeightNuXSecCCQE>, 4> ReWeightEngines;

  // configurable verbosity as an example of some arbitrary systprovider
  // configuration
  int verbosity_level;
};
```

### Tool Configuration

```c++
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
```

### Parameter Headers Configuration

```c++
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
```

### Event Responses

```c++
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

```

### Generating the Parameter Headers File

An example Tool Configuration for this systprovider can be found in [SkeleWeighgter.ToolConfig.fcl](fcl/SkeleWeighgter.ToolConfig.fcl),
it contains:

```
 MySkeleProvider_toolconfig: {
    tool_type: "SkeleWeighter"
    instance_name: "myskele"

    p1_variation_descriptor: "(-3,3,1)"
    p2_variation_descriptor: "(-3,3,1)"
    p3_variation_descriptor: "(-3,3,1)"
    p4_variation_descriptor: "(-3,3,1)"

    verbosity_level: 5
  }

syst_providers: [MySkeleProvider_toolconfig]
```

It can be processed to a parameter header file like: `GenerateSystProviderConfigNuSyst -c ./fcl/SkeleWeighter.ToolConfig.fcl` which produces the below output:

```
BEGIN_PROLOG
generated_systematic_provider_configuration: { 
	SkeleWeighter_myskele: { 
		instance_name: myskele 

		p1: { isSplineable: true paramVariations: [-3, -2, -1, 0, 1, 2, 3] prettyName: p1 systParamId: 0 } 
		p2: { isSplineable: true paramVariations: [-3, -2, -1, 0, 1, 2, 3] prettyName: p2 systParamId: 1 } 
		p3: { isSplineable: true paramVariations: [-3, -2, -1, 0, 1, 2, 3] prettyName: p3 systParamId: 2 } 
		p4: { isSplineable: true paramVariations: [-3, -2, -1, 0, 1, 2, 3] prettyName: p4 systParamId: 3 } 
		parameter_headers: [p1, p2, p3, p4] 
		tool_options: { verbosity_level: 5 } 

		tool_type: SkeleWeighter 
	} 
	syst_providers: [SkeleWeighter_myskele] 
}
END_PROLOG
```

### Fetch CMake 3.21

If you just want a local build of cmake for this package and do not want to update
your system cmake, you can try the below recipe.

```
CMAKEVERSION=3.21.7
cd /path/to/repo
wget https://github.com/Kitware/CMake/archive/refs/tags/v${CMAKEVERSION}.tar.gz
tar -zxf v${CMAKEVERSION}.tar.gz
cd CMake-${CMAKEVERSION}
./bootstrap
make -j 8
export PATH=$(pwd)/bin:${PATH}
which cmake; cmake --version
```