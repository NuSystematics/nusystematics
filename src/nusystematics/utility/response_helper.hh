#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "nusystematics/utility/make_instance.hh"

#include "systematicstools/interface/SystParamHeader.hh"

#include "systematicstools/interpreters/ParamHeaderHelper.hh"

#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"

#include "Framework/EventGen/EventRecord.h"

#include "TFile.h"
#include "TTree.h"

#include <chrono>
#include <string>

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(response_helper_found_no_parameters);

class response_helper : public systtools::ParamHeaderHelper {

  size_t NEvsProcessed;
  std::map<simb_mode_copy, std::map<size_t, std::tuple<double, double, size_t>>>
      ProfileStats;

private:
  constexpr static size_t Order = 5;
  constexpr static size_t NCoeffs = Order + 1;

  size_t ProfilerRate;

  std::string config_file;
  std::vector<std::unique_ptr<IGENIESystProvider_tool>> syst_providers;

public:
  response_helper() : NEvsProcessed(0), ProfilerRate(0) {}
  response_helper(std::string const &fhicl_config_filename) : NEvsProcessed(0) {
    LoadConfiguration(fhicl_config_filename);
  }

  std::vector<std::unique_ptr<IGENIESystProvider_tool>>& GetSystProvider(){
    return syst_providers;
  };

  void LoadProvidersAndHeaders(fhicl::ParameterSet const &ps) {
    syst_providers = systtools::ConfigureISystProvidersFromParameterHeaders<
        IGENIESystProvider_tool>(ps, make_instance);
    
    if (!syst_providers.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected to load some systematic providers from input: "
          << std::quoted(config_file);
    }
    
    systtools::param_header_map_t configuredParameterHeaders =
        systtools::BuildParameterHeaders(syst_providers);
    if (!configuredParameterHeaders.size()) {
      throw response_helper_found_no_parameters()
          << "[ERROR]: Expected systematric providers loaded from input: "
          << std::quoted(config_file) << " to provide some parameter headers.";
    }
    
    SetHeaders(configuredParameterHeaders);
  }

  void LoadConfiguration(std::string const &fhicl_config_filename) {
    config_file = fhicl_config_filename;

    // TODO
    std::unique_ptr<cet::filepath_maker> fm = std::make_unique<cet::filepath_maker>();
    fhicl::ParameterSet ps = fhicl::ParameterSet::make(config_file, *fm);

    LoadProvidersAndHeaders(ps.get<fhicl::ParameterSet>(
        "generated_systematic_provider_configuration"));

    ProfilerRate = ps.get<size_t>("ProfileRate", 0);
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  systtools::event_unit_response_t
  GetEventResponses(genie::EventRecord const &GenieGHep,
                    systtools::paramId_t i) {
    systtools::event_unit_response_t response;
    for (auto &sp : syst_providers) {
      systtools::event_unit_response_t prov_response =
          sp->GetEventResponse(GenieGHep, i);
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }
    return response;
  }

  systtools::event_unit_response_w_cv_t
  GetEventVariationAndCVResponse(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_w_cv_t response;

    simb_mode_copy mode = GetSimbMode(GenieGHep);

    for (size_t sp_it = 0; sp_it < syst_providers.size(); ++sp_it) {
      std::unique_ptr<IGENIESystProvider_tool> const &sp =
          syst_providers[sp_it];

      std::chrono::high_resolution_clock::time_point start;
      if (ProfilerRate) {
        start = std::chrono::high_resolution_clock::now();
      }

      systtools::event_unit_response_w_cv_t prov_response =
          sp->GetEventVariationAndCVResponse(GenieGHep);

      if (ProfilerRate && prov_response.size()) {
        auto end = std::chrono::high_resolution_clock::now();
        auto diff_ms =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        if (!ProfileStats[mode].count(sp_it)) {
          ProfileStats[mode][sp_it] = {0, 0, 0};
        }

        std::get<2>(ProfileStats[mode][sp_it])++;

        double delta = diff_ms - std::get<0>(ProfileStats[mode][sp_it]);

        std::get<0>(ProfileStats[mode][sp_it]) +=
            delta / std::get<2>(ProfileStats[mode][sp_it]);

        std::get<1>(ProfileStats[mode][sp_it]) +=
            delta * (diff_ms - std::get<0>(ProfileStats[mode][sp_it]));
      }
      for (auto &&er : prov_response) {
        response.push_back(std::move(er));
      }
    }

    if (ProfilerRate && NEvsProcessed && !(NEvsProcessed % ProfilerRate)) {
      std::cout << std::endl
                << "[PROFILE]: Event number = " << NEvsProcessed << std::endl;
      for (size_t sp_it = 0; sp_it < syst_providers.size(); ++sp_it) {
        std::cout << "\tSystProvider: "
                  << syst_providers[sp_it]->GetFullyQualifiedName()
                  << std::endl;
        for (auto const &m : ProfileStats) {
          if (!m.second.count(sp_it) || (std::get<2>(m.second.at(sp_it)) < 2)) {
            continue;
          }
          std::cout << "\t\tMode: " << tostr(m.first)
                    << ", NEvs = " << std::get<2>(m.second.at(sp_it))
                    << ", mean: " << std::get<0>(m.second.at(sp_it))
                    << " us, stddev: "
                    << sqrt(std::get<1>(m.second.at(sp_it)) /
                            (std::get<2>(m.second.at(sp_it)) - 1))
                    << " us." << std::endl;
        }
      }
    }
    NEvsProcessed++;

    return response;
  }

  double GetEventWeightResponse(genie::EventRecord const &GenieGHep,
                                systtools::param_value_list_t const &vals) {
    double weight = 1;
    for (auto &sp : syst_providers) {
      weight *= sp->GetEventWeightResponse(GenieGHep, vals);
    }
    return weight;
  }
}; // namespace nusyst
} // namespace nusyst
