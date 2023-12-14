#pragma once

#include "systematicstools/interface/ISystProviderTool.hh"

#include "fhiclcpp/ParameterSet.h"

// GENIE
#include "Framework/EventGen/EventRecord.h"
// Extra includes needed for CheckTune()
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

namespace nusyst {

class IGENIESystProvider_tool : public systtools::ISystProviderTool {
protected:
  // Based on GENIEHelper::FindTune()
  // TODO: reduce code duplication here
  // -- S. Gardiner, 20 December 2018
  void CheckTune(const std::string &tune_name) {

    std::string fhicl_tune_name = tune_name;

    // The default tune name is ${GENIE_XSEC_TUNE}, which
    // should be converted into the value of the corresponding
    // enviornment variable, as is done below.
    if (fhicl_tune_name.front() == '$') {
      // need to remove ${}'s
      std::string tuneEnvVar = fhicl_tune_name;
      std::string rmchars("$(){} ");
      // std::remove_if removes characters in [first,last) that are found
      //   within the rmchars string. It returns returns a past-the-end
      //   iterator for the new end of the range [funky!]
      // std::string::erase actually trims the string
      tuneEnvVar.erase(std::remove_if(tuneEnvVar.begin(), tuneEnvVar.end(),
                                      [&rmchars](const char &c) -> bool {
                                        return rmchars.find(c) !=
                                               std::string::npos;
                                      }),
                       tuneEnvVar.end());

      const char *tune = std::getenv(tuneEnvVar.c_str());
      if (tune) {
        fhicl_tune_name = std::string(tune);
      } else {
        throw systtools::invalid_ToolConfigurationFHiCL()
            << "can't resolve TuneName: " << fhicl_tune_name;
      }
    }

    // If the XSecSplineList returns a non-empty string as the current tune
    // name, then genie::RunOpt::BuildTune() has already been called.
    std::string current_tune = genie::XSecSplineList::Instance()->CurrentTune();
    if (current_tune.empty()) {
      // Constructor automatically calls grunopt->Init();
      genie::RunOpt *grunopt = genie::RunOpt::Instance();
      grunopt->SetTuneName(fhicl_tune_name);
      grunopt->BuildTune();
    } else {
      // It has already been built, so just check consistency
      if (fhicl_tune_name != current_tune) {
        throw systtools::invalid_ToolConfigurationFHiCL()
            << "Requested GENIE tune \"" << fhicl_tune_name
            << "\" does not match previously built tune \"" << current_tune
            << '\"';
      }
    }
  }

public:
  IGENIESystProvider_tool(fhicl::ParameterSet const &ps)
      : ISystProviderTool(ps), fGENIEModuleLabel(ps.get<std::string>(
                                   "genie_module_label", "generator")) {

    std::string tune_name =
        ps.get<std::string>("TuneName", "${GENIE_XSEC_TUNE}");
    this->CheckTune(tune_name);
  }

  NEW_SYSTTOOLS_EXCEPT(invalid_response);

  /// Calculates configured response for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &) = 0;

  /// Calculates configured response for a given vector of GHep record
  std::unique_ptr<systtools::EventResponse>
  GetEventResponses(std::vector<std::unique_ptr<genie::EventRecord>> const &gheps){

    std::unique_ptr<systtools::EventResponse> er =
        std::make_unique<systtools::EventResponse>();

    for (size_t eu_it = 0; eu_it < gheps.size(); ++eu_it) {
      er->push_back(GetEventResponse(*gheps[eu_it]));
    }
    return er;

  };

  systtools::event_unit_response_w_cv_t
  GetEventVariationAndCVResponse(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_w_cv_t responseandCV;

    systtools::event_unit_response_t prov_response =
        GetEventResponse(GenieGHep);

    // Foreach param
    for (systtools::ParamResponses &pr : prov_response) {
      // Get CV resp
      systtools::SystParamHeader const &hdr =
          GetParam(GetSystMetaData(), pr.pid);

      if (pr.responses.size() != hdr.paramVariations.size()) {
        throw invalid_response()
            << "[ERROR]: Parameter: " << hdr.prettyName << ", with "
            << hdr.paramVariations.size() << " parameter variations, returned "
            << pr.responses.size() << " responses.";
      }

      double CVResp = hdr.isWeightSystematicVariation ? 1 : 0;
      size_t NVars = hdr.paramVariations.size();

      double cv_param_val = 0;
      if (hdr.centralParamValue != systtools::kDefaultDouble) {
        cv_param_val = hdr.centralParamValue;
      }
      for (size_t idx = 0; idx < NVars; ++idx) {
        if (fabs(cv_param_val - hdr.paramVariations[idx]) <=
            std::numeric_limits<float>::epsilon()) {
          CVResp = pr.responses[idx];
          break;
        }
      }


/*
//TODO
      // if we didn't find it, the CVResp stays as 1/0 depending on whether it
      // is a weight or not.
      for (size_t idx = 0; idx < NVars; ++idx) {
        if (hdr.isWeightSystematicVariation) {
          pr.responses[idx] /= CVResp;
        } else {
          pr.responses[idx] -= CVResp;
        }
      }
*/

      responseandCV.push_back({pr.pid, CVResp, pr.responses});
    } // end for parameter response

    return responseandCV;
  }

  /// Calculates the response to a single parameter for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &, systtools::paramId_t) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement systtools::event_unit_response_t "
           "GetEventResponse(genie::EventRecord &, systtools::paramId_t).";
  }

  /// Calculates the multiplicatively combined responses for a given set of
  /// parameter--value pairs.
  ///
  /// \note This convenience method should only be used for weight responses.
  virtual double GetEventWeightResponse(genie::EventRecord const &,
                                        systtools::param_value_list_t const &) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement double "
           "GetEventWeightResponse(genie::EventRecord "
           "&,systtools::param_value_list_t const &).";
  }

  std::string fGENIEModuleLabel;
};
} // namespace nusyst
