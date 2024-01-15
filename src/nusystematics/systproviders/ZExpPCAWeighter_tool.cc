#include "nusystematics/systproviders/ZExpPCAWeighter_tool.hh"

#include "RwFramework/GSyst.h"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "TRandom3.h"

#include "Eigen/Dense"

#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>

#define private public
#define protected public
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
using namespace genie;
using namespace genie::rew;
fstream file;

using namespace nusyst;
using namespace systtools;
using namespace std;

// global variables
namespace PRD_93_113015 {
Eigen::Vector4d a_14_cv{2.3, -0.6, -3.8, 2.3};
Eigen::Vector4d a_14_errors{0.13, 1, 2.5, 2.7};
Eigen::MatrixXd Covariance_Matrix{{0.0169, 0.0455, -0.22035, 0.214461},
                                  {0.0455, 1., -2.245, 0.9909},
                                  {-0.22035, -2.245, 6.25, -4.62375},
                                  {0.214461, 0.9909, -4.62375, 7.29}};
Eigen::Array4d central_values_a_from_genie = {2.30, -0.60, -3.80, 2.30};
Eigen::Array4d central_values_afrom_errors_genie_code = {0.14, 0.67, 1,
                                                         0.75}; //-----in genie
Eigen::Array4d central_values_afrom_errors_genie = {0.13, 1, 2.5, 2.7};
} // namespace PRD_93_113015

Eigen::MatrixXd
GetPCAFromCovarianceMatrix(Eigen::MatrixXd const &Covariance_Matrix) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Covariance_Matrix);
  if (eigensolver.info() != Eigen::Success)
    abort();
  Eigen::MatrixXd eigenvector_colmatrix = eigensolver.eigenvectors();
  Eigen::VectorXd eigenvalues_sqrt = (eigensolver.eigenvalues()).array().sqrt();
  for (size_t i = 0; i < 4; ++i) {
    eigenvector_colmatrix.col(i) *= eigenvalues_sqrt[i];
  }
  return eigenvector_colmatrix;
}

// Function that takes uncorellated b values, and transforms them into a values
// using PCA
Eigen::VectorXd ChangeBasisBParams(std::vector<double> const &BVals,
                                   Eigen::MatrixXd const &PCA) {
  Eigen::VectorXd avals = Eigen::VectorXd::Zero(4);
  for (size_t i = 0; i < 4; ++i) {
    avals += PCA.col(i) * BVals[i];
  }
  return avals;
}

// Function to rescale the a parameters for GENIE ReWeight
Eigen::Array4d ScaleAparamsforGenie(Eigen::VectorXd &avals) {
  Eigen::VectorXd avalues_for_genie =
      (((PRD_93_113015::a_14_cv + avals).array() /
        PRD_93_113015::central_values_a_from_genie.array()) -
       1.0) /
      PRD_93_113015::central_values_afrom_errors_genie_code.array();
  return avalues_for_genie;
}

std::array<std::string, 4> b_params_for_zexp_names = {
    "b1", "b2", "b3",
    "b4"}; // the  b params to be rotated using PCA (in fcl file)

// constructor passes up configuration object to base class for generic tool
// initialization and initialises our local copies of paramIds to unconfigured
//  flag values
ZExpPCAWeighter::ZExpPCAWeighter(
    fhicl::ParameterSet const &params) // fn takes the parameters in fcl file
    : IGENIESystProvider_tool(
          params), // IGENIESystProvider_tool takes parameters and sets tune and
                   // event records etc. in GENIE
      pidx_Params{kParamUnhandled<size_t>, kParamUnhandled<size_t>,
                  kParamUnhandled<size_t>, kParamUnhandled<size_t>} {
} // The ParamHeaders id of the four free parameters provided by
  // thissystprovider

SystMetaData ZExpPCAWeighter::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                                paramId_t firstId) {
  SystMetaData smd;
  SystParamHeader responseParam;

  // loop through the four named parameters that this tool provides
  for (std::string const &pname : b_params_for_zexp_names) {
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
  tool_options.put("verbosity_level",
                   ps.get<int>("verbosity_level",
                               0)); // put tool config options in papam geaders

  return smd;
}

bool ZExpPCAWeighter::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  genie::RunOpt *grunopt = genie::RunOpt::Instance();
  grunopt->EnableBareXSecPreCalc(true);

  // loop through the named parameters (i) that this tool provides, check that
  // they are configured, and grab their id in the current systmetadata and set
  // up and pre-calculations/configurations required.
  for (size_t i = 0; i < b_params_for_zexp_names.size(); ++i) {

    if (!HasParam(md, b_params_for_zexp_names[i])) { // do we have the b params
                                                     // in the header file?
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter "
                  << b_params_for_zexp_names[i]
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }
    pidx_Params[i] = GetParamIndex(
        md,
        b_params_for_zexp_names[i]); // the id of the ith param is the paaram
                                     // index which is defined by headers+name

    if (verbosity_level > 1) {
      std::cout << "[INFO]: Have parameter "
                << b_params_for_zexp_names[i] // if params are present, start
                                              // config steps
                << " in SystMetaData with ParamId: " << pidx_Params[i]
                << ". Configuring." << std::endl;
    }

    auto param_md = md[pidx_Params[i]]; // param_md
    CVs[i] =
        param_md.centralParamValue; // the initial CV of parameter i (i=1,2,3,4)
    Variations[i] =
        param_md.paramVariations; // the variation of parameter i (i=1,2,3,4)

    // in the concrete version of this example we're going to configure a GENIE
    // GReWeightNuXSecCCQE instance for each configured variation.
    for (auto v : Variations[i]) { // loop through the variations of each param
                                   // e.g (-3,-2...,2,3)std in fcl file
      static std::array<genie::rew::GSyst_t, 4> const zexp_dials{
          genie::rew::kXSecTwkDial_ZExpA1CCQE,
          genie::rew::kXSecTwkDial_ZExpA2CCQE,
          genie::rew::kXSecTwkDial_ZExpA3CCQE,
          genie::rew::kXSecTwkDial_ZExpA4CCQE};

      // instantiate a new weight engine instance for this variation
      ReWeightEngines_new[i].push_back(
          std::make_unique<genie::rew::GReWeightNuXSecCCQE>());
      ReWeightEngines_new[i].back()->SetMode(GReWeightNuXSecCCQE::kModeZExp);

      // Set the b variations
      std::vector<double> bvariations = {0, 0, 0, 0};
      bvariations[i] = v;
      // Get the b parameter variations, using the variations in the parameter
      // md and the rotation function
      Eigen::VectorXd a_variations = ChangeBasisBParams(
          bvariations,
          GetPCAFromCovarianceMatrix(PRD_93_113015::Covariance_Matrix));

      // Now use the b_variations to get the a values
      auto myaparameters = ScaleAparamsforGenie(a_variations);

      if (verbosity_level > 2) {
        std::cout << "b[" << i << "] = " << v << std::endl;
      }
      // loop through the new a params and use them in GENIE
      for (size_t aval_j = 0; aval_j < 4; ++aval_j) {
        if (verbosity_level > 2) {
          std::cout << "  a[" << aval_j << "] param = " << a_variations[aval_j]
                    << ", scaled for GENIE = " << myaparameters[aval_j]
                    << std::endl;
        }
        ReWeightEngines_new[i].back()->SetSystematic(zexp_dials[aval_j],
                                                     myaparameters[aval_j]);
      }
      ReWeightEngines_new[i].back()->Reconfigure();
      if (verbosity_level > 2) {
        std::cout << "Done Reconfigure()" << std::endl;
      }

      if (verbosity_level > 2) {
        std::cout << "[LOUD]: Configured GReWeightNuXSecCCQE instance for "
                     "GENIE dial: "
                  << genie::rew::GSyst::AsString(zexp_dials[i])
                  << " at variation: " << v << std::endl;
      }
    }
  }
  return true; // returning cleanly
}

event_unit_response_t
ZExpPCAWeighter::GetEventResponse(genie::EventRecord const &ev) {

  event_unit_response_t resp;
  SystMetaData const &md = GetSystMetaData();

  // return early if this event isn't one we provide responses for
  if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
      !ev.Summary()->ProcInfo().IsWeakCC() ||
      ev.Summary()->ExclTag().IsCharmEvent()) {
    return resp;
  }

  // loop through and calculate weights
  for (size_t b_i = 0; b_i < b_params_for_zexp_names.size(); ++b_i) {
    // this parameter wasn't configured, nothing to do
    if (pidx_Params[b_i] == kParamUnhandled<size_t>) {
      continue;
    }

    // initialize the response array with this paramId
    resp.push_back({md[pidx_Params[b_i]].systParamId, {}});

    // loop through variations for this parameter
    for (size_t v_it = 0; v_it < Variations[b_i].size(); ++v_it) {

      // put the response weight for this variation of this parameter into the
      // response object
      resp.back().responses.push_back(
          ReWeightEngines_new[b_i][v_it]->CalcWeight(ev));

      if (verbosity_level > 3) {
        std::cout << "[DEBG]: For parameter " << md[pidx_Params[b_i]].prettyName
                  << "The  central value of the parameter is  " << CVs[b_i]
                  << " at variation[" << v_it
                  << "] std = " << Variations[b_i][v_it]
                  << " calculated weight: " << resp.back().responses.back()
                  << std::endl;
      }
    }
  }
  return resp;
}
