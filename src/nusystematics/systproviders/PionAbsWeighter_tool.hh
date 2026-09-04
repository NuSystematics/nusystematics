#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"
#include "systematicstools/utility/YAMLSystParamHeaderUtility.hh"

//#include "RwCalculators/GReWeightNuXSecCCQE.h"
#include "nusystematics/responsecalculators/PionAbsResponseCalculator.hh"

#include "Framework/Messenger/Messenger.h"

#include "TFile.h"
#include "TTree.h"

#include <array>
#include <memory>
#include <string>

#include <cmath>

/*
 * \class PionAbsWeighter
 *
 * \brief A weighting tool that applies the ratio of final state nucleon distributions in pion absorption events from GENIE hA to INCL or G4
 *        --> Numerator: final state nucleon multiplicity in pion abs event from INCL / G4
 *        --> Denominator: final state nucleon multiplicity in pion abs event from GENIE hA
 *
 * \created November 11, 2025
 *
 * \authors Matt King <mking9@uchicago.edu>
 *          John Plows <kplows@liverpool.ac.uk>
 *          Gray Putnam <putnam@fnal.gov>
 */

class PionAbsWeighter : public nusyst::IGENIESystProvider_tool {

  std::unique_ptr<nusyst::PionAbsResponseCalculator> PionAbsCalculator;

public:
  explicit PionAbsWeighter(YAML::Node const &);
  ~PionAbsWeighter();
  //'First' configuration step: tool configuration
  // - takes 'arbitrary' YAML configuration and
  //   produces SystMetaData object which can later be used to configure a
  //   specific set of parameter values to be calculated
  systtools::SystMetaData BuildSystMetaData(YAML::Node const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  // - Reads the preconstructed SystMetaData produced by BuildSystMetaData
  //   to configure an instance of this class to calculate weights
  // - Recieves a copy of the tool_options instance constructed by
  //   BuildSystMetaData as an argument
  bool SetupResponseCalculator(YAML::Node const &);

  // Used to pass arbitrary YAML options from the tool configuration to the
  //   parameter headers.
  YAML::Node GetExtraToolOptions() { return tool_options; }

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);
  
  // Can add as much or as little stateful information here for use when
  // representing this instance as a string.
  std::string AsString() { return "PionAbsWeighter"; }

private:
  // arbitrary additional configuration from the tool configuration/parameter
  // headers can be storeds here
  YAML::Node tool_options;
  
  // Vector to hold the name of the descriptors.
  // Will be loaded into the appropriate members in BuildSystMetaData()
  std::vector<std::string> descriptors = {
    "MultiplicityDiffG4", /// Reweight based on the ratio of nucleon difference (p - n) of G4 to hA
    "MultiplicityDiffINCL", /// Reweight based on the ratio of nucleon difference (p - n) of INCL to hA
    "QuasiDeuteronFraction" /// Reweight based on uncertainty in the absorption mechanism (quasi-deuteron vs mulit-nucleon) for low-KE pion abs
  };

  //size_t ResponseParameterIdx;
  std::vector<size_t> ResponseParameterIndices;

  // The ParamHeaders id of the four free parameters provided by this
  // systprovider
  //std::array<size_t, 4> pidx_Params;

  // The configured variations to precalculate for the four parameters
  std::array<double, 4> CVs;
  std::array<std::vector<double>, 4> Variations;

  // Concrete example is more clear with some actual implementation, we will use
  // some GENIE reweight engines.

  // configurable verbosity as an example of some arbitrary systprovider
  // configuration
  int verbosity_level;

  //Helper functions
  std::vector<int> IndicesOfAbsNucleons(genie::EventRecord const &, int);
  int DiffAbsNucAboveThresh(genie::EventRecord const &, int, double);
  std::string LabelBase(genie::EventRecord const &, int);
  double ReweighthADiff(genie::EventRecord const &, int, double,bool);
  double ReweightQDFraction(genie::EventRecord const &, int, double);
  double ScalableErf(double, double, double);
}; //class PionAbsWeighter
