#pragma once

/*
 * \class PionAbsResponseCalculator
 *
 * \brief Looks up from a histogram in (KEpi,multiplicity_diff) space the bin for an event, and applies it.
 *        Based on QEInterferenceResponseCalculator.hh for boilerplate code
 *
 * \created November 11, 2025
 *
 * \authors Matt King <mking9@uchicago.edu>
 *          John Plows <kplows@liverpool.ac.uk>
 *          Gray Putnam <putnam@fnal.gov>
 */

#include <algorithm>
#include <string>
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/string_parsers.hh"
#include "systematicstools/utility/exceptions.hh"
#include "yaml-cpp/yaml.h"
#include "TH1D.h"
#include "TH3D.h"

NEW_SYSTTOOLS_EXCEPT(invalid_PionAbs_INPUTHIST);
NEW_SYSTTOOLS_EXCEPT(invalid_PionAbs_FILEPATH);
NEW_SYSTTOOLS_EXCEPT(invalid_PionAbs_BIN);

namespace nusyst {

  class PionAbsResponseCalculator {

  public:
    PionAbsResponseCalculator(YAML::Node const & InputParams) {
      LoadInputHistograms(InputParams);
    }

    ~PionAbsResponseCalculator() {}

    void LoadInputHistograms(YAML::Node const & InputParams);

    double GetWeight(bool is_INCL, bool doSum, int pion_pdg, double KEpi_GeV, int mult);

    std::string GetCalculatorName() const { return "PionAbsResponseCalculator"; }

  private:
  
    std::vector<TH2D> ratio_histogram_collection;

  }; // class QEinterferenceResponseCalculator

  inline double PionAbsResponseCalculator::GetWeight( bool is_INCL, bool doSum, int pion_pdg,
							     double KEpi_GeV, int mult) {

    bool debug = false;

    if (debug) std::cout<<"In GetWeight for the event"<<std::endl;

    // Check if the pion pdg is handled, return 1.0 if not
    std::string pion_string;
    if (pion_pdg == 211) pion_string="PiPlus";
    else if (pion_pdg == -211) pion_string="PiMinus";
    else if (pion_pdg == 111)  pion_string="Pi0";
    else return 1.0;

    if (debug) std::cout<<"Pion is good."<<std::endl;

    // If outside the kinematic range, return 1.0
    // We cover pion KEs from 0-1 GeV
    if (KEpi_GeV > 1.0) return 1.0;

    if (debug) std::cout<<"Pion KE is good: "<<KEpi_GeV<<std::endl;
    //Grab the correct histogram
    //Vector of 6 histograms -- order in fcl:
    //Pi Plus G4
    //Pi Plus INCL
    //Pi Minus G4
    //Pi Minus INCL
    //Pi 0 G4
    //pi 0 INCL
    std::string gen_string;
    if (is_INCL) gen_string = "INCL";
    else gen_string = "G4";

   std::string diff_string;
   if (doSum) diff_string = "sum";
   else diff_string = "diff";

   std::string histo_name = Form("%s_%sReweight%s",pion_string.c_str(),gen_string.c_str(),diff_string.c_str());

   if (debug) std::cout<<"Looking for histo named "<<histo_name<<std::endl;

   TH2D this_ratio_histogram;
   //loop, GetTitle
   for (const TH2D& ratio_histogram : ratio_histogram_collection) {
     if (debug) std::cout<<"Current histo title: ratio_histogram.GetName(): "<<ratio_histogram.GetName()<<std::endl;
     if (ratio_histogram.GetName() == histo_name) {this_ratio_histogram = ratio_histogram; break;}
    }
   // Throw a fit if names don't match.
   if (debug) std::cout<<"Title of ratio_histogram: "<<this_ratio_histogram.GetName()<<std::endl;
   
   //Calculate the correct bin index for the energy and multiplicity

   //round KEpi_GeV to the nearest 0.1 GeV. Anything below 0.1 maps to 0.1.
    std::string keStr;
    float rounded_KE = std::round(KEpi_GeV * 10.0f) / 10.0f;

    // Clamp to [0.1, 1.0]
    if (rounded_KE < 0.1f) rounded_KE = 0.1f;
    if (rounded_KE > 1.0f) rounded_KE = 1.0f;

    if (std::fabs(rounded_KE - 1.0) < 1e-6) {
       keStr = "1.0";
       } else {
       keStr = Form("%g", rounded_KE);
       } 

    if (debug) std::cout<<"rounded_KE: "<<rounded_KE<<std::endl;
    int x_bin_idx = static_cast<int>(std::round(10.0 * rounded_KE - 1.0));

    //Get multiplicity bin
    int nBinsMult  = doSum ? 60 : 120;
    double yMinMult = doSum ? 0.0 : -60.0;
    double yMaxMult = 60.0;

    double yBinWidth = (yMaxMult - yMinMult) / nBinsMult;

    int y_bin_idx = static_cast<int>(std::floor((mult - yMinMult) / yBinWidth));

    // Clamp to underflow/overflow range if needed
    if (y_bin_idx < 0) y_bin_idx = 0; 
    if (y_bin_idx > nBinsMult + 1) y_bin_idx = nBinsMult + 1;

    // Seek out the bin and return it
    if (debug) std::cout<<"x_bin_idx: "<<x_bin_idx<<"; y_bin_idx: "<<y_bin_idx<<std::endl;
    double weight = this_ratio_histogram.GetBinContent(x_bin_idx+1, y_bin_idx+1);

    if (debug) std::cout<<"Directly calling the weight bin content: "<<weight<<std::endl;

    weight = std::isnan(weight) ? 1.0 : weight; // Guards against NaNs in histogram

    if (debug) std::cout<<"Weight_to_return from GetWeight after killing NaNs: "<<weight<<std::endl;

    return (weight != 0.0) ? weight : 1.0; // default response if zero

  } // GetWeight()

  inline void PionAbsResponseCalculator::LoadInputHistograms(YAML::Node const & config)
  {
    //Grab the input file
    std::string default_input_file;
    if (config["input_file"]) {
      default_input_file = config["input_file"].as<std::string>();
    }

    // Obtain each reweight histogram
    for( const YAML::Node & val_config : config["inputs"] ) {
      std::string input_hist   = val_config["input_hist"].as<std::string>();
      std::string input_file   = default_input_file;
      if (val_config["input_file"]) {
        input_file = val_config["input_file"].as<std::string>();
      }

      input_file = systtools::expand_env_vars(input_file);

      TH2D ratio_histogram;
      try {
	ratio_histogram = *(GetHistogram<TH2D>(input_file, input_hist ));
      }
      catch (...) {
	throw invalid_PionAbs_INPUTHIST() << "[ERROR]: Could not load histograms for: " << input_hist;
      }
      ratio_histogram_collection.push_back( ratio_histogram );

    } // for each input histogram
  } // LoadInputHistograms

} // namespace nusyst
