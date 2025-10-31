#pragma once

/*
 * \class QEInterferenceResponseCalculator
 *
 * \brief Looks up from a histogram in (Enu, Q0, Q3) space the bin for an event, and applies it.
 *        Based on CCQERPAReweightCalculator.hh for boilerplate code
 *
 * \created September 22, 2025
 *
 * \authors John Plows <kplows@liverpool.ac.uk>
 *          Gray Putnam <putnam@fnal.gov>
 */

#include <algorithm>
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"
#include "fhiclcpp/ParameterSet.h"
#include "TH1D.h"
#include "TH3D.h"

NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_INPUTHIST);
NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_FILEPATH);
NEW_SYSTTOOLS_EXCEPT(invalid_QEIntf_BIN);

namespace nusyst {

  class QEInterferenceResponseCalculator {

  public: 
    QEInterferenceResponseCalculator(fhicl::ParameterSet const & InputParams) {
      LoadInputHistograms(InputParams);
    }

    ~QEInterferenceResponseCalculator() {}

    void LoadInputHistograms(fhicl::ParameterSet const & InputParams);

    double GetWeight(int pdg, int current, double Enu_GeV, double Q0_GeV, double Q3_GeV);

    std::string GetCalculatorName() const { return "QEInterferenceResponseCalculator"; }

  private:

    // The bins in each variable are in principle different. Fill these out per flavour
    std::vector<std::pair<int, std::vector<double>>> Enu_bin_collection;
    std::vector<std::pair<int, std::vector<double>>>  Q0_bin_collection;
    std::vector<std::pair<int, std::vector<double>>>  Q3_bin_collection;
    std::vector<std::pair<int, TH3D>> ratio_histogram_collection;
    
  }; // class QEinterferenceResponseCalculator

  inline double QEInterferenceResponseCalculator::GetWeight( int pdg, int current,
							     double Enu_GeV,
							     double Q0_GeV, double Q3_GeV ) {

    // Check if the flavour is handled, return 1.0 if not
    std::vector<std::pair<int, TH3D>>::iterator it_hst = ratio_histogram_collection.begin();
    int code = 10 * std::abs(pdg) + current;
    code = (pdg < 0) ? -1 * code : code;
    for( it_hst; it_hst != ratio_histogram_collection.end(); ++it_hst ) {
      if( (*it_hst).first == code ) break;
    }
    if( it_hst == ratio_histogram_collection.end() ) {
      return 1.0;
    }

    int collection_index = it_hst - ratio_histogram_collection.begin();
    TH3D ratio_histogram = ratio_histogram_collection[collection_index].second;
    std::vector<double> Enu_bins = Enu_bin_collection[collection_index].second;
    std::vector<double>  Q0_bins =  Q0_bin_collection[collection_index].second;
    std::vector<double>  Q3_bins =  Q3_bin_collection[collection_index].second;
    
    // If outside the kinematic range, return 1.0
    if( Enu_GeV < *(Enu_bins.begin()) || Enu_GeV > *(Enu_bins.end()-2) ||
	Q0_GeV < *(Q0_bins.begin()) || Q0_GeV > *(Q0_bins.end()-2) ||
	Q3_GeV < *(Q3_bins.begin()) || Q3_GeV > *(Q3_bins.end()-2) ) {
      return 1.0;
    }

    // Seek out the bin and return it

    std::vector<double>::iterator it_enu = Enu_bins.begin();
    std::vector<double>::iterator it_q0  = Q0_bins.begin();
    std::vector<double>::iterator it_q3  = Q3_bins.begin();
    
    while( it_enu < Enu_bins.end() && Enu_GeV >= *it_enu ){ ++it_enu; } --it_enu;
    while( it_q0 < Q0_bins.end() && Q0_GeV >= (*it_q0) ){ ++it_q0; } --it_q0;
    while( it_q3 < Q3_bins.end() && Q3_GeV >= (*it_q3) ){ ++it_q3; } --it_q3;

    int bin_enu = it_enu - Enu_bins.begin();
    int bin_q0  = it_q0  -  Q0_bins.begin();
    int bin_q3  = it_q3  -  Q3_bins.begin();

    double weight = ratio_histogram.GetBinContent( bin_enu, bin_q0, bin_q3 );
    weight = std::isnan(weight) ? 1.0 : weight; // Guards against NaNs in histogram if outside ACHILLES coverage
    return (weight != 0.0) ? weight : 1.0; // default response if zero

  } // GetWeight()

  inline void QEInterferenceResponseCalculator::LoadInputHistograms(fhicl::ParameterSet const & ps)
  {
    const std::string & default_input_file = ps.get<std::string>("input_file", "");
    const std::vector<std::string> known_flavours = { "numu", "numubar", "nue", "nuebar" };
    const std::vector<std::string> known_currents = { "cc", "nc" };
    
    // Obtain for each flavour-current combination the name of the input histogram, and load it in.
    for( fhicl::ParameterSet const & val_config :
	   ps.get<std::vector<fhicl::ParameterSet>>("inputs") ) {
      std::string config_name = val_config.get<std::string>("name"); // nu(mu, e)(bar)_(cc, nc)
      std::string input_hist   = val_config.get<std::string>("input_hist"); // name of the hist in file
      std::string input_enu_bins = val_config.get<std::string>("input_enu_bins");
      std::string input_q0_bins = val_config.get<std::string>("input_q0_bins");
      std::string input_q3_bins = val_config.get<std::string>("input_q3_bins");
      std::string input_file   = val_config.get<std::string>("input_file", default_input_file);

      std::string flavour_name = config_name.substr( 0, config_name.find("_") );
      std::string current_name = config_name.substr( config_name.find("_")+1 );

      if( std::find( known_flavours.begin(), known_flavours.end(), flavour_name.c_str() ) ==
	  known_flavours.end() ) {
	throw invalid_QEIntf_INPUTHIST() << "[ERROR]: Unknown input flavour " << flavour_name.c_str()
					 << " ; check inputs:name " << config_name.c_str();
      }
      if( std::find( known_currents.begin(), known_currents.end(), current_name.c_str() ) ==
	  known_currents.end() ) {
	throw invalid_QEIntf_INPUTHIST() << "[ERROR]: Unknown input current " << current_name.c_str()
					 << " ; check inputs:name " << config_name.c_str();
      }

      // If input_file is not given as an absolute path seatch ${NUSYSTEMATICS_FQ_DIR}/data/
      if( input_file.find("/") != 0 ) {
	if( std::getenv("nusystematics_ROOT") == "" ) {
	  throw invalid_QEIntf_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set!\n"
					  << "Given path: " << input_file << "\n"
					  << "Expect absolute path (starts with '/'), or file relative to ${nusystematics_ROOT}/data/";
	}
	input_file = std::string(std::getenv("nusystematics_ROOT")) + "/data/" + input_file;
      } // absolute path not given

      int nu_pdg = -1;
      int current = -1;
      if( flavour_name == "numu" ) nu_pdg = 14;
      else if( flavour_name == "numubar" ) nu_pdg = -14;
      else if( flavour_name == "nue" ) nu_pdg = 12;
      else if( flavour_name == "nuebar" ) nu_pdg = -12;
      if ( current_name == "cc" || current_name == "CC" ) current = 0;
      else if ( current_name == "nc" || current_name == "NC" ) current = 1;

      int collection_code = 10 * std::abs(nu_pdg) + current;
      collection_code = (nu_pdg < 0) ? -1 * collection_code : collection_code;

      std::vector<double> Enu_bins, Q0_bins, Q3_bins;

      TH3D ratio_histogram;
      TH1D hEnu_bins, hQ0_bins, hQ3_bins;
      try {
	ratio_histogram = *(GetHistogram<TH3D>( input_file, input_hist ));
	hEnu_bins = *(GetHistogram<TH1D>( input_file, input_enu_bins ));
	hQ0_bins  = *(GetHistogram<TH1D>( input_file, input_q0_bins  ));
	hQ3_bins  = *(GetHistogram<TH1D>( input_file, input_q3_bins  ));
      }
      catch (...) {
	throw invalid_QEIntf_INPUTHIST() << "[ERROR]: Could not load histograms for: "
					 << input_hist << ", " << input_enu_bins
					 << ", " << input_q0_bins << ", " << input_q3_bins;
      }

      for( Int_t ib = 1; ib <= hEnu_bins.GetNbinsX()-1; ib++ )
	Enu_bins.emplace_back( hEnu_bins.GetBinLowEdge(ib) );
      for( Int_t ib = 1; ib <= hQ0_bins.GetNbinsX()-1; ib++ ) 
	Q0_bins.emplace_back( hQ0_bins.GetBinLowEdge(ib) );
      for( Int_t ib = 1; ib <= hQ3_bins.GetNbinsX()-1; ib++ ) 
	Q3_bins.emplace_back( hQ3_bins.GetBinLowEdge(ib) );

      // Add the high edges
      Enu_bins.emplace_back( hEnu_bins.GetBinLowEdge(hEnu_bins.GetNbinsX()) + 
			     hEnu_bins.GetBinWidth(hEnu_bins.GetNbinsX()) );
      Q0_bins.emplace_back( hQ0_bins.GetBinLowEdge(hQ0_bins.GetNbinsX()) + 
			    hQ0_bins.GetBinWidth(hQ0_bins.GetNbinsX()) );
      Q3_bins.emplace_back( hQ3_bins.GetBinLowEdge(hQ3_bins.GetNbinsX()) + 
			    hQ3_bins.GetBinWidth(hQ3_bins.GetNbinsX()) );

      // Add these to the collection
      Enu_bin_collection.emplace_back(std::pair<int,
				      std::vector<double>>( collection_code, Enu_bins ));
      Q0_bin_collection.emplace_back(std::pair<int,
				     std::vector<double>>( collection_code, Q0_bins ));
      Q3_bin_collection.emplace_back(std::pair<int,
				     std::vector<double>>( collection_code, Q3_bins ));
      ratio_histogram_collection.emplace_back(std::pair<int, 
					      TH3D>( collection_code, ratio_histogram ));

    } // for each flavour in inputs
  } // LoadInputHistograms

} // namespace nusyst
