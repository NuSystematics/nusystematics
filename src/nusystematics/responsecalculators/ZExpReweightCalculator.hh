#ifndef nusystematics_RESPONSE_CALCULATORS_ZExpReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_ZExpReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_ZExpReweightCalculator_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_ZExpReweightCalculator_FILEPATH);

namespace nusyst {

  class ZExpReweightCalculator{

  protected:

    std::unique_ptr<TH2D> map_reweight;

  public:

    ZExpReweightCalculator(fhicl::ParameterSet const &InputManifest) : doDebug(false) {
      LoadInputHistograms(InputManifest);
    }
    ~ZExpReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetZExpReweight(double Enu, double Q2, double parameter_value);

    std::string GetCalculatorName() const { return "ZExpReweightCalculator"; }

    bool doDebug;

  };

  inline double ZExpReweightCalculator::GetZExpReweight(double Enu, double Q2, double parameter_value){

    if(doDebug){
      printf("[ZExpReweightCalculator::GetZExpReweight] (Enu, Q2) = (%1.2f, %1.2f)\n", Enu, Q2);
    }

    // DEFINE HOW WE GET THE REWEIGHT
    // Example below for ROOT Interpolate function
    double this_rw = map_reweight->Interpolate(Enu, Q2);

    if(doDebug){
      printf("[ZExpReweightCalculator::GetZExpReweight] this_rw = %1.2f\n", this_rw);
    }

    // parameter_value = 0 gives CV
    // paremeter_value = 1 gives new Zexp
    // other than that, linerly interpolate

    double weight = 1. * (1. - parameter_value ) + this_rw * parameter_value;



    return weight;

  }

  inline void ZExpReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string input_file = ps.get<std::string>("input_file");
    if(input_file.find("/")!=0){
      std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
      if(tmp_NUSYSTEMATICS_ROOT==""){
        throw invalid_ZExpReweightCalculator_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
      }
      input_file = tmp_NUSYSTEMATICS_ROOT+"/data/"+input_file;
    }
    std::string hist_name = ps.get<std::string>("hist_name");

    map_reweight = std::unique_ptr<TH2D>( GetHistogram<TH2D>(input_file, hist_name) );

    doDebug = ps.get<bool>("debug", false);

  }


} // namespace nusyst

#endif
