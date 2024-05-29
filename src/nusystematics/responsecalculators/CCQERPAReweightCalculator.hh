#ifndef nusystematics_RESPONSE_CALCULATORS_CCQERPAReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCQERPAReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_RPA_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_RPA_FILEPATH);

namespace nusyst {

  class CCQERPAReweightCalculator{

    enum ENuRange {
      LowE = 0,
      HighE = 1,
    };

  protected:

    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithRPAXSec;
    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithoutRPAXSec;

    std::map<int, double> x_FirstBinCenter, x_LastBinCenter;
    std::map<int, double> y_FirstBinCenter, y_LastBinCenter;
    std::map<int, double> z_FirstBinCenter, z_LastBinCenter;

    double ENuBoundary;

  public:

    CCQERPAReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~CCQERPAReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetRPAReweight(double Enu_GeV, std::array<double, 2> bin_kin, double parameter_value);

    std::string GetCalculatorName() const { return "CCQERPAReweightCalculator"; }

  };

  inline double CCQERPAReweightCalculator::GetRPAReweight(double Enu_GeV, std::array<double, 2> bin_kin, double parameter_value){
  //inline double CCQERPAReweightCalculator::GetRPAReweight(double Enu_GeV, double kin_Y, double kin_Z, double parameter_value){

    int enu_range = (Enu_GeV<ENuBoundary) ? 0 : 1;

    //printf("[CCQERPAReweightCalculator::GetRPAReweight] (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f), enu_range = %d\n", Enu_GeV, bin_kin[0], kin_Z, enu_range);

    static double Enu_GeV_epsil = 1E-6;
    double Enu_GeV_ForInterp = Enu_GeV;
    Enu_GeV_ForInterp = std::max( Enu_GeV_ForInterp, x_FirstBinCenter[enu_range] + Enu_GeV_epsil );
    Enu_GeV_ForInterp = std::min( Enu_GeV_ForInterp, x_LastBinCenter[enu_range] - Enu_GeV_epsil );

    static double kin_Y_epsil = 1E-6;
    double kin_Y_ForInterp = bin_kin[0];
    kin_Y_ForInterp = std::max( kin_Y_ForInterp, y_FirstBinCenter[enu_range] + kin_Y_epsil );
    kin_Y_ForInterp = std::min( kin_Y_ForInterp, y_LastBinCenter[enu_range] - kin_Y_epsil );

    static double kin_Z_epsil = 1E-6;
    double kin_Z_ForInterp = bin_kin[1];
    kin_Z_ForInterp = std::max( kin_Z_ForInterp, z_FirstBinCenter[enu_range] + kin_Z_epsil );
    kin_Z_ForInterp = std::min( kin_Z_ForInterp, z_LastBinCenter[enu_range] - kin_Z_epsil );

    //printf("[CCQERPAReweightCalculator::GetRPAReweight] -> (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);
    double xsec_WithRPA = map_ENuRange_to_WithRPAXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp); // CV
    double xsec_WithoutRPA = map_ENuRange_to_WithoutRPAXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);

    //printf("[CCQERPAReweightCalculator::GetRPAReweight] xsec (With RPA, Without RPA) = (%1.3f, %1.3e)\n", xsec_WithRPA, xsec_WithoutRPA);

    if(xsec_WithRPA==0.){
/*
      printf("[CCQERPAReweightCalculator::GetRPAReweight] Zero cross section for\n");
      printf("[CCQERPAReweightCalculator::GetRPAReweight] (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f), enu_range = %d\n", Enu_GeV, bin_kin[0], bin_kin[1], enu_range);
      printf("[CCQERPAReweightCalculator::GetRPAReweight] -> (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);
*/
      return 1.;
    }

    double weight = ( xsec_WithRPA * (1.-parameter_value) + xsec_WithoutRPA * parameter_value ) / xsec_WithRPA;

    //std::cout << "[CCQERPAReweightCalculator] weight = " << weight << std::endl;

    if(weight!=weight){

      printf("[CCQERPAReweightCalculator::GetRPAReweight] Nan weight for\n"); 
      printf("[CCQERPAReweightCalculator::GetRPAReweight] (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f), enu_range = %d\n", Enu_GeV, bin_kin[0], bin_kin[1], enu_range);
      printf("[CCQERPAReweightCalculator::GetRPAReweight] -> (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);
      weight = 1.;

    }

    return weight;

  }

  inline void CCQERPAReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = ps.get<std::string>("input_file", "");
    ENuBoundary = ps.get<double>("ENuBoundary");
    printf("[CCQERPAReweightCalculator::GetRPAReweight] ENuBoundary = %1.2f\n", ENuBoundary);

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file); // If specified per hist, replace it

      // if it does not start with "/", find it under ${NUSYSTEMATICS_FQ_DIR}/data/
      if(input_file.find("/")!=0){
        std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
        if(tmp_NUSYSTEMATICS_ROOT==""){
          throw invalid_CCQE_RPA_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
        }
        input_file = tmp_NUSYSTEMATICS_ROOT+"/data/"+input_file;
      }

      if(hName=="LowE_WithRPA"){
        map_ENuRange_to_WithRPAXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="LowE_WithoutRPA"){
        map_ENuRange_to_WithoutRPAXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="HighE_WithRPA"){
        map_ENuRange_to_WithRPAXSec[1] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="HighE_WithoutRPA"){
        map_ENuRange_to_WithoutRPAXSec[1] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }

    }

    for(int enu_range=0; enu_range<=1; enu_range++){

      const auto& h_WithRPAXsec = map_ENuRange_to_WithRPAXSec[enu_range];

      const auto& XAxis_WithRPAXsec = h_WithRPAXsec->GetXaxis();
      const auto& YAxis_WithRPAXsec = h_WithRPAXsec->GetYaxis();
      const auto& ZAxis_WithRPAXsec = h_WithRPAXsec->GetZaxis();

      x_FirstBinCenter[enu_range] = XAxis_WithRPAXsec->GetBinCenter(1);
      x_LastBinCenter[enu_range] = XAxis_WithRPAXsec->GetBinCenter( XAxis_WithRPAXsec->GetNbins() );

      y_FirstBinCenter[enu_range] = YAxis_WithRPAXsec->GetBinCenter(1);
      y_LastBinCenter[enu_range] = YAxis_WithRPAXsec->GetBinCenter( YAxis_WithRPAXsec->GetNbins() );

      z_FirstBinCenter[enu_range] = ZAxis_WithRPAXsec->GetBinCenter(1);
      z_LastBinCenter[enu_range] = ZAxis_WithRPAXsec->GetBinCenter( ZAxis_WithRPAXsec->GetNbins() );

      printf("@@ Enu range :%d\n", enu_range);
      printf("@@ - x-range: [%1.3f, %1.3f]\n", x_FirstBinCenter[enu_range], x_LastBinCenter[enu_range]);
      printf("@@ - y-range: [%1.3f, %1.3f]\n", y_FirstBinCenter[enu_range], y_LastBinCenter[enu_range]);
      printf("@@ - z-range: [%1.3f, %1.3f]\n", z_FirstBinCenter[enu_range], z_LastBinCenter[enu_range]);

    }

  }


} // namespace nusyst

#endif
