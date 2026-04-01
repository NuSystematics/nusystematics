#ifndef nusystematics_RESPONSE_CALCULATORS_CCQETemplateReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCQETemplateReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_Template_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_Template_FILEPATH);

namespace nusyst {

  // Utility function to find the closest bin center for a value
  inline double find_closest_bin_center(double val, const std::vector<double>& bin_centers) {
    double min_dist = std::abs(val - bin_centers[0]);
    double closest = bin_centers[0];
    for (const auto& bc : bin_centers) {
      double dist = std::abs(val - bc);
      if (dist < min_dist) {
        min_dist = dist;
        closest = bc;
      }
    }
    return closest;
  }

  class CCQETemplateReweightCalculator{

    enum ENuRange {
      LowE = 0,
      HighE = 2, // don't reweight for high energy
    };

  protected:

    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithTemplateXSec;
    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithoutTemplateXSec;

    std::map<int, double> x_FirstBinCenter, x_LastBinCenter;
    std::map<int, double> y_FirstBinCenter, y_LastBinCenter;
    std::map<int, double> z_FirstBinCenter, z_LastBinCenter;

    double ENuBoundary;

  public:

    CCQETemplateReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~CCQETemplateReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetTemplateReweight(double Enu_GeV, std::array<double, 2> bin_kin, double parameter_value);

    std::string GetCalculatorName() const { return "CCQETemplateReweightCalculator"; }

  };

  inline double CCQETemplateReweightCalculator::GetTemplateReweight(double Enu_GeV, std::array<double, 2> bin_kin, double parameter_value){

    int enu_range = (Enu_GeV<ENuBoundary) ? 0 : 1;

    if (enu_range == 0){

      static double epsil = 1E-6;
      double Enu_GeV_ForInterp = Enu_GeV;
      Enu_GeV_ForInterp = std::max( Enu_GeV_ForInterp, x_FirstBinCenter[enu_range] + epsil );
      Enu_GeV_ForInterp = std::min( Enu_GeV_ForInterp, x_LastBinCenter[enu_range] - epsil );

      double kin_Y_ForInterp = bin_kin[0];
      kin_Y_ForInterp = std::max( kin_Y_ForInterp, y_FirstBinCenter[enu_range] + epsil );
      kin_Y_ForInterp = std::min( kin_Y_ForInterp, y_LastBinCenter[enu_range] - epsil );

      double kin_Z_ForInterp = bin_kin[1];
      kin_Z_ForInterp = std::max( kin_Z_ForInterp, z_FirstBinCenter[enu_range] + epsil );
      kin_Z_ForInterp = std::min( kin_Z_ForInterp, z_LastBinCenter[enu_range] - epsil );

      /*
      printf("[CCQETemplateReweightCalculator::GetTemplateReweight] -> (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);
      printf("[CCQETemplateReweightCalculator::GetTemplateReweight] xsec (With Template, Without Template) = (%1.3f, %1.3e)\n", xsec_WithTemplate, xsec_WithoutTemplate);
      */
      double xsec_WithTemplate = map_ENuRange_to_WithTemplateXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp); // CV
      double xsec_WithoutTemplate = map_ENuRange_to_WithoutTemplateXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);

      double weight = ( xsec_WithoutTemplate * (1.-parameter_value) + xsec_WithTemplate * parameter_value ) / xsec_WithoutTemplate;
      // std::cout << "[CCQETemplateReweightCalculator] weight = " << weight << std::endl;

      // If weight is nan, find the closest grid point that isn't nan
      // TODO: is there some existing function I can use for this?
      if(weight!=weight){

        double best_weight = 1.0;
        bool found_valid = false;

        // Try the closest grid point first
        std::vector<double> x_centers, y_centers, z_centers;
        {
          int nx = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetXaxis()->GetNbins();
          int ny = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetYaxis()->GetNbins();
          int nz = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetZaxis()->GetNbins();
          for (int ix = 1; ix <= nx; ++ix)
            x_centers.push_back(map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetXaxis()->GetBinCenter(ix));
          for (int iy = 1; iy <= ny; ++iy)
            y_centers.push_back(map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetYaxis()->GetBinCenter(iy));
          for (int iz = 1; iz <= nz; ++iz)
            z_centers.push_back(map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetZaxis()->GetBinCenter(iz));
        }

        double test_x = find_closest_bin_center(Enu_GeV_ForInterp, x_centers);
        double test_y = find_closest_bin_center(kin_Y_ForInterp, y_centers);
        double test_z = find_closest_bin_center(kin_Z_ForInterp, z_centers);

        double x_for_interp;
        double y_for_interp;
        double z_for_interp;
        x_for_interp = std::max( test_x, x_FirstBinCenter[enu_range] + epsil );
        x_for_interp = std::min( x_for_interp, x_LastBinCenter[enu_range] - epsil );
        y_for_interp = std::max( test_y, y_FirstBinCenter[enu_range] + epsil );
        y_for_interp = std::min( y_for_interp, y_LastBinCenter[enu_range] - epsil );
        z_for_interp = std::max( test_z, z_FirstBinCenter[enu_range] + epsil );
        z_for_interp = std::min( z_for_interp, z_LastBinCenter[enu_range] - epsil );

        double xsec_WithTemplate_grid = map_ENuRange_to_WithTemplateXSec[enu_range]->Interpolate(x_for_interp, y_for_interp, z_for_interp);
        double xsec_WithoutTemplate_grid = map_ENuRange_to_WithoutTemplateXSec[enu_range]->Interpolate(x_for_interp, y_for_interp, z_for_interp);

        double grid_weight = ( xsec_WithoutTemplate_grid * (1.-parameter_value) + xsec_WithTemplate_grid * parameter_value ) / xsec_WithoutTemplate_grid;
        if(grid_weight == grid_weight) { 
          best_weight = grid_weight;
          found_valid = true;
        } 
        else {
          // If still nan, search all grid points for the closest non-nan
          double min_dist = std::numeric_limits<double>::max();
          int nx = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetXaxis()->GetNbins();
          int ny = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetYaxis()->GetNbins();
          int nz = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetZaxis()->GetNbins();
          for (int ix = 1; ix <= nx; ++ix) {
            double x = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetXaxis()->GetBinCenter(ix);
            for (int iy = 1; iy <= ny; ++iy) {
              double y = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetYaxis()->GetBinCenter(iy);
              for (int iz = 1; iz <= nz; ++iz) {
                double z = map_ENuRange_to_WithoutTemplateXSec[enu_range]->GetZaxis()->GetBinCenter(iz);
                double dist = std::sqrt(
                  std::pow(x-Enu_GeV_ForInterp,2) +
                  std::pow(y-kin_Y_ForInterp,2) +
                  std::pow(z-kin_Z_ForInterp,2)
                );
                
                x_for_interp = std::max( x, x_FirstBinCenter[enu_range] + epsil );
                x_for_interp = std::min( x_for_interp, x_LastBinCenter[enu_range] - epsil );
                y_for_interp = std::max( y, y_FirstBinCenter[enu_range] + epsil );
                y_for_interp = std::min( y_for_interp, y_LastBinCenter[enu_range] - epsil );
                z_for_interp = std::max( z, z_FirstBinCenter[enu_range] + epsil );
                z_for_interp = std::min( z_for_interp, z_LastBinCenter[enu_range] - epsil );

                double xsec_WithTemplate_try = map_ENuRange_to_WithTemplateXSec[enu_range]->Interpolate(x_for_interp, y_for_interp, z_for_interp);
                double xsec_WithoutTemplate_try = map_ENuRange_to_WithoutTemplateXSec[enu_range]->Interpolate(x_for_interp, y_for_interp, z_for_interp);
                double try_weight = ( xsec_WithoutTemplate_try * (1.-parameter_value) + xsec_WithTemplate_try * parameter_value ) / xsec_WithoutTemplate_try;
                if(try_weight == try_weight) { 
                  if(dist < min_dist) {
                    min_dist = dist;
                    best_weight = try_weight;
                    found_valid = true;
                  }
                }
              }
            }
          }
        }
        weight = best_weight;
        if(!found_valid) {
          // fallback
          weight = 1.;
        }
      }

      // clip
      weight = std::min(weight, 100.);

      return weight;
    }

  else if (enu_range == 1){
    return 1.;
    }

    // fallback
    return 1.;
  }

  inline void CCQETemplateReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = ps.get<std::string>("input_file", "");
    ENuBoundary = ps.get<double>("ENuBoundary");
    printf("[CCQETemplateReweightCalculator::GetTemplateReweight] ENuBoundary = %1.2f\n", ENuBoundary);

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file); // If specified per hist, replace it

      // if it does not start with "/", find it under ${NUSYSTEMATICS_FQ_DIR}/data/
      if(input_file.find("/")!=0){
        std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
        if(tmp_NUSYSTEMATICS_ROOT==""){
          throw invalid_CCQE_Template_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
        }
        input_file = tmp_NUSYSTEMATICS_ROOT+"/data/"+input_file;
      }

      if(hName=="LowE_WithTemplate"){
        std::cout << "--------LOADED HISTOGRAM: " << hName << std::endl;
        map_ENuRange_to_WithTemplateXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="LowE_WithoutTemplate"){
        std::cout << "--------LOADED HISTOGRAM: " << hName << std::endl;
        map_ENuRange_to_WithoutTemplateXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
    }

    for(int enu_range=0; enu_range<1; enu_range++){

      const auto& h_WithTemplateXsec = map_ENuRange_to_WithTemplateXSec[enu_range];

      const auto& XAxis_WithTemplateXsec = h_WithTemplateXsec->GetXaxis();
      const auto& YAxis_WithTemplateXsec = h_WithTemplateXsec->GetYaxis();
      const auto& ZAxis_WithTemplateXsec = h_WithTemplateXsec->GetZaxis();

      x_FirstBinCenter[enu_range] = XAxis_WithTemplateXsec->GetBinCenter(1);
      x_LastBinCenter[enu_range] = XAxis_WithTemplateXsec->GetBinCenter( XAxis_WithTemplateXsec->GetNbins() );

      y_FirstBinCenter[enu_range] = YAxis_WithTemplateXsec->GetBinCenter(1);
      y_LastBinCenter[enu_range] = YAxis_WithTemplateXsec->GetBinCenter( YAxis_WithTemplateXsec->GetNbins() );

      z_FirstBinCenter[enu_range] = ZAxis_WithTemplateXsec->GetBinCenter(1);
      z_LastBinCenter[enu_range] = ZAxis_WithTemplateXsec->GetBinCenter( ZAxis_WithTemplateXsec->GetNbins() );

      printf("@@ Enu range :%d\n", enu_range);
      printf("@@ - x-range: [%1.3f, %1.3f]\n", x_FirstBinCenter[enu_range], x_LastBinCenter[enu_range]);
      printf("@@ - y-range: [%1.3f, %1.3f]\n", y_FirstBinCenter[enu_range], y_LastBinCenter[enu_range]);
      printf("@@ - z-range: [%1.3f, %1.3f]\n", z_FirstBinCenter[enu_range], z_LastBinCenter[enu_range]);

    }

  }

} // namespace nusyst

#endif
