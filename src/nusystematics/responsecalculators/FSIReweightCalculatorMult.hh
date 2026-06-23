#ifndef nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculatorMult_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculatorMult_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_FSIMult_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_FSIMult_FILEPATH);
using namespace std;
namespace nusyst {

  class FSIReweightCalculatorMult{

    enum ENuRange {
      LowE = 0,
      HighE = 20,
    };

  protected:

    TH2D *hist_nom_protonPlus;
    TH2D *hist_alt_protonPlus;
    TH2D *hist_nom_neutron;
    TH2D *hist_alt_neutron;
    TH2D *hist_nom_piPlus;
    TH2D *hist_alt_piPlus;
    TH2D *hist_nom_pi0;
    TH2D *hist_alt_pi0;
    TH2D *hist_nom_piMinus;
    TH2D *hist_alt_piMinus;


    TH2D *hist_nomDiff_protonPlus;
    TH2D *hist_altDiff_protonPlus;
    TH2D *hist_nomDiff_neutron;
    TH2D *hist_altDiff_neutron;
    TH2D *hist_nomDiff_piPlus;
    TH2D *hist_altDiff_piPlus;
    TH2D *hist_nomDiff_pi0;
    TH2D *hist_altDiff_pi0;
    TH2D *hist_nomDiff_piMinus;
    TH2D *hist_altDiff_piMinus;


  public:

    FSIReweightCalculatorMult(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~FSIReweightCalculatorMult(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetFSIReweightMultSum(double KEini, double nucleons, double parameter_value, int parpdg);
    double GetFSIReweightMultDiff(double KEini, double nucleons, double parameter_value, int parpdg);

    std::string GetCalculatorName() const { return "FSIReweightCalculatorMult"; }

  };

  inline double FSIReweightCalculatorMult::GetFSIReweightMultSum(double KEini, double nucleons, double parameter_value, int parpdg){
    TH2D *hist_nom, *hist_alt;
    if (parpdg == 2212) {
      hist_nom = hist_nom_protonPlus;
      hist_alt = hist_alt_protonPlus;
    }
    else if (parpdg == 2112) {
      hist_nom = hist_nom_neutron;
      hist_alt = hist_alt_neutron;
    }
    else if (parpdg == 211) {
      hist_nom = hist_nom_piPlus;
      hist_alt = hist_alt_piPlus;
    }
    else if (parpdg == 111) {
      hist_nom = hist_nom_pi0;
      hist_alt = hist_alt_pi0;
    }
    else if (parpdg == -211) {
      hist_nom = hist_nom_piMinus;
      hist_alt = hist_alt_piMinus;
    }
    else {
      return 1.;
    }
    //std::cout<<nucleons<<std::endl;
    int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_nucleons = hist_nom->GetYaxis()->FindBin(nucleons);
    double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_nucleons); // CV
    double weight_alt = hist_alt->GetBinContent(idx_KEini, idx_nucleons);
    //cout<<"idx_KEini "<<idx_KEini<<"; idx_nucleons "<<idx_nucleons<<endl;
    //cout<<"weight_nom "<<weight_nom<<endl;
    //cout<<"weight_alt "<<weight_alt<<endl;

    if(weight_nom==0.){
      return 1.;
    }

    double weight = ( weight_nom * (1.-parameter_value) + weight_alt * parameter_value ) / weight_nom;
    //cout<<"weight "<<weight<<endl;
     if(weight<0.01){
    //cout<<"weight_nom==0."<<endl;
    return 0.01;
  }
    return weight;

  }

inline double FSIReweightCalculatorMult::GetFSIReweightMultDiff(double KEini, double nucleons, double parameter_value, int parpdg){
    TH2D *hist_nom, *hist_alt;
    if (parpdg == 2212) {
      hist_nom = hist_nomDiff_protonPlus;
      hist_alt = hist_altDiff_protonPlus;
    }
    else if (parpdg == 2112) {
      hist_nom = hist_nomDiff_neutron;
      hist_alt = hist_altDiff_neutron;
    }
    else if (parpdg == 211) {
      hist_nom = hist_nomDiff_piPlus;
      hist_alt = hist_altDiff_piPlus;
    }
    else if (parpdg == 111) {
      hist_nom = hist_nomDiff_pi0;
      hist_alt = hist_altDiff_pi0;
    }
    else if (parpdg == -211) {
      hist_nom = hist_nomDiff_piMinus;
      hist_alt = hist_altDiff_piMinus;
    }
    else {
      return 1.;
    }

  int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_nucleons = hist_nom->GetYaxis()->FindBin(nucleons);
  double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_nucleons); // CV
  double weight_alt = hist_alt->GetBinContent(idx_KEini,  idx_nucleons);
  //cout<<"idx_KEini "<<idx_KEini<<"; idx_nucleons "<<idx_nucleons<<endl;
  //cout<<"weight_nom "<<weight_nom<<endl;
  //cout<<"weight_alt "<<weight_alt<<endl;

  if(weight_nom==0.){
    //cout<<"weight_nom==0."<<endl;
    return 1.;
  }

  double weight = ( weight_nom * (1.-parameter_value) + weight_alt * parameter_value ) / weight_nom;
  //cout<<"weight "<<weight<<endl;
  if(weight<0.01){
    //cout<<"weight_nom==0."<<endl;
    return 0.01;
  }

  return weight;

}

  inline void FSIReweightCalculatorMult::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = ps.get<std::string>("input_file", "");

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file); // If specified per hist, replace it

      // if it does not start with "/", find it under ${NUSYSTEMATICS_FQ_DIR}/data/
      if(input_file.find("/")!=0){
        std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
        if(tmp_NUSYSTEMATICS_ROOT==""){
          throw invalid_FSIMult_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
        }
        input_file = tmp_NUSYSTEMATICS_ROOT+"/data/"+input_file;
      }

      if(hName=="hist_nom_protonPlus"){
        hist_nom_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_alt_protonPlus"){
        hist_alt_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nom_neutron"){
        hist_nom_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_alt_neutron"){
        hist_alt_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nom_piPlus"){
        hist_nom_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_alt_piPlus"){
        hist_alt_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nom_pi0"){
        hist_nom_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_alt_pi0"){
        hist_alt_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nom_piMinus"){
        hist_nom_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_alt_piMinus"){
        hist_alt_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }

      if(hName=="hist_nomDiff_protonPlus"){
        hist_nomDiff_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_altDiff_protonPlus"){
        hist_altDiff_protonPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nomDiff_neutron"){
        hist_nomDiff_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_altDiff_neutron"){
        hist_altDiff_neutron = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nomDiff_piPlus"){
        hist_nomDiff_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_altDiff_piPlus"){
        hist_altDiff_piPlus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nomDiff_pi0"){
        hist_nomDiff_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_altDiff_pi0"){
        hist_altDiff_pi0 = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_nomDiff_piMinus"){
        hist_nomDiff_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_altDiff_piMinus"){
        hist_altDiff_piMinus = GetHistogram<TH2D>(input_file, input_hist);
      }


    }
  }
} // namespace nusyst

#endif
