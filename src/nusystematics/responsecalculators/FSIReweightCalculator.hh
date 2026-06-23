#ifndef nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"
#include "systematicstools/utility/string_parsers.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_FSI_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_FSI_FILEPATH);
using namespace std;
namespace nusyst {

  class FSIReweightCalculator{

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
    TH3D *hist_nom_2p;
    TH3D *hist_alt_2p;


  public:

    FSIReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~FSIReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetFSIReweight(double KEini, double Ebias, double parameter_value, int parpdg);
    double GetFSIReweight_2par(double KEini_0, double KEini_1, double Ebias, double parameter_value, int parpdg);

    std::string GetCalculatorName() const { return "FSIReweightCalculator"; }

  };

  inline double FSIReweightCalculator::GetFSIReweight(double KEini, double Ebias, double parameter_value, int parpdg){
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
    int idx_KEini = hist_nom->GetXaxis()->FindBin(KEini);
    int idx_Ebias = hist_nom->GetYaxis()->FindBin(Ebias);
    double weight_nom = hist_nom->GetBinContent(idx_KEini, idx_Ebias); // CV
    double weight_alt = hist_alt->GetBinContent(idx_KEini, idx_Ebias);
    //cout<<"idx_KEini "<<idx_KEini<<"; idx_Ebias "<<idx_Ebias<<endl;
    //cout<<"weight_nom "<<weight_nom<<endl;
    //cout<<"weight_alt "<<weight_alt<<endl;

    if(weight_nom==0.){
      return 1.;
    }

    double weight = ( weight_nom * (1.-parameter_value) + weight_alt * parameter_value ) / weight_nom;
    //cout<<"weight "<<weight<<endl;
     if(weight<0.001){
    //cout<<"weight_nom==0."<<endl;
    return 0.001;
  }
    return weight;

  }

inline double FSIReweightCalculator::GetFSIReweight_2par(double KEini_0, double KEini_1, double Ebias, double parameter_value, int parpdg){
  TH3D *hist_nom, *hist_alt;
  
  hist_nom = hist_nom_2p;
  hist_alt = hist_alt_2p;

  int idx_KEini_0 = hist_nom->GetXaxis()->FindBin(KEini_0);
  int idx_KEini_1 = hist_nom->GetYaxis()->FindBin(KEini_1);
  int idx_Ebias = hist_nom->GetZaxis()->FindBin(Ebias);
  double weight_nom = hist_nom->GetBinContent(idx_KEini_0, idx_KEini_1, idx_Ebias); // CV
  double weight_alt = hist_alt->GetBinContent(idx_KEini_0, idx_KEini_1, idx_Ebias);
  //cout<<"idx_KEini "<<idx_KEini<<"; idx_Ebias "<<idx_Ebias<<endl;
  //cout<<"weight_nom "<<weight_nom<<endl;
  //cout<<"weight_alt "<<weight_alt<<endl;

  if(weight_nom==0.){
    //cout<<"weight_nom==0."<<endl;
    return 1.;
  }

  double weight = ( weight_nom * (1.-parameter_value) + weight_alt * parameter_value ) / weight_nom;
  //cout<<"weight "<<weight<<endl;

  return weight;

}

  inline void FSIReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = systtools::expand_env_vars( ps.get<std::string>("input_file", "") );

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = systtools::expand_env_vars( val_config.get<std::string>("input_file", default_root_file) ); // If specified per hist, replace it

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

    }
  }
} // namespace nusyst

#endif
