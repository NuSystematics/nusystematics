#ifndef nusystematics_RESPONSE_CALCULATORS_EMISSENGINE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_EMISSENGINE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetEmissCorrTailRW(double Emiss_preFSI, double parameter_value) {

    if (Emiss_preFSI < 0.01 && Emiss_preFSI > 0.){
      // TH: flat weight applied to correlated tail to shift it up/down
	  return parameter_value;
    }
    else {
      // TH: don't reweight events not in the correlated tail
      return 1;
    }
  }

  inline double GetEmissLinearRW(double Emiss_preFSI, double parameter_value){

    if (Emiss_preFSI > 0.01){
      // TH: linear re-weight to increase mean of removal energy distribution
	  double weight = 1 + parameter_value * Emiss_preFSI;
      if (weight > 0){
        return weight;
      }
      else {
		// TH: no negative weights!
        return 0;
      }
    }
    else{
      // TH: don't reweight events in the correlated tail
      return 1;
    }
  }

  inline double GetEmissTrigRW(double Emiss_preFSI, double parameter_value){
    // TH: dial to move the peak of the removal energy distribution
	//     Do this by weighting down the peak at 10 MeV and weighting up a higher energy region using a cosine weight distribution 
    if (Emiss_preFSI > 0.01 && Emiss_preFSI < 0.025){
	  double weight = 1 + parameter_value * TMath::Cos( 0.1 * TMath::Pi() * (Emiss_preFSI/0.001));
      if (weight > 0.2){
        return weight;
      }
      else{
		// TH: cap weights at 0.2 to avoid discontinuity/zero weights at ~10 MeV
        return 0.2;
      }
    }
    else {
      // TH: don't re-weight events in the correlated tail or at high Emiss
      return 1;
    }
    
  }
                                  

} // namespace nusyst


#endif
