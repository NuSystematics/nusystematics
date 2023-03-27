#ifndef nusystematics_RESPONSE_CALCULATORS_RPAENGINE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_RPAENGINE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetRPA_LowETransfer_0_RW(double q0, double parameter_value) {
    if (q0 > 0.01 && q0 <= 0.06){
      return parameter_value;
    }
    else {
      return 1;
    }  
  }

  inline double GetRPA_LowETransfer_1_RW(double q0, double parameter_value) {
    if (q0 > 0.06 && q0 <= 0.11){
      return parameter_value;
    }
    else {
      return 1;
    }  
  }

  inline double GetRPA_LowETransfer_2_RW(double q0, double parameter_value) {
    if (q0 > 0.11 && q0 <= 0.16){
      return parameter_value;
    }
    else {
      return 1;
    }  
  }

  inline double GetRPA_LowETransfer_3_RW(double q0, double parameter_value) {
    if (q0 > 0.16 && q0 <= 0.21){
      return parameter_value;
    }
    else {
      return 1;
    }  
  }

  inline double GetRPA_HighETransfer_0_RW(double q0, double parameter_value) {
    if (q0 > 0.21){
      return parameter_value;
    }
    else {
      return 1;
    }  
  }                  

} // namespace nusyst


#endif
