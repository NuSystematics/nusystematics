#ifndef nusystematics_RESPONSE_CALCULATORS_EMISSENGINE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_EMISSENGINE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetCorrTailRW(double Emiss, 
                                     double parameter_value) {

    if (Emiss < 0.01 && Emiss > 0.){
      return parameter_value;
    }
    else {
      // TH: don't reweight events not in the correlated tail
      return 1;
    }
                                     }

} // namespace nusyst


#endif