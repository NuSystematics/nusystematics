#ifndef nusystematics_RESPONSE_CALCULATORS_RPAENGINE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_RPAENGINE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetRPA_RW(double q0, double parameter_value) {
    
    // TH: blank for now, just a template
    double weight = 1 + q0*parameter_value;
    return weight;    
  }
                                  

} // namespace nusyst


#endif
