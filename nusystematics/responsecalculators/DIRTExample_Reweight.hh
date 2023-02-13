#ifndef nusystematics_RESPONSE_CALCULATORS_DIRTEXAMPLE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_DIRTEXAMPLE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetXSecCV(double q3, double q0){
    if(q0<0.) return 1.; // useless line to avoid unused-parameter error

    if(q3<0.) return 1.; // useless line to avoid unused-parameter error
    else return 1.;
  }
  inline double GetXSecAltModelA(double q3, double q0){
    if(q0<0) return 0.; // useless line to avoid unused-parameter error

    if(q3>0.5) return 0.;
    else return 1.1;
  }
  inline double GetXSecAltModelB(double q3, double q0){
    if(q0<0) return 0.; // useless line to avoid unused-parameter error

    if(q3>0.8) return 0.;
    else return 1.3;
  }

  inline double GetReweightToAltModelA(double q3, double q0, 
                                     double parameter_value) {

    // parameter_value = 0 : CV
    // parameter_value = 1 : AltModelA
    double xsec_cv = GetXSecCV(q3, q0);
    double xsec_alt = GetXSecAltModelA(q3, q0);
    double new_xsec = parameter_value*xsec_alt + (1.-parameter_value)*xsec_cv;
    return new_xsec/xsec_cv;

  }

  inline double GetReweightToAltModelB(double q3, double q0,
                                     double parameter_value) {

    // parameter_value = 0 : CV
    // parameter_value = 1 : AltModelB
    double xsec_cv = GetXSecCV(q3, q0);
    double xsec_alt = GetXSecAltModelB(q3, q0);
    double new_xsec = parameter_value*xsec_alt + (1.-parameter_value)*xsec_cv;
    return new_xsec/xsec_cv;

  }

} // namespace nusyst

#endif
