#pragma once

#include "systematicstools/interface/SystMetaData.hh"

#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace nusyst {

systtools::SystMetaData
ConfigureQEParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                            fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureMECParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                            fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureNCELParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                              fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureRESParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                             fhicl::ParameterSet &tool_options);
systtools::SystMetaData
ConfigureCOHParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                             fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureDISParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                             fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureFSIParameterHeaders(fhicl::ParameterSet const &, systtools::paramId_t,
                             fhicl::ParameterSet &tool_options);

systtools::SystMetaData
ConfigureOtherParameterHeaders(fhicl::ParameterSet const &,
                               systtools::paramId_t);

} // namespace nusyst
