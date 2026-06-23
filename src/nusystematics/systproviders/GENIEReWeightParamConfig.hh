#pragma once

#include "systematicstools/interface/SystMetaData.hh"

#include "yaml-cpp/yaml.h"

#include <string>

namespace nusyst {

systtools::SystMetaData
ConfigureQEParameterHeaders(YAML::Node const &, systtools::paramId_t,
                            YAML::Node &tool_options);

systtools::SystMetaData
ConfigureMECParameterHeaders(YAML::Node const &, systtools::paramId_t,
                            YAML::Node &tool_options);

systtools::SystMetaData
ConfigureNCELParameterHeaders(YAML::Node const &, systtools::paramId_t,
                              YAML::Node &tool_options);

systtools::SystMetaData
ConfigureRESParameterHeaders(YAML::Node const &, systtools::paramId_t,
                             YAML::Node &tool_options);
systtools::SystMetaData
ConfigureCOHParameterHeaders(YAML::Node const &, systtools::paramId_t,
                             YAML::Node &tool_options);

systtools::SystMetaData
ConfigureDISParameterHeaders(YAML::Node const &, systtools::paramId_t,
                             YAML::Node &tool_options);

systtools::SystMetaData
ConfigureFSIParameterHeaders(YAML::Node const &, systtools::paramId_t,
                             YAML::Node &tool_options);

systtools::SystMetaData
ConfigureOtherParameterHeaders(YAML::Node const &,
                               systtools::paramId_t);

} // namespace nusyst
