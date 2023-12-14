####### Interface library
SET(IFCE_IMPLFILES)

SET(IFCE_HDRFILES
  IGENIESystProvider_tool.hh)


add_library(nusystematics_interface INTERFACE)
add_library(nusyst::interface ALIAS nusystematics_interface)

set_target_properties(nusystematics_interface PROPERTIES 
  PUBLIC_HEADER "${IFCE_HDRFILES}"
  EXPORT_NAME interface )

target_link_libraries(nusystematics_interface INTERFACE nusyst::commondeps GENIE3::All ROOT::Core systtools::interface)

install(TARGETS nusystematics_interface
    EXPORT nusyst-targets
    PUBLIC_HEADER DESTINATION include/nusystematics/interface COMPONENT Development)
