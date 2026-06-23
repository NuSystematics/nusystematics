####### utility library
SET(UTLY_IMPLFILES)

SET(UTLY_HDRFILES
  GENIEUtils.hh
  channel_classification.hh
  enumclass2int.hh
  exceptions.hh
  genie_metadata.hh
  simbUtility.hh
  make_instance.hh
  response_helper.hh
  silence_genie.hh
  KinVarUtils.hh
)


add_library(nusystematics_utility INTERFACE)
add_library(nusyst::utility ALIAS nusystematics_utility)

set_target_properties(nusystematics_utility PROPERTIES 
  PUBLIC_HEADER "${UTLY_HDRFILES}"
  EXPORT_NAME utility )

target_link_libraries(nusystematics_utility INTERFACE nusyst::interface)

install(TARGETS nusystematics_utility
    EXPORT nusyst-targets
    PUBLIC_HEADER DESTINATION include/nusystematics/utility COMPONENT Development)
