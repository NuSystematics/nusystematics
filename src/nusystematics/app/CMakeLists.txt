LIST(APPEND TARGETS_TO_BUILD 
GenerateSystProviderConfigNuSyst
DumpConfiguredTweaksNuSyst
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ} 
    nusyst::systproviders
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
  )

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)
