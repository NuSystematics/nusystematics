LIST(APPEND TARGETS_TO_BUILD 
GenerateSystProviderConfigNuSyst
DumpConfiguredTweaksNuSyst
PlotSystVariationsNuSyst
DeclaredDialTestNuSyst
DumpDialResponseNuSyst
GenerateAllDialsConfigNuSyst
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ}
    nusyst::systproviders
    nusyst::utility
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
  )
  # Bake the install prefix in at compile time so InventoryDefaultPath()
  # can resolve the kitchen-sink cache against the binary's OWN install,
  # rather than against env vars (nusystematics_ROOT / NUSYST) that can be
  # overwritten by other software's setup scripts.
  target_compile_definitions(${targ} PRIVATE
    NUSYST_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)

# Unified-CLI shell dispatcher: a single `nusyst` entrypoint that forwards
# subcommands (config / inventory / tweaks / response / plots) to the binaries
# above. Installed as a program so the executable bit is preserved.
install(PROGRAMS nusyst DESTINATION bin/)
