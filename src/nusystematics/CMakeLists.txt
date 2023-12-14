add_subdirectory(interface)
add_subdirectory(responsecalculators)
add_subdirectory(systproviders)
add_subdirectory(utility)
add_subdirectory(app)

add_library(nusystematics_all INTERFACE)
target_link_libraries(nusystematics_all INTERFACE nusyst::systproviders)

set_target_properties(nusystematics_all PROPERTIES EXPORT_NAME all)

add_library(nusyst::all ALIAS nusystematics_all)

install(TARGETS 
          nusystematics_all
    EXPORT nusyst-targets)
