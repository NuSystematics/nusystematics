get_filename_component(nusystematicsDependencies_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

include(${nusystematicsDependencies_CMAKE_DIR}/CPM.cmake)
CPMFindPackage(
    NAME CMakeModules
    GIT_TAG stable
    GITHUB_REPOSITORY NuHepMC/CMakeModules
    DOWNLOAD_ONLY
)
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

include(message)

find_package(systematicstools 23.06 REQUIRED)
find_package(GENIE3 REQUIRED)

if(DEFINED ROOT_CXX_STANDARD)
  if(NOT DEFINED CMAKE_CXX_STANDARD OR ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
  endif()
endif()