cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

set(CMAKE_CXX_STANDARD 17)

project(nusystematics VERSION 02.00.01 LANGUAGES CXX)

#Changes default install path to be a subdirectory of the build dir.
#Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${CMAKE_SYSTEM_NAME}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE DebWithRelInfo)
endif()

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(CPM)

find_package(ROOT 6.10 REQUIRED)

CPMFindPackage(
    NAME systematicstools
    GIT_TAG develop
    GITHUB_REPOSITORY NuSystematics/systematicstools
)

if(DEFINED ENV{NuHepMCModules_ROOT})
  include($ENV{NuHepMCModules_ROOT}/NuHepMCModules.cmake)
else()
  CPMAddPackage(
      NAME CMakeModules
      GIT_TAG stable
      GITHUB_REPOSITORY NuHepMC/CMakeModules
      DOWNLOAD_ONLY
  )
  include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)
endif()

CPMFindPackage(
    NAME Eigen3
    VERSION 3.4.0
    GITLAB_REPOSITORY "libeigen/eigen"
    GIT_TAG "3.4.0"
    OPTIONS "EIGEN_BUILD_DOC no" "BUILD_TESTING no"
)
## CMakeLists of Eigen3 include GNUInstallDirs, at it detects the system and change CMAKE_INSTALL_LIBDIR to lib64
set(CMAKE_INSTALL_LIBDIR "lib")

find_package(GENIE3 REQUIRED)
if(NOT TARGET GENIE3::All)
  message(FATAL_ERROR "Expected find_package(GENIE3 REQUIRED) call to set up target GENIE3::All.")
endif()

###### Compiler set up
add_library(nusystematics_dependencies INTERFACE)
target_link_libraries(nusystematics_dependencies INTERFACE 
  systtools::all)
target_include_directories(nusystematics_dependencies INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include> )
set_target_properties(nusystematics_dependencies PROPERTIES EXPORT_NAME commondeps)
add_library(nusyst::commondeps ALIAS nusystematics_dependencies)

install(TARGETS 
          nusystematics_dependencies
    EXPORT nusyst-targets)

add_subdirectory(src/nusystematics)

add_subdirectory(data)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/nusystematicsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(
  "${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/nusystematicsConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/nusystematicsConfig.cmake"
  INSTALL_DESTINATION cmake
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

install(EXPORT nusyst-targets
        NAMESPACE nusyst::
        DESTINATION lib/cmake/nusystematics)
install(FILES "${PROJECT_BINARY_DIR}/nusystematicsConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/nusystematicsConfig.cmake"
        DESTINATION lib/cmake/nusystematics)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/setup.nusystematics.sh.in 
  ${PROJECT_BINARY_DIR}/setup.nusystematics.sh @ONLY)
install(PROGRAMS ${PROJECT_BINARY_DIR}/setup.nusystematics.sh DESTINATION bin)

