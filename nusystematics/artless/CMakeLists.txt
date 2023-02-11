###### SystematicsTools set up
if (NOT DEFINED SYSTTOOLS_DIR)
  include(ExternalProject)

  ExternalProject_Add(systematicstools
    PREFIX "${PROJECT_BINARY_DIR}/systematicstools"
    GIT_REPOSITORY https://github.com/jedori0228/systematicstools.git
    GIT_TAG feature/jskim_artless
    UPDATE_DISCONNECTED 1
    CMAKE_ARGS -DUSEART=0
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})

  set(SYSTTOOLS_PATH ${CMAKE_INSTALL_PREFIX})
  set(EXTERNAL_SYSTTOOLS 1)
else()
  set(SYSTTOOLS_PATH ${SYSTTOOLS_DIR})
  set(EXTERNAL_SYSTTOOLS 0)
endif()

if(CMAKE_LINK_FLAGS STREQUAL "")
  SET(CMAKE_LINK_FLAGS "-L${SYSTTOOLS_PATH}/lib")
else()
  SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} -L${SYSTTOOLS_PATH}/lib")
endif()

set(SYSTTOOLS_LIBS
  systematicstools_interface
  systematicstools_interpreters
  systematicstools_systproviders
  systematicstools_utility)

###### ROOT set up

if (NOT DEFINED ENV{ROOTSYS})
  message (FATAL_ERROR "$ROOTSYS is not defined, please set up ROOT first.")
else()
  message(STATUS "Using ROOT installed at $ENV{ROOTSYS}")
  set(CMAKE_ROOTSYS $ENV{ROOTSYS})
endif()

#Get cflags from ROOT
execute_process (COMMAND root-config
  --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
#Get libdir from ROOT
execute_process (COMMAND root-config
  --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
#Get version from ROOT
execute_process (COMMAND root-config
  --version OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
#Get features from ROOT
execute_process (COMMAND root-config
  --features OUTPUT_VARIABLE ROOT_FEATURES OUTPUT_STRIP_TRAILING_WHITESPACE)

SET(ROOT_LIBS)
LIST(APPEND ROOT_LIBS
  Core
  RIO
  XMLIO
  Net
  Hist
  Graf
  Graf3d
  Gpad
  Tree
  Rint
  Postscript
  Matrix
  Physics
  MathCore
  MathMore
  Thread
  EG
  Geom
  GenVector
  EGPythia6
  Pythia6)

message (STATUS "[ROOT]: root-config --version: ${ROOT_VERSION} ")
message (STATUS "[ROOT]: root-config --cflags : ${ROOT_CXX_FLAGS} ")
message (STATUS "[ROOT]: libs in use          : ${ROOT_LIBS} ")

###### GENIE setup
  include(${CMAKE_CURRENT_SOURCE_DIR}/GENIESetup.cmake)
###### Compiler set up
#-Wl,--no-as-needed used to ensuring dynamic linking of GENIE algorithm libraries so that TClass::GetClass has access to the dictionaries at runtime
set(EXTRA_LINK_FLAGS "-Wl,--no-as-needed -L${CMAKE_INSTALL_PREFIX}/lib ${STR_GENIE_LINK_DIRS} -L${ROOT_LIBDIR}")
if(CMAKE_LINK_FLAGS STREQUAL "")
  SET(CMAKE_LINK_FLAGS "${EXTRA_LINK_FLAGS}")
else()
  SET(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${EXTRA_LINK_FLAGS}")
endif()

set(EXTRA_CXX_FLAGS "-DNO_ART -Wall -Wextra -Werror ${ROOT_CXX_FLAGS} ${GENIE_CXX_FLAGS} -Wno-misleading-indentation")
if(CMAKE_CXX_FLAGS STREQUAL "")
  SET(CMAKE_CXX_FLAGS ${EXTRA_CXX_FLAGS})
else()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()


if(NOT EXTERNAL_SYSTTOOLS)
  include_directories(${SYSTTOOLS_DIR}/include)
endif()
include_directories(${GENIE_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_INSTALL_PREFIX}/include)

####### utility library
SET(UTIL_HDRFILES
  ${CMAKE_SOURCE_DIR}/nusystematics/utility/enumclass2int.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/utility/exceptions.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/utility/GENIEUtils.hh)

INSTALL(FILES ${UTIL_HDRFILES} DESTINATION include/nusystematics/utility)

####### responsecalculators library
SET(RC_HDRFILES
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/BeRPA.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/MKSinglePiTemplate_ReWeight.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/MINERvA2p2hq0q3.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/MINERvAq0q3Weighting_data.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/MINERvARPAq0q3_ReWeight.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/responsecalculators/DIRTExample_Reweight.hh)

INSTALL(FILES ${RC_HDRFILES} DESTINATION include/nusystematics/responsecalculators)

####### systproviders library
SET(SP_IMPLFILES
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/GENIEReWeight_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/GENIEReWeightEngineConfig.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/GENIEReWeightParamConfig.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MKSinglePiTemplate_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MINERvAq0q3Weighting_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/NOvAStyleNonResPionNorm_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/BeRPAWeight_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MiscInteractionSysts_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MINERvAE2p2h_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/EbLepMomShift_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/FSILikeEAvailSmearing_tool.cc
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/DIRTExample_tool.cc
)

SET(SP_HDRFILES
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/GENIEReWeight_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MKSinglePiTemplate_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/GENIEResponseParameterAssociation.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/NOvAStyleNonResPionNorm_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/BeRPAWeight_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MiscInteractionSysts_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/MINERvAE2p2h_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/EbLepMomShift_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/FSILikeEAvailSmearing_tool.hh
  ${CMAKE_SOURCE_DIR}/nusystematics/systproviders/DIRTExample_tool.hh
)

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
  add_library(nusystematics_systproviders STATIC ${SP_IMPLFILES})
else(CMAKE_BUILD_TYPE MATCHES RELEASE)
  add_library(nusystematics_systproviders SHARED ${SP_IMPLFILES})
endif()

if(EXTERNAL_SYSTTOOLS)
  add_dependencies(nusystematics_systproviders systematicstools)
endif()

INSTALL(TARGETS nusystematics_systproviders DESTINATION lib)
INSTALL(FILES ${SP_HDRFILES} DESTINATION include/nusystematics/systproviders)

####### GENIE RwExample
add_executable(ReWeightGHep ReWeightGHep.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(ReWeightGHep systematicstools)
endif()

target_link_libraries(ReWeightGHep nusystematics_systproviders)

set_target_properties(ReWeightGHep PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})
target_link_libraries(ReWeightGHep ${SYSTTOOLS_LIBS})
target_link_libraries(ReWeightGHep ${GENIE_LIBS})
target_link_libraries(ReWeightGHep ${ROOT_LIBS})

INSTALL(TARGETS ReWeightGHep DESTINATION bin)

####### GenerateSystProvider app
add_executable(GenerateSystProviderConfigNuSyst ${CMAKE_SOURCE_DIR}/nusystematics/artless/GenerateSystProviderConfigNuSyst.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(GenerateSystProviderConfigNuSyst systematicstools)
endif()
set_target_properties(GenerateSystProviderConfigNuSyst PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(GenerateSystProviderConfigNuSyst nusystematics_systproviders)
target_link_libraries(GenerateSystProviderConfigNuSyst ${SYSTTOOLS_LIBS})
target_link_libraries(GenerateSystProviderConfigNuSyst ${GENIE_LIBS})
target_link_libraries(GenerateSystProviderConfigNuSyst ${ROOT_LIBS})

INSTALL(TARGETS GenerateSystProviderConfigNuSyst DESTINATION bin)

####### DumpPrecalculatedPolyResponse app
add_executable(DumpPrecalculatedPolyResponse ${CMAKE_SOURCE_DIR}/nusystematics/artless/DumpPrecalculatedPolyResponse.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(DumpPrecalculatedPolyResponse systematicstools)
endif()
set_target_properties(DumpPrecalculatedPolyResponse PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(DumpPrecalculatedPolyResponse nusystematics_systproviders)
target_link_libraries(DumpPrecalculatedPolyResponse ${SYSTTOOLS_LIBS})
target_link_libraries(DumpPrecalculatedPolyResponse ${GENIE_LIBS})
target_link_libraries(DumpPrecalculatedPolyResponse ${ROOT_LIBS})

INSTALL(TARGETS DumpPrecalculatedPolyResponse DESTINATION bin)

####### TestPolyResponse5 app
add_executable(TestPolyResponse5 ${CMAKE_SOURCE_DIR}/nusystematics/artless/TestPolyResponse5.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(TestPolyResponse5 systematicstools)
endif()
set_target_properties(TestPolyResponse5 PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(TestPolyResponse5 nusystematics_systproviders)
target_link_libraries(TestPolyResponse5 ${SYSTTOOLS_LIBS})
target_link_libraries(TestPolyResponse5 ${GENIE_LIBS})
target_link_libraries(TestPolyResponse5 ${ROOT_LIBS})

INSTALL(TARGETS TestPolyResponse5 DESTINATION bin)

####### FindISystProviderNuSyst app
add_executable(FindISystProviderNuSyst ${CMAKE_SOURCE_DIR}/nusystematics/artless/FindISystProviderNuSyst.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(FindISystProviderNuSyst systematicstools)
endif()
set_target_properties(FindISystProviderNuSyst PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(FindISystProviderNuSyst nusystematics_systproviders)
target_link_libraries(FindISystProviderNuSyst ${SYSTTOOLS_LIBS})
target_link_libraries(FindISystProviderNuSyst ${GENIE_LIBS})
target_link_libraries(FindISystProviderNuSyst ${ROOT_LIBS})

INSTALL(TARGETS FindISystProviderNuSyst DESTINATION bin)

####### DumpConfiguredTweaksNuSyst app
add_executable(DumpConfiguredTweaksNuSyst ${CMAKE_SOURCE_DIR}/nusystematics/app/DumpConfiguredTweaksNuSyst.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(DumpConfiguredTweaksNuSyst systematicstools)
endif()
set_target_properties(DumpConfiguredTweaksNuSyst PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(DumpConfiguredTweaksNuSyst nusystematics_systproviders)
target_link_libraries(DumpConfiguredTweaksNuSyst ${SYSTTOOLS_LIBS})
target_link_libraries(DumpConfiguredTweaksNuSyst ${GENIE_LIBS})
target_link_libraries(DumpConfiguredTweaksNuSyst ${ROOT_LIBS})

INSTALL(TARGETS DumpConfiguredTweaksNuSyst DESTINATION bin)

####### DumpConfiguredTweaksNuSyst app
add_executable(BindingEnergyFlatTreeMaker ${CMAKE_SOURCE_DIR}/nusystematics/inputgenerationtools/BindingEnergyFlatTreeMaker.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(BindingEnergyFlatTreeMaker systematicstools)
endif()
set_target_properties(BindingEnergyFlatTreeMaker PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(BindingEnergyFlatTreeMaker nusystematics_systproviders)
target_link_libraries(BindingEnergyFlatTreeMaker ${SYSTTOOLS_LIBS})
target_link_libraries(BindingEnergyFlatTreeMaker ${GENIE_LIBS})
target_link_libraries(BindingEnergyFlatTreeMaker ${ROOT_LIBS})

INSTALL(TARGETS BindingEnergyFlatTreeMaker DESTINATION bin)


####### Validation app
add_executable(TweaksNuSyst_Validate ${CMAKE_SOURCE_DIR}/nusystematics/scripts/validation/TweaksNuSyst_Validate.cc)
if(EXTERNAL_SYSTTOOLS)
  add_dependencies(TweaksNuSyst_Validate systematicstools)
endif()
set_target_properties(TweaksNuSyst_Validate PROPERTIES LINK_FLAGS ${CMAKE_LINK_FLAGS})

target_link_libraries(TweaksNuSyst_Validate ${ROOT_LIBS})

INSTALL(TARGETS TweaksNuSyst_Validate DESTINATION bin)

####### interface
INSTALL(FILES ${CMAKE_SOURCE_DIR}/nusystematics/interface/IGENIESystProvider_tool.hh DESTINATION include/nusystematics/interface)

####### response_helper
INSTALL(FILES ${CMAKE_SOURCE_DIR}/nusystematics/artless/response_helper.hh DESTINATION include/nusystematics/artless)

####### fhicl files
file(GLOB FCL ${CMAKE_SOURCE_DIR}/nusystematics/fcl/*.fcl)
install(FILES ${FCL} DESTINATION ${CMAKE_INSTALL_PREFIX}/fcl)

####### environment scripts
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.sh" DESTINATION
  ${CMAKE_INSTALL_PREFIX})
