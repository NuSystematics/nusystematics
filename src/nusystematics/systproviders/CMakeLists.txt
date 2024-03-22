####### systproviders library
SET(SYSTP_IMPLFILES
  BeRPAWeight_tool.cc
  EbLepMomShift_tool.cc
  FSILikeEAvailSmearing_tool.cc
  MINERvAE2p2h_tool.cc
  MINERvAq0q3Weighting_tool.cc
  MKSinglePiTemplate_tool.cc
  MiscInteractionSysts_tool.cc
  NOvAStyleNonResPionNorm_tool.cc
  DIRT2_Emiss_tool.cc
  )

if(GENIEReWeight_ENABLED)
  LIST(APPEND SYSTP_IMPLFILES
    GENIEReWeightEngineConfig.cc
    GENIEReWeightParamConfig.cc
    GENIEReWeight_tool.cc
    SkeleWeighter_tool.cc
    ZExpPCAWeighter_tool.cc
    ResIso_tool.cc)
endif()

SET(SYSTP_HDRFILES
  BeRPAWeight_tool.hh
  EbLepMomShift_tool.hh
  FSILikeEAvailSmearing_tool.hh
  MINERvAE2p2h_tool.hh
  MINERvAq0q3Weighting_tool.hh
  MKSinglePiTemplate_tool.hh
  MiscInteractionSysts_tool.hh
  NOvAStyleNonResPionNorm_tool.hh
  DIRT2_Emiss_tool.hh
)

if(GENIEReWeight_ENABLED)
  LIST(APPEND SYSTP_HDRFILES
    GENIEReWeightEngineConfig.hh
    GENIEReWeightParamConfig.hh
    GENIEReWeight_tool.hh
    GENIEResponseParameterAssociation.hh
    SkeleWeighter_tool.hh
    ZExpPCAWeighter_tool.hh
    ResIso_tool.hh)
endif()

add_library(nusystematics_systproviders SHARED ${SYSTP_IMPLFILES})
add_library(nusyst::systproviders ALIAS nusystematics_systproviders)

set_target_properties(nusystematics_systproviders PROPERTIES 
  PUBLIC_HEADER "${SYSTP_HDRFILES}"
  EXPORT_NAME systproviders )

target_link_libraries(nusystematics_systproviders PUBLIC nusyst::responsecalculators systtools::all Eigen3::Eigen)

install(TARGETS nusystematics_systproviders
    EXPORT nusyst-targets
    PUBLIC_HEADER DESTINATION include/nusystematics/systproviders COMPONENT Development)
