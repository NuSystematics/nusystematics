####### responsecalculators library
SET(RSPC_IMPLFILES)

SET(RSPC_HDRFILES
  BeRPA.hh
  MINERvA2p2hEnergyDependenceScale.hh
  MKSinglePiTemplate_ReWeight.hh
  nuenumu_xsec_ratio.hh
  C12ToAr40_2p2hScaling.hh
  MINERvA2p2hq0q3.hh
  SPPLowQ2Suppression.hh
  EnuBinnedTemplateResponseCalculator.hh
  MINERvARPAq0q3_ReWeight.hh
  TemplateResponseCalculatorBase.hh
  FSILikeEAvailSmearing.hh
  MINERvAq0q3Weighting_data.hh
  nuenuebar_xsec_ratio.hh
  DIRT2_EmissEngine_Reweight.hh
)


add_library(nusystematics_responsecalculators INTERFACE)
add_library(nusyst::responsecalculators ALIAS nusystematics_responsecalculators)

set_target_properties(nusystematics_responsecalculators PROPERTIES 
  PUBLIC_HEADER "${RSPC_HDRFILES}"
  EXPORT_NAME responsecalculators )

target_link_libraries(nusystematics_responsecalculators INTERFACE nusyst::utility)

install(TARGETS nusystematics_responsecalculators
    EXPORT nusyst-targets
    PUBLIC_HEADER DESTINATION include/nusystematics/responsecalculators COMPONENT Development)
