# Here we define various varaibles that can be used inside c++ modules

# 1) GENIE version in more detail

string(REPLACE "." ";" VERSION_LIST ${GENIE_VERSION})
list(GET VERSION_LIST 0 MAJOR)
list(GET VERSION_LIST 1 MINOR)
list(GET VERSION_LIST 2 PATCH)
math(EXPR GENIE_VERSION_CODE "${MAJOR} * 10000 + ${MINOR} * 100 + ${PATCH}")
message(STATUS "GENIE_VERSION_CODE: ${GENIE_VERSION_CODE}")

# AR25 FSI dials
# TODO These modules are not yet merged into tagged Reweight
# https://github.com/GENIE-MC/Reweight/pull/44
# https://github.com/GENIE-MC/Reweight/pull/45
# Assuming these are activated for GENIE>=3.08.00 for now (01/12/2026),
# but need to be updated
if(GENIE_VERSION_CODE GREATER_EQUAL 30800)
  set(BUILD_AR25_FSI_DIALS TRUE)
else()
  set(BUILD_AR25_FSI_DIALS FALSE)
endif()

message(STATUS "BUILD_AR25_FSI_DIALS: ${BUILD_AR25_FSI_DIALS}")
add_definitions(-DBUILD_AR25_FSI_DIALS=${BUILD_AR25_FSI_DIALS})
