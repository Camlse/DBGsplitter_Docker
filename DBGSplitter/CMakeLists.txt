project(DBGSplitter)

cmake_minimum_required(VERSION 3.1.0)

################################################################################
# The version number of YOUR tool.
################################################################################
SET (gatb-tool_VERSION_MAJOR 0)
SET (gatb-tool_VERSION_MINOR 0)
SET (gatb-tool_VERSION_PATCH 5)

IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# SUPPORTED KMER SIZES
################################################################################

# One can uncomment this line and set the wanted values
#set (KSIZE_LIST "32   64   96  128  160  192  224  256")
set (KSIZE_LIST "32   64   96  128")

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_INCLUDE_EXAMPLES  0)

# GATB CORE
include (GatbCore)

################################################################################
# TOOLS
################################################################################

ADD_SUBDIRECTORY(tools)

################################################################################
# Packaging
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "LBBE")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")  

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# Packaging the source ; we ignore unwanted files 
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
    "^${GATB_CORE_HOME}/.project"
    "^${GATB_CORE_HOME}/.gitignore"
    "^${GATB_CORE_HOME}/doc"
    "^${GATB_CORE_HOME}/DELIVERY.md"
)

# Packaging the binary ; we want to include some additional files
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/LICENCE                     DESTINATION .)
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/THIRDPARTIES.md  DESTINATION .)
INSTALL (PROGRAMS   ${CMAKE_CURRENT_SOURCE_DIR}/scripts/differential_analysis/analyse_diff_eve_cpx.R                     DESTINATION ./bin/)
INSTALL (PROGRAMS   ${CMAKE_CURRENT_SOURCE_DIR}/scripts/differential_analysis/run_diff_ana.R                     DESTINATION ./bin/)
INSTALL (PROGRAMS   ${CMAKE_CURRENT_SOURCE_DIR}/scripts/get_repeat_overlapping_nodes_using_ref/pipeline.sh                     DESTINATION ./bin/)
INSTALL (PROGRAMS   ${CMAKE_BINARY_DIR}/tools/run_complex_event_pipeline.py                     DESTINATION ./bin/)
INSTALL (DIRECTORY   ${PROJECT_SOURCE_DIR}/tests  DESTINATION .)

include (CPack)
