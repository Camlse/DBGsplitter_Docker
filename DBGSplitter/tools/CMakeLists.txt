################################################################################
#  TOOLS
################################################################################

MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET (dirlist "")
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
          LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    SET(${result} ${dirlist})
ENDMACRO()

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# We add the gatb-core include directory
include_directories (${gatb-core-includes})

# We add the pstream include directory
set (PSTREAMS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/pstreams/)
include_directories (${PSTREAMS_SOURCE_DIR})

# We add the whereami include directory
set (WHEREAMI_SOURCE_DIR ${PROJECT_SOURCE_DIR}/whereami/src)
include_directories (${WHEREAMI_SOURCE_DIR})

# We add the BOOST include directory
set(Boost_USE_STATIC_LIBS   ON) #This makes boost plug the libraries we specify as static libraries (not dynamic ones - the user does not even need to have boost installed in their pc)
FIND_PACKAGE(Boost COMPONENTS filesystem system serialization timer chrono REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})
link_directories(${Boost_LIBRARY_DIR})

# We get the list of subdirectories names
#SUBDIRLIST(PROGRAMS ${PROJECT_SOURCE_DIR}/tools)
#MESSAGE("PROGRAMS VAR: ${PROGRAMS}")

# We create a target for each
#FOREACH (program ${PROGRAMS})
set (COMMON_SOURCE_DIR ${PROJECT_SOURCE_DIR}/tools/DBGSplitter/src)
include_directories (${COMMON_SOURCE_DIR})
file (GLOB_RECURSE  ProjectFiles  ${COMMON_SOURCE_DIR}/*.cpp)
file (GLOB_RECURSE  WhereamiFiles  ${WHEREAMI_SOURCE_DIR}/*.c)
add_executable(step1 ${PROJECT_SOURCE_DIR}/tools/DBGSplitter/step1/main.cpp ${ProjectFiles} ${WhereamiFiles})
target_link_libraries(step1 ${gatb-core-libraries} ${Boost_LIBRARIES})
add_executable(step2 ${PROJECT_SOURCE_DIR}/tools/DBGSplitter/step2/main.cpp ${ProjectFiles} ${WhereamiFiles})
target_link_libraries(step2 ${gatb-core-libraries} ${Boost_LIBRARIES})
add_executable(step3 ${PROJECT_SOURCE_DIR}/tools/DBGSplitter/step3/main.cpp ${ProjectFiles} ${WhereamiFiles})
target_link_libraries(step3 ${gatb-core-libraries} ${Boost_LIBRARIES})
add_executable(step4 ${PROJECT_SOURCE_DIR}/tools/DBGSplitter/step4/main.cpp ${ProjectFiles} ${WhereamiFiles})
target_link_libraries(step4 ${gatb-core-libraries} ${Boost_LIBRARIES})

#ENDFOREACH (program)


#copy some files to the tools folder
configure_file(${PROJECT_SOURCE_DIR}/scripts/get_repeat_overlapping_nodes_using_ref/pipeline.sh ${PROJECT_BINARY_DIR}/tools/ COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/scripts/differential_analysis/analyse_diff_eve_cpx.R ${PROJECT_BINARY_DIR}/tools/ COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/scripts/differential_analysis/run_diff_ana.R ${PROJECT_BINARY_DIR}/tools/ COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/tools/DBGSplitter/run_complex_event_pipeline.py ${PROJECT_BINARY_DIR}/tools/)


################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS step1 step2 step3 step4 DESTINATION bin)
