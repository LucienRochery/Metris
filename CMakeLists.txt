# FetchContent_MakeAvailable is latest constraint:
cmake_minimum_required(VERSION 3.14.0)
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set(CMAKE_INCLUDE_CURRENT_DIR ON)


# 20 required for concepts
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF) # Otherwise get -std=gnu++17 messes with intel

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Get the name of the binary directory
STRING( TOUPPER ${CMAKE_BINARY_DIR} BIN_DIR_NAME )
GET_FILENAME_COMPONENT( BIN_DIR_NAME ${BIN_DIR_NAME} NAME_WE )

# Suppress nlopt warnings. NEW policy does not work with what nlopt is trying to do.
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27)
  set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
endif()

# Version >= 3.30 Boost related
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.30.0)
  cmake_policy(SET CMP0167 NEW)
endif()

# This must come before project or we get an infinite loop.
include(cmake/Compiler.cmake)

project(Metris)

include(cmake/BuildType.cmake)
include(cmake/MetrisFlags.cmake)
 

set(METRIS_DEPS_LIBRARIES "")
set(METRIS_DEPS_INCLUDE_DIRS ${PROJECT_SOURCE_DIR})


find_package(Git)
execute_process(COMMAND "${GIT_EXECUTABLE}" config --get remote.origin.url
                WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" 
                RESULT_VARIABLE dum
                OUTPUT_VARIABLE GITURL
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)

add_compile_definitions(METRIS_GIT_URL="${GITURL}")
message("Metris using METRIS_GIT_URL = ${GITURL}")
if(GITURL MATCHES "github.com")
  set(LIGHT_BUILD ON)
  message("-- Disabled codegen (max degree 2 supported).")
  set(REQ_CODEGEN OFF)
  set(CODEGEN_FILES ${PROJECT_SOURCE_DIR}/src/codegen_ccoef.hxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02_d_coord.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02_d_pt.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.02.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.02_d_coord.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef_d.hxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_lag2bez.hxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_lagrange.hxx)
  foreach(f ${CODEGEN_FILES})
    if(NOT EXISTS "${f}")
      message(FATAL_ERROR "Codegen file ${f} does not exist")
    endif()
  endforeach()
else()
  message("-- Metris with codegen enabled ")
  set(LIGHT_BUILD OFF)
  set(REQ_CODEGEN ON)
endif()


add_compile_definitions(SURREAL_LAZY)

if(USE_TRACELIBS)
  message(FATAL_ERROR "-- Compiling with Boost stacktrace libraries")
else()
  add_compile_definitions(NO_BOOST_EXCEPT)
endif()

set(CMAKE_MACOSX_RPATH 1)

include(${PROJECT_SOURCE_DIR}/cmake/Dependencies.cmake)


#find_package(MPI REQUIRED)
#set(CMAKE_CXX_COMPILER mpicxx)
#message("Have CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
#message("Using CXX = ${CMAKE_CXX_COMPILER}")

get_filename_component(SHORT_COMPILER_NAME ${CMAKE_C_COMPILER} NAME_WE)



if(METRIS_MAX_DEG)
  add_compile_definitions(METRIS_MAX_DEG=${METRIS_MAX_DEG})
else()
  set(METRIS_MAX_DEG 2)
  add_compile_definitions(METRIS_MAX_DEG=${METRIS_MAX_DEG})
endif()

if(ARRAY2D_EXPENSIVE)
  add_compile_definitions(ARRAY2D_EXPENSIVE)
endif()

if(USE_MULTIPRECISION)
  add_compile_definitions(USE_MULTIPRECISION)
endif()

if(USE_ABSL)
  add_compile_definitions(USE_ABSL)
endif()



#file (GLOB_RECURSE SRC_HEADERS CONFIGURE_DEPENDS "src/*.cxx")
#file (GLOB_RECURSE SRC_HEADERS CONFIGURE_DEPENDS "src/*.hxx")
#file (GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "bunit/*.cxx")
#file (GLOB_RECURSE CGEN_SOURCES CONFIGURE_DEPENDS "tools/*.cxx")

#set (SRC_INCLUDE_DIRS "")
#foreach (h ${SRC_HEADERS})
#  get_filename_component(d ${h} PATH)
#  list (APPEND SRC_INCLUDE_DIRS ${d})
#endforeach()
#list(REMOVE_DUPLICATES SRC_INCLUDE_DIRS)

#find_package(LAPACK REQUIRED)
#find_package(PkgConfig REQUIRED)
#find_package(NLopt REQUIRED)
#find_package(Boost REQUIRED)
#find_package(LPlib REQUIRED)


# 2 --- Main program
add_subdirectory(libs)
add_subdirectory(src)
add_subdirectory(SANS)
#add_subdirectory(libs/libkdtree)


add_executable(metris "progs/main.cxx")
add_executable(egads2lite "progs/EGADS2lite.cxx")

set_target_properties(metris PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)



target_link_libraries(metris PUBLIC libMetris 
                                    SANS 
                                    ${ABSL_LIBRARIES}
                                    ${LAPACK_LIBRARIES}
                                    ${EGADS_LIBRARIES} 
                                    ${EGADSLITE_LIBRARIES}
                                    #nlopt
#                                   LP.3
                                    Boost::headers # yes this is weird
                                    ${Boost_LIBRARIES}
                                    ${Boost_PROGRAM_OPTIONS_LIBRARY}
                                    ${PETSC_LIBRARIES}
                                    ${MPI_LIBRARIES}
                                    ${BOOST_TRACELIBS}
                                    ${METRIS_DEPS_LIBRARIES}
                                    )

#if(USE_PETSC STREQUAL "True")
#  target_link_libraries(metris PUBLIC PkgConfig::PETSC ${MPI_C_LIBRARIES})
#  #target_link_libraries(metris PUBLIC PETSc)
#endif()

target_link_libraries(egads2lite PUBLIC libMetris 
                                        SANS 
                                        ${ABSL_LIBRARIES}
                                        ${LAPACK_LIBRARIES}
                                        ${EGADS_LIBRARIES} 
                                        ${EGADSLITE_LIBRARIES}
                                        #nlopt
#                                       LP.3
                                        Boost::headers # yes this is weird
                                        ${Boost_LIBRARIES}
                                        ${Boost_PROGRAM_OPTIONS_LIBRARY}
                                        ${PETSC_LIBRARIES}
                                        ${MPI_LIBRARIES}
                                        ${BOOST_TRACELIBS})

#For some reason, not specifying the full path does not work for EXISTS
if(EXISTS ${CMAKE_SOURCE_DIR}/progs/scatter_rh.cxx)
  add_executable(scatter_rh "progs/scatter_rh.cxx")
  target_link_libraries(scatter_rh PUBLIC libMetris 
                                          #libegads libegadslite
                                          #egads egadslite
                                          ${EGADS_LIBRARIES} 
                                          ${EGADSLITE_LIBRARIES}
                                          Boost::headers # yes this is weird
                                          ${Boost_LIBRARIES}
                                          ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(scatter_rh PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisFlags(scatter_rh PUBLIC)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/progs/extract_surf.cxx)
  add_executable(extract_surf "progs/extract_surf.cxx")
  target_link_libraries(extract_surf PUBLIC libMetris 
                                          #libegads libegadslite
                                          #egads egadslite
                                          ${EGADS_LIBRARIES} ${EGADSLITE_LIBRARIES}
                                          Boost::headers # yes this is weird
                                          ${Boost_LIBRARIES}
                                          ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(extract_surf PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisFlags(extract_surf PUBLIC)
endif()

# Note on includes: the src lib target (libMetris) already includes publicly everything that we need
# The main.cxx itself only needs main.hxx additionally. 
target_include_directories(metris PRIVATE "progs/main.hxx"
                                  PUBLIC ${METRIS_DEPS_INCLUDE_DIRS}
                                  )
setMetrisFlags(metris PUBLIC)

target_include_directories(egads2lite PRIVATE  "."
                           ${ESP_ROOT}/include)
setMetrisFlags(egads2lite PUBLIC)


# --- Unit tests

# Fetch/set up cases directory.
if(DEFINED METRIS_CASES_DIR)
  message("-- Using provided METRIS_CASES_DIR = ${METRIS_CASES_DIR}")
elseif(DEFINED ENV{METRIS_CASES_DIR})
  set(METRIS_CASES_DIR $ENV{METRIS_CASES_DIR})
  message("-- Using env var METRIS_CASES_DIR = ${METRIS_CASES_DIR}")
else()
  set(METRIS_CASES_DIR "${CMAKE_SOURCE_DIR}/examples")
  message("-- Using default METRIS_CASES_DIR = root/cases = ${METRIS_CASES_DIR}")
endif()
set(METRIS_CASES_DIR ${METRIS_CASES_DIR} CACHE PATH "Directory where Metris cases are found.")
message("-- Created directory ${METRIS_CASES_DIR} (unless it already existed)")
make_directory(${METRIS_CASES_DIR})

enable_testing()

file (GLOB TEST_SOURCES CONFIGURE_DEPENDS "bunit/*.cxx")

# -- Boost/test
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# -- Get the required meshes
set(METRIS_UNIT_TEST_MESHES "")
foreach(test_file ${TEST_SOURCES})
  file(READ ${test_file} CMAKE_FILE_CONTENT)

  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+)[\" ]" MATCHES_NOEXT ${CMAKE_FILE_CONTENT})
  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.mesh)" MATCHES_PMESH ${CMAKE_FILE_CONTENT})
  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.meshb)" MATCHES_PMESHB ${CMAKE_FILE_CONTENT})

  foreach(match IN LISTS MATCHES_NOEXT)
      string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+)[\" ]" "\\1" meshname "${match}")
      set(meshname "${meshname}.meshb")
      set(meshname "${METRIS_CASES_DIR}/${meshname}")
      if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
        list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
      endif()
  endforeach()

  foreach(match IN LISTS MATCHES_PMESH)
      string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.mesh)" "\\1" meshname "${match}")
      set(meshname "${METRIS_CASES_DIR}/${meshname}")
      if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
        list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
      endif()
  endforeach()

  foreach(match IN LISTS MATCHES_PMESHB)
      string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.meshb)" "\\1" meshname "${match}")
      set(meshname "${METRIS_CASES_DIR}/${meshname}")
      if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
        list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
      endif()
  endforeach()
endforeach()


add_custom_command(OUTPUT ${METRIS_UNIT_TEST_MESHES}
                   DEPENDS metris
                   COMMAND ${METRIS_CASES_DIR}/gen_cases.sh)
add_custom_target(deploy_cases DEPENDS ${METRIS_UNIT_TEST_MESHES})

foreach(test_file ${TEST_SOURCES})
  get_filename_component(test_target ${test_file} NAME_WE) # Get name e.g. unit_dbg1 out of $PATHTOFILE/unit_dbg1.cxx
  set(test_target_build "${test_target}_build")
  add_executable(${test_target_build} ${test_file})
  set_target_properties(${test_target_build} PROPERTIES OUTPUT_NAME ${test_target})
  setMetrisFlags(${test_target_build} PUBLIC)
  target_compile_definitions(${test_target_build} PUBLIC METRIS_CASES_DIR="${METRIS_CASES_DIR}")
  target_link_libraries(${test_target_build} PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_link_libraries(${test_target_build} PUBLIC libMetris 
                                              SANS 
                                              ${ABSL_LIBRARIES}
                                              ${LAPACK_LIBRARIES}
                                              ${EGADS_LIBRARIES} 
                                              ${EGADSLITE_LIBRARIES}
                                              Boost::headers # yes this is weird
                                              ${Boost_LIBRARIES}
                                              ${Boost_PROGRAM_OPTIONS_LIBRARY}
                                              ${PETSC_LIBRARIES}
                                              ${MPI_LIBRARIES}
                                              ${BOOST_TRACELIBS})
  #target_include_directories(${test_target} PUBLIC "."  ${PETSC_INCLUDE}  ${ABSL_INCLUDE_DIRS})


  add_custom_target(${test_target}
      DEPENDS ${test_target_build}
      COMMAND ${test_target}
  )

  # Ensure that foo_build is built before foo
  add_dependencies(${test_target} ${test_target_build} deploy_cases)
endforeach()
    

install(TARGETS metris DESTINATION ${CMAKE_BINARY_DIR})

