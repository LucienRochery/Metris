# Parameters:
# - METRIS_BUILD_CORES : build cores for unit test deployment, default n-1
# - METRIS_MAX_DEG : maximum degree to support, default 2


# FetchContent_MakeAvailable is latest constraint:
cmake_minimum_required(VERSION 3.14.0)
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# 20 required for concepts
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF) # Otherwise get -std=gnu++17 messes with intel



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

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  cmake_policy(SET CMP0135 NEW)
endif()


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Install destination" FORCE)
  message(STATUS "Using default install prefix ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "Using provided install prefix ${CMAKE_INSTALL_PREFIX}")
endif()

# This must come before project or we get an infinite loop.
include(cmake/Compiler.cmake)

project(Metris VERSION 0.5.0)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

include(cmake/BuildType.cmake)
include(cmake/MetrisFlags.cmake)
 

set(METRIS_DEPS_LIBRARIES "")
# INTERNAL_INCLUDE_DIRS paths are all relative to ${PROJECT_SOURCE_DIR}. The directory structure is replicated
# in install/include 
set(METRIS_INTERNAL_INCLUDE_DIRS "")
# These are system-wide install paths only:
set(METRIS_EXTERNAL_INCLUDE_DIRS "")


find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} config --get remote.origin.url
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
                OUTPUT_VARIABLE GITURL
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${GIT_EXECUTABLE} branch --show-current
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
                OUTPUT_VARIABLE GIT_BRANCH
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)

add_compile_definitions(METRIS_GIT_URL="${GITURL}")
add_compile_definitions(METRIS_GIT_COMMIT_HASH="${GIT_COMMIT_HASH}")
message("Metris using METRIS_GIT_URL = ${GITURL}")
message("In branch: ${GIT_BRANCH}")
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



add_subdirectory(libs)
add_subdirectory(SANS)
# 2 --- Main program
add_subdirectory(src)
#add_subdirectory(libs/libkdtree)


add_executable(metris "progs/main.cxx")
add_executable(egads2lite EXCLUDE_FROM_ALL "progs/EGADS2lite.cxx")

set_target_properties(metris PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


target_link_libraries(metris PUBLIC libMetris 
                                    Boost::headers # yes this is weird
                                    ${Boost_LIBRARIES}
                                    )

#if(METRIS_USE_PETSC STREQUAL "True")
#  target_link_libraries(metris PUBLIC PkgConfig::PETSC ${MPI_C_LIBRARIES})
#  #target_link_libraries(metris PUBLIC PETSc)
#endif()

target_link_libraries(egads2lite PUBLIC libMetris 
                                        Boost::headers # yes this is weird
                                        ${Boost_LIBRARIES}
                                        )

#For some reason, not specifying the full path does not work for EXISTS
if(EXISTS ${CMAKE_SOURCE_DIR}/progs/scatter_rh.cxx)
  add_executable(scatter_rh EXCLUDE_FROM_ALL "progs/scatter_rh.cxx")
  target_link_libraries(scatter_rh PUBLIC libMetris 
                                          Boost::headers # yes this is weird
                                          ${Boost_LIBRARIES}
                                          )
  target_include_directories(scatter_rh PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisCXXFlags(scatter_rh PUBLIC)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/progs/extract_surf.cxx)
  add_executable(extract_surf EXCLUDE_FROM_ALL "progs/extract_surf.cxx")
  target_link_libraries(extract_surf PUBLIC libMetris 
                                            Boost::headers # yes this is weird
                                            ${Boost_LIBRARIES}
                                            )
  target_include_directories(extract_surf PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisCXXFlags(extract_surf PUBLIC)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/progs/offset_coord.cxx)
  add_executable(offset_coord EXCLUDE_FROM_ALL "progs/offset_coord.cxx")
  target_link_libraries(offset_coord PUBLIC libMetris 
                                            Boost::headers # yes this is weird
                                            ${Boost_LIBRARIES}
                                            )
  target_include_directories(offset_coord PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisCXXFlags(offset_coord PUBLIC)
endif()

# Note on includes: the src lib target (libMetris) already includes publicly everything that we need
# The main.cxx itself only needs main.hxx additionally. 
target_include_directories(metris PRIVATE "progs/main.hxx"
                                  PUBLIC ${METRIS_DEPS_INCLUDE_DIRS}
                                  )
setMetrisCXXFlags(metris PUBLIC)

target_include_directories(egads2lite PRIVATE  "."
                           ${ESP_ROOT}/include)
setMetrisCXXFlags(egads2lite PUBLIC)





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
set(CTEST_OUTPUT_ON_FAILURE 1)

file (GLOB TEST_SOURCES CONFIGURE_DEPENDS "bunit/*.cxx")
file (GLOB_RECURSE TEST_REG_SOURCES CONFIGURE_DEPENDS "bunit/regression/*.cxx")


# -- Boost/test
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# -- For building in FindUnitTestMeshes and setup_unit_build_and_deploy
if(NOT DEFINED METRIS_BUILD_CORES)
  include(ProcessorCount)
  ProcessorCount(NCORES)
  if(NCORES EQUAL 0)
    set(NCORES 1)
  endif()
  math(EXPR METRIS_BUILD_CORES "${NCORES} - 1")
  if(METRIS_BUILD_CORES LESS 1)
    set(METRIS_BUILD_CORES 1)
  endif()
  set(METRIS_BUILD_CORES ${METRIS_BUILD_CORES} CACHE STRING "Number of cores to use when building Metris for unit test mesh generation")
endif()

# -- Get the required meshes by reading bunit/test_ files
# Uses TEST_SOURCES 
include(cmake/FindUnitTestMeshes.cmake)
# We can now concatenate TEST_SOURCES and TEST_REG_SOURCES
set(TEST_SOURCES ${TEST_SOURCES} ${TEST_REG_SOURCES})

#add_custom_command(OUTPUT ${METRIS_UNIT_TEST_MESHES}
#                   DEPENDS metris
#                   COMMAND ${CMAKE_COMMAND} -E env 
#                           METRIS_CASES_DIR=${METRIS_CASES_DIR}
#                           ${METRIS_CASES_DIR}/gen_cases.sh)
#add_custom_target(deploy_cases DEPENDS ${METRIS_UNIT_TEST_MESHES})

set(METRIS_TARGET_UNIT_TESTS_BUILD "")
set(METRIS_TARGET_REGRESSION_TESTS_BUILD "")

# Dummy test that runs the required builds
add_test(
  NAME setup_unit_build_and_deploy
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target deploy_cases metris_unit_tests_build --parallel ${METRIS_BUILD_CORES}
)
set_tests_properties(setup_unit_build_and_deploy PROPERTIES
                     LABELS "unit"
                     FIXTURES_SETUP "unit_build_and_deploy"
                     )

add_test(
  NAME setup_regression_build_and_deploy
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target metris_regression_tests_build --parallel ${METRIS_BUILD_CORES}
)
set_tests_properties(setup_regression_build_and_deploy PROPERTIES
                     LABELS "regression"
                     FIXTURES_SETUP "regression_build_and_deploy"
                     )

foreach(test_file ${TEST_SOURCES})
  get_filename_component(test_target ${test_file} NAME_WLE) # Get name e.g. unit_dbg1 out of $PATHTOFILE/unit_dbg1.cxx
  set(test_target_build "${test_target}_build")
  add_executable(${test_target_build} EXCLUDE_FROM_ALL ${test_file})
  set_target_properties(${test_target_build} PROPERTIES OUTPUT_NAME ${test_target})
  setMetrisCXXFlags(${test_target_build} PUBLIC)
  target_compile_definitions(${test_target_build} PUBLIC METRIS_CASES_DIR="${METRIS_CASES_DIR}"
                                                         METRIS_ROOT_DIR="${CMAKE_SOURCE_DIR}"
                                                         )
                                                        # METRIS_REGRESSION_DIR="${METRIS_REGRESSION_DIR}"
                                                         
  target_link_libraries(${test_target_build} PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_link_libraries(${test_target_build} PUBLIC libMetris 
                                                    ${METRIS_DEPS_LIBRARIES}
                                                    Boost::headers # yes this is weird
                                                    ${Boost_LIBRARIES}
                                                    )
  #target_include_directories(${test_target} PUBLIC "."  ${PETSC_INCLUDE}  ${ABSL_INCLUDE_DIRS})


  #add_custom_target(${test_target}
  #    DEPENDS ${test_target_build}
  #    COMMAND $<TARGET_FILE:${test_target_build}> 
  #)

  

  # Note: CMake uses regexp, not  bash-style globs. 
  # Hence _* would match to any number of '_', but not to any number of any characters...
  if(test_target MATCHES "^test_.*")
    list(APPEND METRIS_TARGET_UNIT_TESTS_BUILD ${test_target_build})
    add_test(NAME ${test_target}
             COMMAND $<TARGET_FILE:${test_target_build}>
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            )

    # Doing add_test(toto) doesn't create a toto target 
    # that we can add deps to. Instead, we should use fixtures.
    set_tests_properties(${test_target} PROPERTIES
                         FIXTURES_REQUIRED "unit_build_and_deploy"
                         LABELS "unit"
                         )
  elseif(test_target MATCHES "^reg_.*")
    list(APPEND METRIS_TARGET_REGRESSION_TESTS_BUILD ${test_target_build})
    add_test(NAME ${test_target}
             COMMAND $<TARGET_FILE:${test_target_build}>
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            )
    set_tests_properties(${test_target} PROPERTIES
                         FIXTURES_REQUIRED "regression_build_and_deploy"
                         LABELS "regression"
                         )
  endif()

  ## Ensure that foo_build is built before foo
  #add_dependencies(${test_target}  deploy_cases)
endforeach()

# A target to manually compile the tests, no longer necessary
add_custom_target(metris_unit_tests_build
                  DEPENDS ${METRIS_TARGET_UNIT_TESTS_BUILD}
                  )
add_custom_target(metris_regression_tests_build
                  DEPENDS ${METRIS_TARGET_REGRESSION_TESTS_BUILD}
                  )



# --- Install libMetris

install(TARGETS libMetris libmeshb lplib3 libkdtree alglib libSANS_Metris EXPORT libMetrisTargets
  RUNTIME  DESTINATION   ${CMAKE_INSTALL_BINDIR}
           COMPONENT     runtime OPTIONAL
  LIBRARY  DESTINATION   ${CMAKE_INSTALL_LIBDIR}
           COMPONENT     runtime OPTIONAL
           #NAMELINK_SKIP
  #ARCHIVE  DESTINATION   ${CMAKE_INSTALL_LIBDIR}
           #COMPONENT     Y
  INCLUDES DESTINATION   ${CMAKE_INSTALL_INCLUDEDIR})


# Install all headers. We don't have any specific public headers, users can do anything (as we do in unit tests at times, for instance)
install(
  FILES ${PROJECT_SOURCE_DIR}/src/Metris.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src/
)
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src/
  FILES_MATCHING
    PATTERN "*.hxx"
    PATTERN ".*" EXCLUDE
)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/SANS/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SANS/
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN ".*" EXCLUDE
)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/libs/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libs/
  FILES_MATCHING
    PATTERN "*.h*"
    PATTERN ".*" EXCLUDE
)

install(EXPORT libMetrisTargets
  FILE libMetrisTargets.cmake
  NAMESPACE Metris::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Metris
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfigVersion.cmake"
    #VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)



set(METRIS_DEPS_EXPORT_CODE "")
foreach(dep_entry ${METRIS_CONFIG_DEPENDENCIES})
  message( "dep_entry = ${dep_entry}")
  string(REPLACE "|" ";" dep_parts "${dep_entry}")
  list(GET dep_parts 0 dep_type)
  list(GET dep_parts 1 dep_name)
  list(GET dep_parts 2 dep_components)
  message( "dep_type = ${dep_type}, dep_name = ${dep_name}, dep_components = ${dep_components}")
  
  if(dep_type STREQUAL "find_package")
    if(dep_components)
      string(APPEND METRIS_DEPS_EXPORT_CODE "find_dependency(${dep_name} COMPONENTS ${dep_components})\n")
    else()
      string(APPEND METRIS_DEPS_EXPORT_CODE "find_dependency(${dep_name})\n")
    endif()
  elseif(dep_type STREQUAL "FetchContent")
    # For FetchContent dependencies, we need a different approach
    string(APPEND METRIS_DEPS_EXPORT_CODE "# ${dep_name} was built via FetchContent and is included in the exported targets\n")
  elseif(dep_type STREQUAL "provided")
    # Parse the space separated components into two items:
    string(REPLACE " " ";" dep_components "${dep_components}")
    list(GET dep_components 0 DEP_LIBRARIES)
    list(GET dep_components 1 DEP_INCLUDE_DIRS)

    # For provided dependencies, we can just export the variables
    string(APPEND METRIS_DEPS_EXPORT_CODE "list(APPEND METRIS_LIBRARIES ${DEP_LIBRARIES})\n")
    string(APPEND METRIS_DEPS_EXPORT_CODE "list(APPEND METRIS_INCLUDE_DIRS ${DEP_INCLUDE_DIRS})\n")
  else()
    message(FATAL_ERROR "Unknown dependency type: ${dep_type} for ${dep_name}")
  endif()
endforeach()

#message(FATAL_ERROR "METRIS_DEPS_EXPORT_CODE = ${METRIS_DEPS_EXPORT_CODE}")

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/MetrisConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Metris
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfigVersion.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/FindNLopt.cmake"  # Install your custom finder
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Metris
)

message("=================================================================")
message(STATUS "Metris configure script end, using "
        "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}, \n"
        "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}, \n"
        "CMAKE_INSTALL_INCLUDEDIR = ${CMAKE_INSTALL_INCLUDEDIR}, \n"
        "CMAKE_INSTALL_BINDIR = ${CMAKE_INSTALL_BINDIR}, \n"
        "METRIS_DEPS_LIBRARIES = ${METRIS_DEPS_LIBRARIES}, \n"
        "METRIS_INTERNAL_INCLUDE_DIRS = ${METRIS_INTERNAL_INCLUDE_DIRS}, \n"
        "METRIS_EXTERNAL_INCLUDE_DIRS = ${METRIS_EXTERNAL_INCLUDE_DIRS}")
message("=================================================================")
