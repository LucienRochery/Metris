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

option(METRIS_INSTALL "Install Metris: disable if managing some of Metris deps" ON)

# Get the name of the binary directory
string(TOUPPER ${CMAKE_BINARY_DIR} BIN_DIR_NAME)
get_filename_component(BIN_DIR_NAME ${BIN_DIR_NAME} NAME_WE)

# Suppress nlopt warnings. NEW policy does not work with what nlopt is trying to do.
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27)
  set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
endif()

# Version >= 3.30 Boost related
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.30.0)
  # Consumers may use older CMake/Boost setups that expect the
  # legacy FindBoost/module behavior. Preserve legacy behavior to
  # remain compatible with outdated consumer configs.
  cmake_policy(SET CMP0167 OLD)
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
 

# These are system-wide install paths only:
set(METRIS_EXTERNAL_INCLUDE_DIRS "")
add_library(metris_deps INTERFACE)


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

# Include CMake utilities for RPATH, version checking, etc.
include(${PROJECT_SOURCE_DIR}/cmake/MetrisCMakeUtils.cmake)

# Set up RPATH configuration
metris_setup_rpath()

# Initialize dependency export list (cached) so helper functions can append to it
if(NOT DEFINED METRIS_CONFIG_DEPENDENCIES)
  set(METRIS_CONFIG_DEPENDENCIES "" CACHE INTERNAL "Metris dependencies list (internal)")
endif()
include(${PROJECT_SOURCE_DIR}/cmake/Dependencies.cmake)

# Get Boost libraries after dependencies are loaded (prefers modern targets)
metris_get_boost_libraries(BOOST_LIBS)

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
# 2 --- Main program
add_subdirectory(src)
#add_subdirectory(libs/libkdtree)


add_executable(metris "progs/main.cxx")
add_executable(egads2lite EXCLUDE_FROM_ALL "progs/EGADS2lite.cxx")

# Set up RPATH for metris executable
metris_setup_target_rpath(metris)

target_link_libraries(metris PUBLIC libMetris ${BOOST_LIBS})

#if(METRIS_USE_PETSC STREQUAL "True")
#  target_link_libraries(metris PUBLIC PkgConfig::PETSC ${MPI_C_LIBRARIES})
#  #target_link_libraries(metris PUBLIC PETSc)
#endif()

# Set up RPATH for egads2lite
metris_setup_target_rpath(egads2lite)

target_link_libraries(egads2lite PUBLIC libMetris ${BOOST_LIBS})

#For some reason, not specifying the full path does not work for EXISTS
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/progs/scatter_rh.cxx)
  add_executable(scatter_rh EXCLUDE_FROM_ALL "progs/scatter_rh.cxx")
  target_link_libraries(scatter_rh PUBLIC libMetris ${BOOST_LIBS})
  target_include_directories(scatter_rh PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisCXXFlags(scatter_rh PUBLIC)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/progs/extract_surf.cxx)
  add_executable(extract_surf EXCLUDE_FROM_ALL "progs/extract_surf.cxx")
  target_link_libraries(extract_surf PUBLIC libMetris ${BOOST_LIBS})
  target_include_directories(extract_surf PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisCXXFlags(extract_surf PUBLIC)
endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/progs/offset_coord.cxx)
  add_executable(offset_coord EXCLUDE_FROM_ALL "progs/offset_coord.cxx")
  target_link_libraries(offset_coord PUBLIC libMetris ${BOOST_LIBS})
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

# Optionally add consumer tests that exercise find_package(Metris)
option(METRIS_ENABLE_CONSUMER_TESTS "Enable consumer CMake tests under cmake/tests" ON)
if(METRIS_ENABLE_CONSUMER_TESTS)
  add_subdirectory(${PROJECT_SOURCE_DIR}/cmake/tests)
endif()





# --- Unit tests

# Fetch/set up cases directory.
if(DEFINED METRIS_CASES_DIR)
  message("-- Using provided METRIS_CASES_DIR = ${METRIS_CASES_DIR}")
elseif(DEFINED ENV{METRIS_CASES_DIR})
  set(METRIS_CASES_DIR $ENV{METRIS_CASES_DIR})
  message("-- Using env var METRIS_CASES_DIR = ${METRIS_CASES_DIR}")
else()
  set(METRIS_CASES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/examples")
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
  if(COMMAND ProcessorCount)
    ProcessorCount(NCORES)
  else()
    message(STATUS "ProcessorCount module not available, defaulting to 1 core")
    set(NCORES 1)
  endif()
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
                                                         METRIS_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
                                                         )
                                                        # METRIS_REGRESSION_DIR="${METRIS_REGRESSION_DIR}"
                                                         
  target_link_libraries(${test_target_build} PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_link_libraries(${test_target_build} PUBLIC libMetris 
                                                    metris_deps
                                                    ${BOOST_LIBS}
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
if(METRIS_INSTALL)
  install(TARGETS libMetris metris_deps EXPORT libMetrisTargets
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

  # Also export the targets to a file in the build tree so consumers can
  # use the build tree (via METRIS_DIR) without requiring an install step.
  # This creates `${CMAKE_CURRENT_BINARY_DIR}/libMetrisTargets.cmake`.
  export(EXPORT libMetrisTargets
         FILE "${CMAKE_CURRENT_BINARY_DIR}/libMetrisTargets.cmake"
         NAMESPACE Metris::)

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfigVersion.cmake"
      #VERSION ${PROJECT_VERSION}
      COMPATIBILITY AnyNewerVersion
  )



  set(METRIS_DEPS_EXPORT_CODE "")
  foreach(dep_entry ${METRIS_CONFIG_DEPENDENCIES})
    message( "dep_entry = ${dep_entry}")
    # Support both old pipe-separated format and new semicolon-separated format
    if(dep_entry MATCHES "\\|")
      # Old format: "TYPE|NAME|COMPONENTS"
      string(REPLACE "|" ";" dep_parts "${dep_entry}")
    else()
      # New format: already a list with semicolons
      set(dep_parts "${dep_entry}")
    endif()
    list(LENGTH dep_parts dep_parts_len)
    if(dep_parts_len GREATER_EQUAL 1)
      list(GET dep_parts 0 dep_type)
    endif()
    if(dep_parts_len GREATER_EQUAL 2)
      list(GET dep_parts 1 dep_name)
    endif()
    if(dep_parts_len GREATER_EQUAL 3)
      list(GET dep_parts 2 dep_components)
    else()
      set(dep_components "")
    endif()
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
      list(LENGTH dep_components comp_len)
      if(comp_len GREATER_EQUAL 1)
        list(GET dep_components 0 DEP_LIBRARIES)
      endif()
      if(comp_len GREATER_EQUAL 2)
        list(GET dep_components 1 DEP_INCLUDE_DIRS)
      endif()

      # For provided dependencies, normally we export variables listing libraries and include dirs.
      # However some header-only or target-based provided deps (notably Eigen) are referenced by
      # the exported targets as imported targets (e.g. Eigen3::Eigen). In that case consumers
      # need a corresponding find_dependency() call in the generated config so CMake defines
      # the imported target when the package is consumed.
      if(dep_name STREQUAL "Eigen3")
        # Ensure consumers will find Eigen and get Eigen3::Eigen target available
        string(APPEND METRIS_DEPS_EXPORT_CODE "find_dependency(Eigen3 REQUIRED)\\n")
      endif()

      # Still export the raw lists for compatibility with older consumers that use variables
      string(APPEND METRIS_DEPS_EXPORT_CODE "list(APPEND METRIS_LIBRARIES ${DEP_LIBRARIES})\\n")
      string(APPEND METRIS_DEPS_EXPORT_CODE "list(APPEND METRIS_INCLUDE_DIRS ${DEP_INCLUDE_DIRS})\\n")
    else()
      message(FATAL_ERROR "Unknown dependency type: ${dep_type} for ${dep_name}")
    endif()
  endforeach()

  #message(FATAL_ERROR "METRIS_DEPS_EXPORT_CODE = ${METRIS_DEPS_EXPORT_CODE}")

  # Robustness: if the generated targets file references an imported target
  # like Eigen3::Eigen but the dependency list did not record a find_package
  # entry (e.g. due to 'provided' handling or formatting), ensure we still
  # emit a find_dependency for Eigen so consumers of the installed package
  # won't fail trying to resolve the imported target.
  set(_targets_file "${CMAKE_CURRENT_BINARY_DIR}/libMetrisTargets.cmake")
  if(EXISTS "${_targets_file}")
    file(READ "${_targets_file}" _targets_contents)
    string(FIND "${_targets_contents}" "Eigen3::Eigen" _found_eigen)
    if(NOT _found_eigen EQUAL -1)
      string(FIND "${METRIS_DEPS_EXPORT_CODE}" "find_dependency(Eigen3" _has_eigen_find)
      if(_has_eigen_find EQUAL -1)
        message(STATUS "libMetrisTargets references Eigen3::Eigen but no find_dependency(Eigen3) was recorded; adding fallback find_dependency(Eigen3 REQUIRED) in package config.")
        string(APPEND METRIS_DEPS_EXPORT_CODE "find_dependency(Eigen3 REQUIRED)\\n")
      endif()
    endif()
    unset(_targets_contents)
  endif()

    # Deduplicate dependency lines while preserving order and removing empty lines
    string(REPLACE "\n" ";" _deps_list_raw "${METRIS_DEPS_EXPORT_CODE}")
    set(_deps_list_clean "")
    foreach(_entry IN LISTS _deps_list_raw)
      if(NOT _entry STREQUAL "")
        list(APPEND _deps_list_clean "${_entry}")
      endif()
    endforeach()
    list(REMOVE_DUPLICATES _deps_list_clean)
    # Rebuild METRIS_DEPS_EXPORT_CODE from cleaned list
    set(METRIS_DEPS_EXPORT_CODE "")
    foreach(_entry IN LISTS _deps_list_clean)
      string(APPEND METRIS_DEPS_EXPORT_CODE "${_entry}\n")
    endforeach()

    configure_package_config_file(
      "${PROJECT_SOURCE_DIR}/cmake/MetrisConfig.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfig.cmake"
      INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Metris
  )
    # Post-process the generated package config produced by
    # configure_package_config_file(). We intentionally avoid writing the
    # template directly (which would leave @PACKAGE_INIT@ unexpanded).
    # Instead, read the generated config and replace the @METRIS_DEPS_EXPORT_CODE@
    # marker only if present.
    set(_generated_cfg "${CMAKE_CURRENT_BINARY_DIR}/MetrisConfig.cmake")
    if(EXISTS "${_generated_cfg}")
      file(READ "${_generated_cfg}" _cfg_contents)
      string(FIND "${_cfg_contents}" "@METRIS_DEPS_EXPORT_CODE@" _has_token)
      if(NOT _has_token EQUAL -1)
        string(REPLACE "@METRIS_DEPS_EXPORT_CODE@" "${METRIS_DEPS_EXPORT_CODE}" _cfg_contents "${_cfg_contents}")
        file(WRITE "${_generated_cfg}" "${_cfg_contents}")
        message(STATUS "Injected METRIS_DEPS_EXPORT_CODE into ${_generated_cfg}")
      endif()
      unset(_cfg_contents)
    else()
      message(WARNING "Expected generated config ${_generated_cfg} not found; no post-processing done")
    endif()

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
          "METRIS_EXTERNAL_INCLUDE_DIRS = ${METRIS_EXTERNAL_INCLUDE_DIRS}")
  message("=================================================================")
endif()