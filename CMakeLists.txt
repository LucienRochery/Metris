cmake_minimum_required(VERSION 3.22)
project(Metris)
set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")


# 20 required for concepts
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

include(MetrisFlags)

find_package(Git)
execute_process(COMMAND "${GIT_EXECUTABLE}" config --get remote.origin.url
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                RESULT_VARIABLE dum
                OUTPUT_VARIABLE GITURL
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)

if(GITURL MATCHES "github.com")
  set(LIGHT_BUILD ON)
  message(WARNING "Disabled codegen!")
  set(REQ_CODEGEN OFF)
  set(DO_UNITTESTS OFF)
  set(CODEGEN_FILES ${PROJECT_SOURCE_DIR}/src/codegen_ccoef.hxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02_d.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef2.02_d_pt.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.02.cxx
                    ${PROJECT_SOURCE_DIR}/src/codegen_ccoef3.02_d.cxx
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
  set(DO_UNITTESTS ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "No build type defined: defaulting to Release")
  set(CMAKE_BUILD_TYPE Release)
endif()


add_compile_definitions(SURREAL_LAZY)

add_compile_definitions(NO_BOOST_EXCEPT)
#if(USE_TRACELIBS)
#  message("-- Compiling with Boost stacktrace libraries")
#else()
#  add_compile_definitions(NO_BOOST_EXCEPT)
#endif()

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

find_package(LAPACK REQUIRED)
find_package(PkgConfig REQUIRED)
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


#message("MPI_LIBRARIES = ${MPI_LIBRARIES}")

target_link_libraries(metris PUBLIC libMetris 
                                    libmeshb
                                    SANS 
                                    ${ABSL_LIBRARIES}
                                    lapack
                                    ${EGADS_LIBRARY} ${EGADSLITE_LIBRARY}
                                    #nlopt
#                                   LP.3
                                    lplib3
                                    Boost::headers # yes this is weird
                                    ${Boost_LIBRARIES}
                                    ${Boost_PROGRAM_OPTIONS_LIBRARY}
                                    ${PETSC_LIBRARIES}
                                    ${MPI_LIBRARIES}
                                    ${BOOST_TRACELIBS}
                                    )

#if(USE_PETSC STREQUAL "True")
#  target_link_libraries(metris PUBLIC PkgConfig::PETSC ${MPI_C_LIBRARIES})
#  #target_link_libraries(metris PUBLIC PETSc)
#endif()

target_link_libraries(egads2lite PUBLIC libMetris 
                                        libmeshb 
                                        ${EGADS_LIBRARY} ${EGADSLITE_LIBRARY}
                                        Boost::headers # yes this is weird
                                        ${Boost_LIBRARIES}
                                        ${Boost_PROGRAM_OPTIONS_LIBRARY})

#For some reason, not specifying the full path does not work for EXISTS
if(EXISTS ${CMAKE_SOURCE_DIR}/progs/scatter_rh.cxx)
  add_executable(scatter_rh "progs/scatter_rh.cxx")
  target_link_libraries(scatter_rh PUBLIC libMetris 
                                          libmeshb 
                                          #libegads libegadslite
                                          #egads egadslite
                                          ${EGADS_LIBRARY} ${EGADSLITE_LIBRARY}
                                          Boost::headers # yes this is weird
                                          ${Boost_LIBRARIES}
                                          ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_include_directories(scatter_rh PRIVATE  "."
                             ${ESP_ROOT}/include)
  setMetrisFlags(scatter_rh PUBLIC)
endif()

# Note on includes: the src lib target (libMetris) already includes publicly everything that we need
# The main.cxx itself only needs main.hxx additionally. 
target_include_directories(metris PRIVATE "progs/main.hxx" "." "SANS"
                                  PUBLIC ${Boost_INCLUDE_DIRS}
                                         ${PETSC_INCLUDE}
                                         ${ABSL_INCLUDE_DIRS}
#                                  PUBLIC ${LPLIB_INCLUDE_DIR}
                                  )
setMetrisFlags(metris PUBLIC)
target_include_directories(egads2lite PRIVATE  "."
                           ${ESP_ROOT}/include)
setMetrisFlags(egads2lite PUBLIC)


# icc having trouble with GTest on my system, figure this out later
#if(NOT("icc" IN_LIST SHORT_COMPILER_NAME OR "gcc-13" IN_LIST SHORT_COMPILER_NAME) AND "Debug" IN_LIST CMAKE_BUILD_TYPE)
if(NOT("icc" IN_LIST SHORT_COMPILER_NAME OR "gcc-13" IN_LIST SHORT_COMPILER_NAME)  AND  DO_UNITTESTS)

  #include(FetchContent)
  #FetchContent_Declare(
  #  googletest
  #  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  #)
  #FetchContent_MakeAvailable(googletest)

  ## --- Testing
  enable_testing()
  #
  ## -- Google tests
  #include(FetchContent)
  #FetchContent_Declare(
  #  googletest
  #  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  #)
  #FetchContent_MakeAvailable(googletest)
  #
  ## - Generate a different test target from each cxx file found in gunit
  #include(GoogleTest)
  #foreach (t ${TEST_SOURCES1})
  #  get_filename_component(d ${t} NAME_WE) # Get name e.g. unit_dbg1 out of $PATHTOFILE/unit_dbg1.cxx
  #  add_executable(${d} ${t})
  #  target_link_libraries(${d} PUBLIC GTest::gtest_main)
  #  target_include_directories(${d} PRIVATE ${googletest_SOURCE_DIR}/googletest/include)
  #  gtest_discover_tests(${d})
  #endforeach()
  
  file (GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "bunit/*.cxx")
  
  # -- Boost/test
  find_package(Boost COMPONENTS unit_test_framework REQUIRED)
  #find_package(NLopt REQUIRED) 
  #if(NLopt_FOUND)
  #  message("NLOP in ${NLopt_INCLUDE_DIRS} ${NLopt_INCLUDE_DIRS}")
  #else()
  #  message(FATAL_ERROR "NLOPT NOT FOUND")
  #endif()
  
  foreach (t ${TEST_SOURCES})
    get_filename_component(d ${t} NAME_WE) # Get name e.g. unit_dbg1 out of $PATHTOFILE/unit_dbg1.cxx
    add_executable(${d} ${t})
    setMetrisFlags(${d} PUBLIC)
    target_link_libraries(${d} PUBLIC libMetris)
    target_link_libraries(${d} PUBLIC ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_link_libraries(${d} PUBLIC ${Boost_PROGRAM_OPTIONS_LIBRARY}
                                      ${BOOST_TRACELIBS}
                                      ${CLP_LIBRARIES}
                                      ${PETSC_LIBRARIES}
                                      ${ABSL_LIBRARIES}
                                      lapack 
                                      #libegads libegadslite 
                                      #egads egadslite
                                      ${EGADS_LIBRARY} ${EGADSLITE_LIBRARY} 
                                      nlopt
                                      libmeshb
                                      pthread)
    if(CLP_FOUND)
      target_link_libraries(${d} PUBLIC ${CLP_LIBRARIES})
    endif()
                          
    target_include_directories(${d} PUBLIC "."  ${PETSC_INCLUDE}  ${ABSL_INCLUDE_DIRS})
  endforeach()
else()
  message("Skipping unit test targets (use clang)")
endif()
    

install(TARGETS metris DESTINATION ${CMAKE_BINARY_DIR})

