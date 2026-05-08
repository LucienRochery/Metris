include(FetchContent)

include(cmake/MetrisFlags.cmake)

# Since some dependencies can use either find_package or FetchContent, we're going to track how they were imported.
# This is ultimately so we can export to the MetrisConfig.cmake file for find_package(Metris) to work correctly.
set(METRIS_CONFIG_DEPENDENCIES "")
# Note: do not provide DEP_COMPONENTS as ;-separated list as this messes up iterating over METRIS_CONFIG_DEPENDENCIES.
# Types can be "find_package", "FetchContent", or "provided".
# If provided, add first the library, then the include dirs.
function(metris_register_dependency IMPORT_TYPE DEP_NAME DEP_COMPONENTS)
  list(APPEND METRIS_CONFIG_DEPENDENCIES "${IMPORT_TYPE}|${DEP_NAME}|${DEP_COMPONENTS}")
  set(METRIS_CONFIG_DEPENDENCIES ${METRIS_CONFIG_DEPENDENCIES} PARENT_SCOPE)
endfunction()

if(USE_GMP)
  if(NOT DEFINED GMP_DIR)
    if(DEFINED ENV{GMP_DIR})
      set(GMP_DIR $ENV{GMP_DIR})
    else()
      message(FATAL_ERROR "Specify GMP_DIR as CMake option or env var")
    endif()
  endif()
  set(USE_MULTIPRECISION ON)
  add_compile_definitions(USE_GMP)
  add_compile_definitions(USE_MULTIPRECISION)
  set(GMP_INCLUDE_DIRS ${GMP_DIR}/include)
  find_library(GMP_LIBRARIES NAMES gmp
               HINTS $ENV{GMP_DIR}
               PATH_SUFFIXES lib)

  target_link_libraries(metris_deps INTERFACE ${GMP_LIBRARIES})
  #list(APPEND METRIS_DEPS_LIBRARIES ${GMP_LIBRARIES})
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${GMP_INCLUDE_DIRS})
  message("-- GMP_INCLUDE_DIRS = ${GMP_INCLUDE_DIRS}")
  message("-- GMP_LIBRARIES    = ${GMP_LIBRARIES}")
endif()


# Profiling tool
if(USE_TRACY)
  FetchContent_Declare(
    Tracy_fetch
    GIT_REPOSITORY https://github.com/wolfpld/tracy.git
    GIT_TAG master
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
  )
  set(TRACY_ENABLE ON)
  FetchContent_MakeAvailable(Tracy_fetch)
  set(TRACY_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/tracy_fetch-src/public/tracy/")
  #message("TRACY_INCLUDE_DIRS = ${TRACY_INCLUDE_DIRS}")
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${TRACY_INCLUDE_DIRS})
  #list(APPEND METRIS_DEPS_LIBRARIES Tracy::TracyClient)
  target_link_libraries(metris_deps INTERFACE Tracy::TracyClient)
  add_compile_definitions(TRACY_ENABLE)
endif()

## Used for mlpack::BallTree
#FetchContent_Declare(
#  mlpack_fetch
#  GIT_REPOSITORY https://github.com/mlpack/mlpack.git
#  GIT_TAG master
#  GIT_SHALLOW TRUE
#  GIT_PROGRESS TRUE
#)
#set(BUILD_CLI_EXECUTABLES OFF CACHE BOOL "" FORCE)
#set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
#set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
#set(BUILD_JULIA_BINDINGS OFF CACHE BOOL "" FORCE)
#set(BUILD_GO_BINDINGS OFF CACHE BOOL "" FORCE)
#set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
#set(BUILD_MARKDOWN_BINDINGS OFF CACHE BOOL "" FORCE)
#set(ARMA_METRIS_USE_LAPACK OFF CACHE BOOL "" FORCE)
#set(ARMA_USE_BLAS OFF CACHE BOOL "" FORCE)
#FetchContent_MakeAvailable(mlpack_fetch)
#set(MLPACK_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/mlpack_fetch-src/src/")
#list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${MLPACK_INCLUDE_DIRS})
#list(APPEND METRIS_DEPS_LIBRARIES Tracy::TracyClient)

if(REQ_CODEGEN)
  if(NOT DEFINED ENV{GINAC_DIR} AND NOT (DEFINED GINAC_LIBRARIES AND DEFINED GINAC_INCLUDE_DIRS))
    message(FATAL_ERROR "Set environment variable GINAC_DIR to directory containing lib/libginac and include/ginac/ginac.h")
  endif()
  if(NOT DEFINED ENV{CLN_DIR} AND NOT DEFINED CLN_INCLUDE_DIRS)
    message(FATAL_ERROR "Set environment variable GINAC_DIR to directory containing lib/libginac and include/ginac/ginac.h")
  endif()

  set(GINAC_DIR $ENV{GINAC_DIR})
  if(NOT DEFINED GINAC_LIBRARIES)
    find_library(GINAC_LIBRARIES NAMES ginac
                 HINTS $ENV{GINAC_DIR}
                 PATH_SUFFIXES lib)
  endif()
  if(NOT DEFINED GINAC_INCLUDE_DIRS)
    set(GINAC_INCLUDE_DIRS "$ENV{GINAC_DIR}/include")
  endif()
  if(NOT DEFINED CLN_INCLUDE_DIRS)
    set(CLN_INCLUDE_DIRS "$ENV{CLN_DIR}/include")
  endif()

  message("-- GINAC_INCLUDE_DIRS = ${GINAC_INCLUDE_DIRS}")
  message("-- CLN_INCLUDE_DIRS   = ${CLN_INCLUDE_DIRS}")
  message("-- GINAC_LIBRARIES    = ${GINAC_LIBRARIES}")

  #list(APPEND METRIS_DEPS_LIBRARIES ${GINAC_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${GINAC_LIBRARIES})
  #list(APPEND METRIS_DEPS_LIBRARIES ${GINAC_LIBRARIES})
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${GINAC_INCLUDE_DIRS} ${CLN_INCLUDE_DIRS})
endif()


set(EIGEN3_INCLUDE_DIRS "$ENV{EIGEN_DIR}")
if (NOT EIGEN3_INCLUDE_DIRS)
  FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG bcce88c99ed687b756b7a537554cb7c1780b816e
    #GIT_TAG master
    #FIND_PACKAGE_ARGS NAMES Eigen3
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(Eigen3)
  set(EIGEN3_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/eigen3-src/")
  if(METRIS_INSTALL)
    install(DIRECTORY ${EIGEN3_INCLUDE_DIRS}/Eigen ${EIGEN3_INCLUDE_DIRS}/unsupported
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  message("EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}")
  #list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIRS})
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>)
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
  
endif()


if(METRIS_USE_LAPACK)
  message("METRIS_USE_LAPACK was set to ON")
  
  add_compile_definitions(METRIS_USE_LAPACK)

  enable_language(Fortran)
  include(cmake/FindLAPACK.cmake)
  message("-- LAPACK_LIBRARIES    = ${LAPACK_LIBRARIES}")
  #list(APPEND METRIS_DEPS_LIBRARIES ${LAPACK_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${LAPACK_LIBRARIES})

  find_package(BLAS REQUIRED)
  message("-- BLAS_LIBRARIES    = ${BLAS_LIBRARIES}")
  #list(APPEND METRIS_DEPS_LIBRARIES ${BLAS_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${BLAS_LIBRARIES})

  find_library(GFORTRAN_LIBRARIES NAMES gfortran
               HINTS /usr/local/lib /opt/homebrew/opt/gcc/lib/gcc/current/
               )
  if(NOT GFORTRAN_LIBRARIES)
    message("-- GFORTRAN_LIBRARIES not found, defaulting to -lgfortran")
   # set(GFORTRAN_LIBRARIES gfortran)
   set(GFORTRAN_LIBRARIES "")
  else()
    message("-- GFORTRAN_LIBRARIES = ${GFORTRAN_LIBRARIES}")
  endif()
  #list(APPEND METRIS_DEPS_LIBRARIES ${GFORTRAN_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${GFORTRAN_LIBRARIES})
endif()




if(NOT DEFINED ESP_ROOT)
  if(DEFINED ENV{ESP_ROOT})
    set(ESP_ROOT $ENV{ESP_ROOT})
  elseif(DEFINED ENV{ESP_DIR})
    set(ESP_ROOT $ENV{ESP_DIR})
  else()
    message(FATAL_ERROR "Set environment variable ESP_ROOT (here $ENV{ESP_ROOT})"
      " or ESP_DIR (here $ENV{ESP_DIR}) to the folder containing include/egads.h")
  endif()
endif()

if(NOT DEFINED CAS_ROOT)
  if(DEFINED ENV{CASROOT})
    set(CAS_ROOT $ENV{CASROOT})
  elseif(DEFINED ENV{CAS_ROOT})
    set(CAS_ROOT $ENV{CAS_ROOT})
  elseif(DEFINED ENV{CAS_DIR})
    set(CAS_ROOT $ENV{CAS_DIR})
  else()
    message(FATAL_ERROR "Set environment variable CASROOT (here $ENV{CASROOT})"
      " or CAS_ROOT (here $ENV{CAS_ROOT}) or CAS_DIR (here $ENV{CAS_DIR})"
      " to the folder containing bin/ include/ lib/ share/ of OpenCascade library to link to.")
  endif()
endif()

#if(NOT DEFINED CAS_REV)
#  if(DEFINED ENV{CASREV})
#    set(CAS_REV $ENV{CASREV})
#  elseif(DEFINED ENV{CAS_REV})
#    set(CAS_REV $ENV{CAS_REV})
#  elseif(DEFINED ENV{CAS_VERSION})
#    set(CAS_REV $ENV{CAS_VERSION})
#  else()
#    message(FATAL_ERROR "Set environment variable CASREV/CAS_REV/CAS_VERSION, usually by sourcing ESP_env.sh")
#  endif()
#endif()
#
#message("Using ESP_ROOT = ${ESP_ROOT}")
#message("Using CAS_ROOT = ${CAS_ROOT}, CAS_REV = ${CAS_REV}")

# Linker still needs path to lib files at compile time
find_library(EGADS_LIBRARIES NAMES egads
             HINTS ${ESP_ROOT}
             PATH_SUFFIXES lib)
find_library(EGADSLITE_LIBRARIES NAMES egadslite
             HINTS ${ESP_ROOT}
             PATH_SUFFIXES lib)
list(APPEND CMAKE_BUILD_RPATH   ${ESP_ROOT}/lib/)
list(APPEND CMAKE_INSTALL_RPATH ${ESP_ROOT}/lib/)
list(APPEND CMAKE_BUILD_RPATH   ${CAS_ROOT}/lib/)
list(APPEND CMAKE_INSTALL_RPATH ${CAS_ROOT}/lib/)

set(EGADS_INCLUDE_DIRS ${ESP_ROOT}/include)

#list(APPEND METRIS_DEPS_LIBRARIES    ${EGADS_LIBRARIES})
#list(APPEND METRIS_DEPS_LIBRARIES    ${EGADSLITE_LIBRARIES})
target_link_libraries(metris_deps INTERFACE ${EGADS_LIBRARIES} 
                                             ${EGADSLITE_LIBRARIES})
list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${EGADS_INCLUDE_DIRS})


#if(${CAS_REV} VERSION_GREATER_EQUAL 7.8) 
#  set(OCC_LIBRARY_NAMES TKBO TKBRep TKBin TKBinL TKBinTObj TKBinXCAF TKBool TKCAF TKCDF TKDE TKDECascade TKDEGLTF TKDEIGES TKDEOBJ TKDEPLY TKDESTEP TKDESTL TKDEVRML TKExpress TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKLCAF TKMath TKMesh TKOffset TKPrim TKRWMesh TKService TKShHealing TKStd TKStdL TKTObj TKTopAlgo TKV3d TKVCAF TKXCAF TKXMesh TKXSBase TKXml TKXmlL TKXmlTObj TKXmlXCAF TKernel)
#  message(STATUS "Using OCC version >= 7.8 library names: ${OCC_LIBRARY_NAMES}")
#else()
#  set(OCC_LIBRARY_NAMES TKBool TKernel TKFeat TKBO TKGeomAlgo TKMath TKOffset TKPrim TKTopAlgo TKBRep TKG2d TKG3d TKGeomBase TKShHealing TKSTEP TKSTEP209 TKSTEPBase TKSTEPAttr TKXSBase TKIGES TKFillet)
#  message(STATUS "Using OCC version <= 7.7 library names: ${OCC_LIBRARY_NAMES}")
#endif()

#foreach(OCClib ${OCC_LIBRARY_NAMES})
#  find_library(CAS_${OCClib} ${OCClib}
#               HINTS ${CAS_ROOT}/lib
#               REQUIRED)
#  list(APPEND CAS_LIBRARIES ${CAS_${OCClib}})
#endforeach()
#
#list(APPEND METRIS_DEPS_LIBRARIES ${CAS_LIBRARIES})



if(USE_CLP)
  message("Using CLP")
  include(FindCLP)
  if(NOT(CLP_FOUND))
    message(WARNING "CLP was not found on this system.")
  else()
    add_compile_definitions(USE_CLP)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${CLP_INCLUDE_DIRS})
    #list(APPEND METRIS_DEPS_LIBRARIES    ${CLP_LIBRARIES})
    target_link_libraries(metris_deps INTERFACE ${CLP_LIBRARIES})
  endif()
endif()

# External libraries to be fetched
find_package(fmt QUIET)
#message(WARNING "Disabled find fmt")
if(fmt_FOUND)
  message(STATUS "fmt lib found: ${fmt_VERSION}")
  metris_register_dependency("find_package" "fmt REQUIRED" "")
else()
  message(STATUS "fmt lib not found, fetching")
  # fmt library: std::format precursor with better performance
  # and pre-C++20 support.
  FetchContent_Declare(
    fmt_fetch
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281
    EXCLUDE_FROM_ALL) # 10.2.1
  FetchContent_MakeAvailable(fmt_fetch)

  FetchContent_GetProperties(fmt_fetch)
  if(NOT fmt_fetch_POPULATED)
    message(FATAL_ERROR "fmt was not fetched correctly.")
  endif()

  if(METRIS_INSTALL)
    install(DIRECTORY ${fmt_fetch_SOURCE_DIR}/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h*")
  endif()

  # This is to suppress the warning that PUBLIC_HEADER exists but not
  # being installed. For some reason, we can't get that to install 
  # correctly, and are installing headers manually. 
  set_target_properties(fmt PROPERTIES PUBLIC_HEADER "")

  if(METRIS_INSTALL)
    install(TARGETS fmt
            EXPORT libMetrisTargets
            LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
            #PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            )
  endif()

  if(NOT TARGET fmt::fmt)
    add_library(fmt::fmt ALIAS fmt)
  endif()

  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${fmt_fetch_SOURCE_DIR}/include>)
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
  metris_register_dependency("FetchContent" "fmt" "")
endif()
#list(APPEND METRIS_DEPS_LIBRARIES fmt::fmt)
target_link_libraries(metris_deps INTERFACE fmt::fmt)

include(FetchContent)

# nlohmann json library
# Skip our own fetch if a parent project (or earlier find_package) has already
# defined the nlohmann_json target; otherwise the second add_subdirectory of
# json's tree collides with the first by re-creating the same target.
if(TARGET nlohmann_json)
  # Supplied by a parent project (e.g. Metris built via add_subdirectory by a
  # consumer that already pulls nlohmann_json). Reuse their target and pull
  # its include path into METRIS_EXTERNAL_INCLUDE_DIRS so libMetris sources
  # (which include "nlohmann/json_fwd.hpp") see the header.
  if(NOT TARGET nlohmann_json::nlohmann_json)
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
  endif()
  get_target_property(_metris_nj_includes nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
  if(_metris_nj_includes)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${_metris_nj_includes})
  endif()
  unset(_metris_nj_includes)
  metris_register_dependency("provided" "json" "")
else()
  FetchContent_Declare(json_fetch URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
  FetchContent_MakeAvailable(json_fetch)
  if(NOT TARGET nlohmann_json::nlohmann_json)
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
  endif()
  if(METRIS_INSTALL)
    install(FILES ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json.hpp
                  ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json_fwd.hpp
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann/
            )
  endif()
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${json_fetch_SOURCE_DIR}/single_include>)
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
  metris_register_dependency("FetchContent" "json" "")
endif()
#target_link_libraries(foo PRIVATE nlohmann_json::nlohmann_json)

if(USE_ABSL)
  message("Enabled absl")
  FetchContent_Declare(
    fetch_absl
    #URL https://github.com/abseil/abseil-cpp/archive/e7fe9ec9ebfc6607765d489b76c9954e0a88c5d4.zip
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
    #GIT_TAG e7fe9ec9ebfc6607765d489b76c9954e0a88c5d4
    #FIND_PACKAGE_ARGS NAMES absl
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(fetch_absl)
  set(ABSL_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/fetch_absl-src/")
  set(ABSL_LIBRARIES absl::hash absl::flat_hash_map)
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${ABSL_INCLUDE_DIRS})
  #list(APPEND METRIS_DEPS_LIBRARIES    ${ABSL_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${ABSL_LIBRARIES})
  message("-- ABSL_INCLUDE_DIRS = ${ABSL_INCLUDE_DIRS}")
  message("-- ABSL_LIBRARIES = ${ABSL_LIBRARIES}")
endif()


##if(USE_MULTIPRECISION)
###  list(APPEND METRIS_BOOST_COMPONENTS headers REQUIRED)
##  if(NOT(Boost_headers_FOUND))
##    message(FATAL_ERROR "Boost headers not found: either clone all Boost libs here or install system wide")
##  endif()
##endif()

#find_package(Boost COMPONENTS program_options)
#if(NOT(Boost_program_options_FOUND))
#  FetchContent_Declare(
#    fetch_program_options
#    GIT_REPOSITORY https://github.com/boostorg/program_options.git
#    GIT_TAG master
#    #GIT_SHALLOW TRUE
#    #FIND_PACKAGE_ARGS NAMES Boost COMPONENTS program_options REQUIRED
#    EXCLUDE_FROM_ALL
#  )
#  FetchContent_MakeAvailable(fetch_program_options)
#  list(APPEND METRIS_DEPS_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
#  metris_register_dependency("FetchContent" "Boost" "program_options")
#else()
#  list(APPEND METRIS_DEPS_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
#  metris_register_dependency("find_package" "Boost" "program_options")
#endif()
#set(METRIS_BOOST_COMPONENTS "${METRIS_BOOST_COMPONENTS}   program_options") #math # exception
find_package(Boost REQUIRED COMPONENTS program_options) #math exception
#list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
#message("-- Boost_PROGRAM_OPTIONS_LIBRARY = ${Boost_PROGRAM_OPTIONS_LIBRARY}")
#message("-- Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
metris_register_dependency("find_package" "Boost REQUIRED" "program_options") #math exception
#list(APPEND METRIS_DEPS_LIBRARIES Boost::program_options)
target_link_libraries(metris_deps INTERFACE Boost::program_options)


# --- NLopt Dependency ---

# A CMake project including Metris should prefer to find its own NLopt installation and provide
# NLOPT_LIBRARIES and NLOPT_INCLUDE_DIRS.
if(DEFINED NLOPT_LIBRARIES AND DEFINED NLOPT_INCLUDE_DIRS)
  message(STATUS "Using provided NLOPT_LIBRARIES and NLOPT_INCLUDE_DIRS.")
  message(STATUS "NLOPT_LIBRARIES = ${NLOPT_LIBRARIES}")
  message(STATUS "NLOPT_INCLUDE_DIRS = ${NLOPT_INCLUDE_DIRS}")
  #list(APPEND METRIS_DEPS_LIBRARIES ${NLOPT_LIBRARIES})
  target_link_libraries(metris_deps INTERFACE ${NLOPT_LIBRARIES})
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS})
  metris_register_dependency("provided" "NLopt" "${NLOPT_LIBRARIES} ${NLOPT_INCLUDE_DIRS}")
elseif(DEFINED NLOPT_DIR OR DEFINED ENV{NLOPT_DIR})
  # Otherwise, an NLOPT_DIR can be passed in which should include the relevant CMake configuration
  # files (nloptConfig.cmake or nlopt-config.cmake).
  if(NOT DEFINED NLOPT_DIR)
    set(NLOPT_DIR $ENV{NLOPT_DIR})
    message(STATUS "Using NLOPT_DIR from environment: ${NLOPT_DIR}")
  endif()
  message(STATUS "Attempting to find external NLopt using hint NLOPT_DIR=${NLOPT_DIR}...")
  # find_package is the most robust way to do this. It will use NLOPT_DIR as a hint.
  find_package(NLopt REQUIRED HINTS "${NLOPT_DIR}")
  
  message(STATUS "Found external NLopt: ${NLOPT_LIBRARIES}")
  #list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)
  target_link_libraries(metris_deps INTERFACE NLopt::nlopt)

  metris_register_dependency("find_package" "NLopt REQUIRED HINTS ${NLOPT_DIR}" "")
else()
  # First try finding the package. If that fails, then clone and build.

  find_package(NLopt)

  if(NLopt_FOUND)
    message(STATUS "Found NLopt libraries: ${NLOPT_LIBRARIES}")
    message(STATUS "Found NLopt include directories: ${NLOPT_INCLUDE_DIRS}")
    #list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)
    target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
    metris_register_dependency("find_package" "NLopt REQUIRED" "")
  else()
    # Lastly, fetch and build our own for standalone builds.
    message(STATUS "find_package(NLopt) failed, cloning NLopt.")
    FetchContent_Declare(
        nlopt_fetch
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        019f61ac7253a537760d9cdd9febd927ec97320c
    )
    
    set(NLOPT_BUILD_TESTS OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(nlopt_fetch)
    
    # After FetchContent, an 'nlopt' target is available.
    # The generated config header (nlopt.h) is in the binary directory of the fetch.
    FetchContent_GetProperties(nlopt_fetch)
    if(NOT nlopt_fetch_POPULATED)
      message(FATAL_ERROR "NLopt was not fetched correctly.")
    endif()
    # Add NLopt library directory to RPATH so the runtime linker can find libnlopt.so
    list(APPEND CMAKE_BUILD_RPATH   ${nlopt_fetch_BINARY_DIR})
    list(APPEND CMAKE_INSTALL_RPATH ${nlopt_fetch_BINARY_DIR})
    # We do this because find_package() uses this naming.
    if(NOT TARGET NLopt::nlopt)
      add_library(NLopt::nlopt ALIAS nlopt)
    endif()
    if(METRIS_INSTALL)
      install(DIRECTORY ${nlopt_fetch_BINARY_DIR}
              DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
              FILES_MATCHING PATTERN "*.h*")
      install(TARGETS nlopt
              EXPORT libMetrisTargets
              LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
              ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
              RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
            # Install NLopt headers if using FetchContent
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${nlopt_fetch_BINARY_DIR}>)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
    # CMAKE_BUILD_WITH_INSTALL_RPATH is TRUE, so build RPATHs use CMAKE_INSTALL_RPATH.
    # Append the FetchContent build dir so executables can find libnlopt.so at runtime.
    list(APPEND CMAKE_INSTALL_RPATH "${nlopt_fetch_BINARY_DIR}")
    #list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)
    target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
    metris_register_dependency("FetchContent" "NLopt" "")
  endif()
endif()
# --- End NLopt Dependency ---


include(FetchContent)


FetchContent_MakeAvailable(${FETCH_LIST})



# This is necessary to make the sanitizer work correctly. Also we want to
# propagate flags.
if(USE_ABSL)
  setMetrisCXXFlags(absl_hash INTERFACE)
  setMetrisCXXFlags(absl_flat_hash_map INTERFACE)
  setMetrisCXXFlags(absl_spinlock_wait INTERFACE)
  setMetrisCXXFlags(absl_int128 INTERFACE)
  setMetrisCXXFlags(absl_exponential_biased INTERFACE)
  setMetrisCXXFlags(absl_log_severity INTERFACE)
  setMetrisCXXFlags(absl_civil_time INTERFACE)
  setMetrisCXXFlags(absl_raw_logging_internal INTERFACE)
  setMetrisCXXFlags(absl_time_zone INTERFACE)
  setMetrisCXXFlags(absl_bad_variant_access INTERFACE)
  setMetrisCXXFlags(absl_debugging_internal INTERFACE)
  setMetrisCXXFlags(absl_cordz_functions INTERFACE)
  setMetrisCXXFlags(absl_bad_optional_access INTERFACE)
  setMetrisCXXFlags(absl_throw_delegate INTERFACE)
  setMetrisCXXFlags(absl_base INTERFACE)
  setMetrisCXXFlags(absl_stacktrace INTERFACE)
  setMetrisCXXFlags(absl_crc_cpu_detect INTERFACE)
  setMetrisCXXFlags(absl_demangle_internal INTERFACE)
  setMetrisCXXFlags(absl_string_view INTERFACE)
  setMetrisCXXFlags(absl_city INTERFACE)
  setMetrisCXXFlags(absl_malloc_internal INTERFACE)
  setMetrisCXXFlags(absl_low_level_hash INTERFACE)
  setMetrisCXXFlags(absl_strings_internal INTERFACE)
  setMetrisCXXFlags(absl_crc_internal INTERFACE)
  setMetrisCXXFlags(absl_graphcycles_internal INTERFACE)
  setMetrisCXXFlags(absl_strings INTERFACE)
  setMetrisCXXFlags(absl_hash INTERFACE)
  setMetrisCXXFlags(absl_symbolize INTERFACE)
  setMetrisCXXFlags(absl_time INTERFACE)
  setMetrisCXXFlags(absl_str_format_internal INTERFACE)
  setMetrisCXXFlags(absl_kernel_timeout_internal INTERFACE)
  setMetrisCXXFlags(absl_crc32c INTERFACE)
  setMetrisCXXFlags(absl_crc_cord_state INTERFACE)
  setMetrisCXXFlags(absl_synchronization INTERFACE)
  setMetrisCXXFlags(absl_cord_internal INTERFACE)
  setMetrisCXXFlags(absl_cordz_handle INTERFACE)
  setMetrisCXXFlags(absl_hashtablez_sampler INTERFACE)
  setMetrisCXXFlags(absl_cordz_info INTERFACE)
  setMetrisCXXFlags(absl_raw_hash_set INTERFACE)
  setMetrisCXXFlags(absl_cord INTERFACE)
endif()


