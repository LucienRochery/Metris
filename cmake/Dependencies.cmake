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

  list(APPEND METRIS_DEPS_LIBRARIES ${GMP_LIBRARIES})
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
  list(APPEND METRIS_DEPS_LIBRARIES Tracy::TracyClient)
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

  list(APPEND METRIS_DEPS_LIBRARIES ${GINAC_LIBRARIES})
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
  install(DIRECTORY ${EIGEN3_INCLUDE_DIRS}/Eigen ${EIGEN3_INCLUDE_DIRS}/unsupported
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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
  list(APPEND METRIS_DEPS_LIBRARIES ${LAPACK_LIBRARIES})

  find_package(BLAS REQUIRED)
  message("-- BLAS_LIBRARIES    = ${BLAS_LIBRARIES}")
  list(APPEND METRIS_DEPS_LIBRARIES ${BLAS_LIBRARIES})

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
  list(APPEND METRIS_DEPS_LIBRARIES ${GFORTRAN_LIBRARIES})
endif()


if(METRIS_USE_PETSC)
  message("PETSC enabled")

  if(PETSC_FOUND AND NOT PETSC_PKGCONFIG)

    message(STATUS "PETSC_FOUND = ${PETSC_FOUND}")
    message(STATUS "PETSC_LIBRARIES = ${PETSC_LIBRARIES}")

  else()

    if(DEFINED PETSC_DIR AND DEFINED PETSC_ARCH)

      message(STATUS "PETSC_DIR = ${PETSC_DIR} and PETSC_ARCH = ${PETSC_ARCH} already set: assuming PETSc target already defined")
      #find_package(MPI)
      set(PETSC ${PETSC_DIR}/${PETSC_ARCH})
      set(PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)
      list(APPEND PETSC_INCLUDE_DIRS ${PETSC_DIR}/${PETSC_ARCH}/include)
      list(APPEND PETSC_INCLUDE_DIRS ${MPI_INCLUDE_PATH})

      message("PETSC_LIBRARIES = ${PETSC_LIBRARIES}")

    else()

      if(NOT(DEFINED ENV{PETSC_DIR}))
        message(FATAL_ERROR "Set environment variables PETSC_DIR")
      endif()
      if(NOT(DEFINED ENV{PETSC_ARCH}))
        message(WARNING "PETSC_ARCH not set.")
      endif()
      set(PETSC_DIR  $ENV{PETSC_DIR})
      set(PETSC_ARCH $ENV{PETSC_ARCH})
      find_package(MPI)
      message("-- MPI_INCLUDE_DIRS = ${MPI_INCLUDE_DIRS}")
      message("-- MPI_LIBRARIES    = ${MPI_LIBRARIES}")
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${MPI_INCLUDE_DIRS})
      list(APPEND METRIS_DEPS_LIBRARIES    ${MPI_LIBRARIES})


      set(PETSC_ARCH_DIR ${PETSC_DIR}/${PETSC_ARCH})
      set(PETSC_INCLUDE_DIRS ${PETSC_DIR}/include)

      list(APPEND PETSC_INCLUDE_DIRS ${PETSC_DIR}/${PETSC_ARCH}/include)
      list(APPEND PETSC_INCLUDE_DIRS ${MPI_INCLUDE_PATH})

      set(ENV{PKG_CONFIG_PATH} ${PETSC_ARCH_DIR}/lib/pkgconfig)
      pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
      if(NOT PETSC_FOUND)
        message(FATAL_ERROR "PETSC NOT FOUND !")
      endif()

      message(" CASE 2 PETSC_LIBRARIES = ${PETSC_LIBRARIES}")
      set(PETSC_LIBRARIES PkgConfig::PETSC ${MPI_C_LIBRARIES} CACHE INTERNAL "PETSC LIBRARIES (internal)")
      message(" Set to PkgConfig + MPI_C_LIBRARIES -> PETSC_LIBRARIES = ${PETSC_LIBRARIES}")
      set(PETSC_PKGCONFIG ON CACHE INTERNAL "Call PkgConfig::PETSC each reconfig")

    endif()

    add_compile_definitions(METRIS_USE_PETSC)

    message("-- PETSC_INCLUDE_DIRS = ${PETSC_INCLUDE_DIRS}")
    message("-- PETSC_LIBRARIES    = ${PETSC_LIBRARIES}")

    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS})
    list(APPEND METRIS_DEPS_LIBRARIES    ${PETSC_LIBRARIES})
  endif()
else()
  message("PETSC disabled")
  set(PETSC_INCLUDE_DIRS "")
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

message("Using ESP_ROOT = ${ESP_ROOT}")
message("Using CAS_ROOT = ${CAS_ROOT}")

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

list(APPEND METRIS_DEPS_LIBRARIES    ${EGADS_LIBRARIES})
list(APPEND METRIS_DEPS_LIBRARIES    ${EGADSLITE_LIBRARIES})
list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${EGADS_INCLUDE_DIRS})

find_library(CAS_TKBool       TKBool
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKernel      TKernel
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKFeat       TKFeat
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKBO         TKBO
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKGeomAlgo   TKGeomAlgo
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKMath       TKMath
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKOffset     TKOffset
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKPrim       TKPrim
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKTopAlgo    TKTopAlgo
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKBRep       TKBRep
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKG2d        TKG2d
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKG3d        TKG3d
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKGeomBase   TKGeomBase
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKShHealing  TKShHealing
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKSTEP       TKSTEP
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKSTEP209    TKSTEP209
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKSTEPBase   TKSTEPBase
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKSTEPAttr   TKSTEPAttr
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKXSBase     TKXSBase
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKIGES       TKIGES
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
find_library(CAS_TKFillet     TKFillet
             HINTS ${CAS_ROOT}/lib
             REQUIRED)
list(APPEND CAS_LIBRARIES 
            ${CAS_TKBool}
            ${CAS_TKernel}
            ${CAS_TKFeat}
            ${CAS_TKBO}
            ${CAS_TKGeomAlgo}
            ${CAS_TKMath}
            ${CAS_TKOffset}
            ${CAS_TKPrim}
            ${CAS_TKTopAlgo}
            ${CAS_TKBRep}
            ${CAS_TKG2d}
            ${CAS_TKG3d}
            ${CAS_TKGeomBase}
            ${CAS_TKShHealing}
            ${CAS_TKSTEP}
            ${CAS_TKSTEP209}
            ${CAS_TKSTEPBase}
            ${CAS_TKSTEPAttr}
            ${CAS_TKXSBase}
            ${CAS_TKIGES}
            ${CAS_TKFillet})
list(APPEND METRIS_DEPS_LIBRARIES ${CAS_LIBRARIES})



if(USE_CLP)
  message("Using CLP")
  include(FindCLP)
  if(NOT(CLP_FOUND))
    message(WARNING "CLP was not found on this system.")
  else()
    add_compile_definitions(USE_CLP)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${CLP_INCLUDE_DIRS})
    list(APPEND METRIS_DEPS_LIBRARIES    ${CLP_LIBRARIES})
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

  install(DIRECTORY ${fmt_fetch_SOURCE_DIR}/include/
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING PATTERN "*.h*")

  # This is to suppress the warning that PUBLIC_HEADER exists but not
  # being installed. For some reason, we can't get that to install 
  # correctly, and are installing headers manually. 
  set_target_properties(fmt PROPERTIES PUBLIC_HEADER "")

  install(TARGETS fmt
          EXPORT libMetrisTargets
          LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
          #PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          )

  if(NOT TARGET fmt::fmt)
    add_library(fmt::fmt ALIAS fmt)
  endif()

  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${fmt_fetch_SOURCE_DIR}/include>)
  list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
  metris_register_dependency("FetchContent" "fmt" "")
endif()
list(APPEND METRIS_DEPS_LIBRARIES fmt::fmt)

include(FetchContent)

# nlohmann json library
FetchContent_Declare(json_fetch URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json_fetch)
if(NOT TARGET nlohmann_json::nlohmann_json)
  add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
endif()
install(FILES ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json.hpp
              ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json_fwd.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann/
        )
list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${json_fetch_SOURCE_DIR}/single_include>)
list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
metris_register_dependency("FetchContent" "json" "")
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
  list(APPEND METRIS_DEPS_LIBRARIES    ${ABSL_LIBRARIES})
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
list(APPEND METRIS_DEPS_LIBRARIES Boost::program_options)


# --- NLopt Dependency ---

# A CMake project including Metris should prefer to find its own NLopt installation and provide
# NLOPT_LIBRARIES and NLOPT_INCLUDE_DIRS.
if(DEFINED NLOPT_LIBRARIES AND DEFINED NLOPT_INCLUDE_DIRS)
  message(STATUS "Using provided NLOPT_LIBRARIES and NLOPT_INCLUDE_DIRS.")
  message(STATUS "NLOPT_LIBRARIES = ${NLOPT_LIBRARIES}")
  message(STATUS "NLOPT_INCLUDE_DIRS = ${NLOPT_INCLUDE_DIRS}")
  list(APPEND METRIS_DEPS_LIBRARIES ${NLOPT_LIBRARIES})
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
  list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)

  metris_register_dependency("find_package" "NLopt REQUIRED HINTS ${NLOPT_DIR}" "")
else()
  # First try finding the package. If that fails, then clone and build.

  find_package(NLopt)

  if(NLopt_FOUND)
    message(STATUS "Found NLopt libraries: ${NLOPT_LIBRARIES}")
    message(STATUS "Found NLopt include directories: ${NLOPT_INCLUDE_DIRS}")
    list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)
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
    # We do this because find_package() uses this naming.
    if(NOT TARGET NLopt::nlopt)
      add_library(NLopt::nlopt ALIAS nlopt)
    endif()
    install(DIRECTORY ${nlopt_fetch_BINARY_DIR}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h*")
    install(TARGETS nlopt
            EXPORT libMetrisTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
            # Install NLopt headers if using FetchContent
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${nlopt_fetch_BINARY_DIR}>)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
    list(APPEND METRIS_DEPS_LIBRARIES NLopt::nlopt)
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


