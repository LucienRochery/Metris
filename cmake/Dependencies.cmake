include(FetchContent)

include(cmake/MetrisFlags.cmake)
include(cmake/MetrisDependencyHelpers.cmake)


# Unified dependency import pattern:
# metris_find_or_fetch_dependency(NAME <name> [COMPONENTS ...] [GIT_REPOSITORY ...] [ENV_VAR ...] [REQUIRED])

# GMP (env or find)
if(USE_GMP)
  metris_find_or_fetch_dependency(NAME GMP ENV_VAR GMP_DIR REQUIRED)
endif()

# Tracy (FetchContent)
if(USE_TRACY)
  metris_find_or_fetch_dependency(NAME Tracy GIT_REPOSITORY https://github.com/wolfpld/tracy.git GIT_TAG master REQUIRED)
endif()

# GINAC/CLN (env)
if(REQ_CODEGEN)
  metris_find_or_fetch_dependency(NAME GINAC ENV_VAR GINAC_DIR REQUIRED)
  metris_find_or_fetch_dependency(NAME CLN ENV_VAR CLN_DIR REQUIRED)
endif()

# Eigen3


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


# --- Eigen3 Dependency ---
# Improved handling: checks for target from parent project, then environment variable,
# then find_package, then FetchContent
metris_find_eigen3()
# --- fmt Dependency ---
# Ensure fmt is available (find or fetch)
metris_find_fmt()
# --- nlohmann/json Dependency ---
# Ensure nlohmann/json is available (find or fetch)
metris_find_json()
# Find Boost program_options so Boost::program_options target exists
find_package(Boost COMPONENTS program_options QUIET)
if(Boost_FOUND)
  message(STATUS "Found Boost program_options")
  target_link_libraries(metris_deps INTERFACE Boost::program_options)
  metris_register_dependency("find_package" "Boost REQUIRED" "program_options")
endif()
# --- End Eigen3 Dependency ---


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
      #list(APPEND METRIS_DEPS_LIBRARIES    ${MPI_LIBRARIES})
      target_link_libraries(metris_deps INTERFACE ${MPI_LIBRARIES})


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
    #list(APPEND METRIS_DEPS_LIBRARIES    ${PETSC_LIBRARIES})
    target_link_libraries(metris_deps INTERFACE ${PETSC_LIBRARIES})
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
# Robustly check for EGADS header
if(NOT EXISTS "${ESP_ROOT}/include/egads.h")
  message(FATAL_ERROR "EGADS header not found: ${ESP_ROOT}/include/egads.h. Set ESP_ROOT to the directory containing include/egads.h.")
endif()

# Propagate EGADS include path to metris_deps
target_include_directories(metris_deps INTERFACE $<BUILD_INTERFACE:${ESP_ROOT}/include> $<INSTALL_INTERFACE:include/egads>)

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




# --- NLopt Dependency ---
# Improved handling for shared dependencies: checks if NLopt::nlopt target already exists
# from parent project, then tries various methods to find or fetch NLopt.
metris_find_nlopt()
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


