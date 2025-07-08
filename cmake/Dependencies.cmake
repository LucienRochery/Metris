include(FetchContent)

include(cmake/MetrisFlags.cmake)


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
  set(GMP_INCLUDE_DIRS ${GMP_DIR}/include)
  find_library(GMP_LIBRARIES NAMES gmp
               HINTS $ENV{GMP_DIR}
               PATH_SUFFIXES lib)

  list(APPEND METRIS_DEPS_LIBRARIES ${GMP_LIBRARIES})
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${GMP_INCLUDE_DIRS})
  message("-- GMP_INCLUDE_DIRS = ${GMP_INCLUDE_DIRS}")
  message("-- GMP_LIBRARIES    = ${GMP_LIBRARIES}")
endif()

if(USE_MULTIPRECISION)
  find_package(Boost COMPONENTS headers)
  if(NOT(Boost_headers_FOUND))
    message(FATAL_ERROR "Boost headers not found: either clone all Boost libs here or install system wide")
    #message(WARNING "find_package(Boost COMPONENTS multiprecision) failed, cloning")
    #FetchContent_Declare(
    #  fetch_multiprecision
    #  GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
    #  GIT_TAG master
    #  EXCLUDE_FROM_ALL
    #)
    #list(APPEND METRIS_DEPS_LIBRARIES Boost::multiprecision)
    #FetchContent_MakeAvailable(fetch_multiprecision)
  endif()
  #list(APPEND METRIS_DEPS_LIBRARIES ${Boost_MULTIPRECISION_LIBRARY})
  #message("-- Boost_LIBRARIES = ${Boost_LIBRARIES}")
  #message("-- Boost_MULTIPRECISION_LIBRARY = ${Boost_MULTIPRECISION_LIBRARY}")
  #message("-- Boost_MULTIPRECISION_LIBRARIES = ${Boost_MULTIPRECISION_LIBRARIES}")
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
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${TRACY_INCLUDE_DIRS})
  list(APPEND METRIS_DEPS_LIBRARIES Tracy::TracyClient)
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

  list(APPEND METRIS_DEPS_LIBRARIES ${GINAC_LIBRARIES})
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${GINAC_INCLUDE_DIRS} ${CLN_INCLUDE_DIRS})
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
  message("EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}")
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIRS})
endif()

message("EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}")


if(USE_LAPACK)
  add_compile_definitions(USE_LAPACK)

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


if(USE_PETSC)
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
      list(APPEND METRIS_DEPS_INCLUDE_DIRS ${MPI_INCLUDE_DIRS})
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

    add_compile_definitions(USE_PETSC)

    message("-- PETSC_INCLUDE_DIRS = ${PETSC_INCLUDE_DIRS}")
    message("-- PETSC_LIBRARIES    = ${PETSC_LIBRARIES}")

    list(APPEND METRIS_DEPS_INCLUDE_DIRS ${PETSC_INCLUDE_DIRS})
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
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${EGADS_INCLUDE_DIRS})
list(APPEND METRIS_DEPS_LIBRARIES    ${EGADSLITE_LIBRARIES})


if(USE_CLP)
  message("Using CLP")
  include(FindCLP)
  if(NOT(CLP_FOUND))
    message(WARNING "CLP was not found on this system.")
  else()
    add_compile_definitions(USE_CLP)
    list(APPEND METRIS_DEPS_INCLUDE_DIRS ${CLP_INCLUDE_DIRS})
    list(APPEND METRIS_DEPS_LIBRARIES    ${CLP_LIBRARIES})
  endif()
endif()

# External libraries to be fetched


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
  list(APPEND METRIS_DEPS_INCLUDE_DIRS ${ABSL_INCLUDE_DIRS})
  list(APPEND METRIS_DEPS_LIBRARIES    ${ABSL_LIBRARIES})
  message("-- ABSL_INCLUDE_DIRS = ${ABSL_INCLUDE_DIRS}")
  message("-- ABSL_LIBRARIES = ${ABSL_LIBRARIES}")
endif()


find_package(Boost COMPONENTS program_options exception math)
if(NOT(Boost_program_options_FOUND))
  FetchContent_Declare(
    fetch_program_options
    GIT_REPOSITORY https://github.com/boostorg/program_options.git
    GIT_TAG master
    #GIT_SHALLOW TRUE
    #FIND_PACKAGE_ARGS NAMES Boost COMPONENTS program_options REQUIRED
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(fetch_program_options)
  list(APPEND METRIS_DEPS_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
else()
  list(APPEND METRIS_DEPS_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif()
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
message("-- Boost_PROGRAM_OPTIONS_LIBRARY = ${Boost_PROGRAM_OPTIONS_LIBRARY}")
message("-- Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")


include(FetchContent)
FetchContent_Declare(cmake_git_version_tracking
  GIT_REPOSITORY https://github.com/andrew-hardin/cmake-git-version-tracking.git
  GIT_TAG 6c0cb87edd029ddfb403a8e24577c144a03605a6
)
FetchContent_MakeAvailable(cmake_git_version_tracking)
set(GITLIB_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/cmake_git_version_tracking-src/")
list(APPEND METRIS_DEPS_INCLUDE_DIRS ${GITLIB_INCLUDE_DIRS})
list(APPEND METRIS_DEPS_LIBRARIES cmake_git_version_tracking)
message("-- GITLIB_INCLUDE_DIRS = ${GITLIB_INCLUDE_DIRS}")



FetchContent_MakeAvailable(${FETCH_LIST})



# This is necessary to make the sanitizer work correctly. Also we want to
# propagate flags.
if(USE_ABSL)
  setMetrisFlags(absl_hash INTERFACE)
  setMetrisFlags(absl_flat_hash_map INTERFACE)
  setMetrisFlags(absl_spinlock_wait INTERFACE)
  setMetrisFlags(absl_int128 INTERFACE)
  setMetrisFlags(absl_exponential_biased INTERFACE)
  setMetrisFlags(absl_log_severity INTERFACE)
  setMetrisFlags(absl_civil_time INTERFACE)
  setMetrisFlags(absl_raw_logging_internal INTERFACE)
  setMetrisFlags(absl_time_zone INTERFACE)
  setMetrisFlags(absl_bad_variant_access INTERFACE)
  setMetrisFlags(absl_debugging_internal INTERFACE)
  setMetrisFlags(absl_cordz_functions INTERFACE)
  setMetrisFlags(absl_bad_optional_access INTERFACE)
  setMetrisFlags(absl_throw_delegate INTERFACE)
  setMetrisFlags(absl_base INTERFACE)
  setMetrisFlags(absl_stacktrace INTERFACE)
  setMetrisFlags(absl_crc_cpu_detect INTERFACE)
  setMetrisFlags(absl_demangle_internal INTERFACE)
  setMetrisFlags(absl_string_view INTERFACE)
  setMetrisFlags(absl_city INTERFACE)
  setMetrisFlags(absl_malloc_internal INTERFACE)
  setMetrisFlags(absl_low_level_hash INTERFACE)
  setMetrisFlags(absl_strings_internal INTERFACE)
  setMetrisFlags(absl_crc_internal INTERFACE)
  setMetrisFlags(absl_graphcycles_internal INTERFACE)
  setMetrisFlags(absl_strings INTERFACE)
  setMetrisFlags(absl_hash INTERFACE)
  setMetrisFlags(absl_symbolize INTERFACE)
  setMetrisFlags(absl_time INTERFACE)
  setMetrisFlags(absl_str_format_internal INTERFACE)
  setMetrisFlags(absl_kernel_timeout_internal INTERFACE)
  setMetrisFlags(absl_crc32c INTERFACE)
  setMetrisFlags(absl_crc_cord_state INTERFACE)
  setMetrisFlags(absl_synchronization INTERFACE)
  setMetrisFlags(absl_cord_internal INTERFACE)
  setMetrisFlags(absl_cordz_handle INTERFACE)
  setMetrisFlags(absl_hashtablez_sampler INTERFACE)
  setMetrisFlags(absl_cordz_info INTERFACE)
  setMetrisFlags(absl_raw_hash_set INTERFACE)
  setMetrisFlags(absl_cord INTERFACE)
endif()


