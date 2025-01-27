include(FetchContent)

include(MetrisFlags)

 
if(REQ_CODEGEN)
  if(NOT DEFINED GINAC_LIBRARIES)
  set(GINAC_LIBRARIES "-L$ENV{GINAC_DIR}/lib/")
  endif()
  if(NOT DEFINED GINAC_INCLUDE_DIRS)
  set(GINAC_INCLUDE_DIRS "$ENV{GINAC_DIR}/include")
  endif()
  if(NOT DEFINED CLN_INCLUDE_DIRS)
  set(CLN_INCLUDE_DIRS "$ENV{CLN_DIR}/include")
  endif()

  message("GINAC_INCLUDE_DIRS = ${GINAC_INCLUDE_DIRS}")
  message("CLN_INCLUDE_DIRS = ${CLN_INCLUDE_DIRS}")
  message("GINAC_LIBRARIES = ${GINAC_LIBRARIES}")

  set(GINAC_LIBRARIES "${GINAC_LIBRARIES} -lginac")

  set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN_DIR}")
  if (NOT EIGEN3_INCLUDE_DIR)
    FetchContent_Declare(
      Eigen3
      #URL https://github.com/abseil/abseil-cpp/archive/e7fe9ec9ebfc6607765d489b76c9954e0a88c5d4.zip
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG d6e3b528b2ae2a55d55749e9ef50b1e12ff34bc3
      #GIT_TAG master
      #FIND_PACKAGE_ARGS NAMES Eigen3
      EXCLUDE_FROM_ALL
    )
    LIST(APPEND FETCH_LIST Eigen3)
    set(EIGEN3_INCLUDE_DIR "${CMAKE_BINARY_DIR}/_deps/Eigen3-src/")
    message("EIGEN3_INCLUDE_DIR = ${EIGEN3_INCLUDE_DIR}")
  endif()
  include_directories( ${EIGEN3_INCLUDE_DIR} )
  #set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN_DIR}" )
  #if (NOT EIGEN3_INCLUDE_DIR)
  #  find_package(Eigen3 REQUIRED)
  #endif()
  #message("EIGEN3_INCLUDE_DIR = ${EIGEN3_INCLUDE_DIR}")
  #message("EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}")
  #include_directories( ${EIGEN3_INCLUDE_DIR} )
endif()



if(NOT(DEFINED ENV{LAPACK_INCLUDE_DIR}))
  message(FATAL_ERROR "Set the environment variable LAPACK_INCLUDE_DIR to the folder containing lapacke.h")
else()
  set(LAPACK_INCLUDE_DIR $ENV{LAPACK_INCLUDE_DIR})
  message("Found LAPACK_INCLUDE_DIR = ${LAPACK_INCLUDE_DIR}")
endif()


if(USE_PETSC STREQUAL "True" 
OR USE_PETSC STREQUAL "ON"
OR USE_PETSC STREQUAL "Yes")
  message(STATUS "PETSC enabled")

  if(PETSC_FOUND AND NOT PETSC_PKGCONFIG)

    message(STATUS "PETSC_FOUND = ${PETSC_FOUND}")
    message("PETSC_LIBRARIES = ${PETSC_LIBRARIES}")

  else()

    if(DEFINED PETSC_DIR AND DEFINED PETSC_ARCH)

      message("PETSC_DIR = ${PETSC_DIR} and PETSC_ARCH = ${PETSC_ARCH} already set: assuming PETSc target already defined")
      #find_package(MPI)
      set(PETSC ${PETSC_DIR}/${PETSC_ARCH})
      set(PETSC_INCLUDE ${PETSC_DIR}/include)
      list(APPEND PETSC_INCLUDE ${PETSC_DIR}/${PETSC_ARCH}/include)
      list(APPEND PETSC_INCLUDE ${MPI_INCLUDE_PATH})

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

      set(PETSC ${PETSC_DIR}/${PETSC_ARCH})
      set(PETSC_INCLUDE ${PETSC_DIR}/include)

      message("Debug PETSC = ${PETSC} PETSC_INCLUDE = ${PETSC_INCLUDE}")
      list(APPEND PETSC_INCLUDE ${PETSC_DIR}/${PETSC_ARCH}/include)
      list(APPEND PETSC_INCLUDE ${MPI_INCLUDE_PATH})
      
      set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
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
    
  endif()
else()
  message("PETSC disabled")
  set(PETSC_INCLUDE "")
endif()




 
if(NOT(DEFINED ENV{ESP_ROOT}))
  message(FATAL_ERROR "Set the environment variable ESP_ROOT to the folder containing include/egads.h")
else()
  set(ESP_ROOT $ENV{ESP_ROOT})
  message("Found ESP_ROOT = ${ESP_ROOT}")
endif()
add_library(libegads SHARED IMPORTED GLOBAL)
add_library(libegadslite SHARED IMPORTED GLOBAL)
find_file(EGADS_LIBRARY NAMES libegads.dylib libegads.so PATHS ${ESP_ROOT}/lib/)
find_file(EGADSLITE_LIBRARY NAMES libegadslite.dylib libegadslite.so PATHS ${ESP_ROOT}/lib/)
message("Using EGADS_LIBRARY = ${EGADS_LIBRARY}")
message("Using EGADSLITE_LIBRARY = ${EGADSLITE_LIBRARY}")
set_target_properties(libegads     PROPERTIES IMPORTED_LOCATION ${EGADS_LIBRARY})
set_target_properties(libegadslite PROPERTIES IMPORTED_LOCATION ${EGADSLITE_LIBRARY})


if(USE_CLP STREQUAL "True" 
OR USE_CLP STREQUAL "ON")
  include(FindCLP)
  if(NOT(CLP_FOUND))
    message(WARNING "CLP was not found on this system.")
    set(CLP_INCLUDE_DIRS "")
    set(CLP_LIBRARIES "")
  else()
    add_compile_definitions(USE_CLP)
  endif()
endif()

# External libraries to be fetched


FetchContent_Declare(
  fetch_absl
  #URL https://github.com/abseil/abseil-cpp/archive/e7fe9ec9ebfc6607765d489b76c9954e0a88c5d4.zip
  GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
  GIT_TAG e7fe9ec9ebfc6607765d489b76c9954e0a88c5d4  
  #GIT_TAG master
  #FIND_PACKAGE_ARGS NAMES absl
  EXCLUDE_FROM_ALL
)
LIST(APPEND FETCH_LIST fetch_absl)
set(ABSL_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/fetch_absl-src/")
message("CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
message("ABSL_INCLUDE_DIRS = ${ABSL_INCLUDE_DIRS}")


#FetchContent_Declare(
#  fetch_boost_hana
#  GIT_REPOSITORY https://github.com/boostorg/hana.git
#  GIT_TAG master   
#  GIT_SHALLOW TRUE
#  FIND_PACKAGE_ARGS NAMES Boost COMPONENTS hana
#  EXCLUDE_FROM_ALL
#)
#LIST(APPEND FETCH_LIST fetch_boost_hana)

#FetchContent_Declare(
#  fetch_nlopt
#  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
#  GIT_TAG master   
#  GIT_SHALLOW TRUE
#  FIND_PACKAGE_ARGS NAMES NLopt 
#  EXCLUDE_FROM_ALL
#)
#LIST(APPEND FETCH_LIST fetch_nlopt)

find_package(Boost COMPONENTS program_options)
if(NOT(Boost_program_options_FOUND))
  FetchContent_Declare(
    fetch_program_options
    GIT_REPOSITORY https://github.com/boostorg/program_options.git
    GIT_TAG master   
    #GIT_SHALLOW TRUE
    #FIND_PACKAGE_ARGS NAMES Boost COMPONENTS program_options REQUIRED
    EXCLUDE_FROM_ALL
  )
  LIST(APPEND FETCH_LIST fetch_program_options)
endif()


#if(USE_TRACELIBS)
#  find_package(Boost COMPONENTS stacktrace_basic  
#                                #stacktrace_backtrace  
#                                stacktrace_addr2line
#                                stacktrace_noop REQUIRED)
#  set(BOOST_TRACELIBS ${Boost_STACKTRACE_BASIC_LIBRARY} 
#                      ${Boost_STACKTRACE_BACKTRACE_LIBRARY} 
#                      ${Boost_STACKTRACE_ADDR2LINE_LIBRARY} 
#                      ${Boost_STACKTRACE_NOOP_LIBRARY} 
#                      pthread dl)
#  
#  find_program(BOOST_STACKTRACE_ADDR2LINE_LOCATION addr2line)
#  if(BOOST_STACKTRACE_ADDR2LINE_LOCATION)
#    message("Found addr2line at ${BOOST_STACKTRACE_ADDR2LINE_LOCATION}")
#    add_compile_definitions(BOOST_STACKTRACE_USE_ADDR2LINE)
#  else()
#    message("Executable addr2line not found.")
#  endif()
#
#else()
  set(BOOST_TRACELIBS "")
#endif()

#FetchContent_Declare(
#  fetch_lapack
#  GIT_REPOSITORY https://github.com/Reference-LAPACK/lapack.git
#  GIT_TAG master   
#  GIT_SHALLOW TRUE
#  FIND_PACKAGE_ARGS NAMES LAPACK
#  EXCLUDE_FROM_ALL
#)
#LIST(APPEND FETCH_LIST fetch_lapack)

#if(USE_PETSC STREQUAL "True")
#  FetchContent_Declare(
#    fetch_PETSc
#    GIT_REPOSITORY https://github.com/petsc/petsc.git
#    GIT_TAG main   
#    GIT_SHALLOW TRUE
#    EXCLUDE_FROM_ALL
#  )
#  FetchContent_MakeAvailable(fetch_PETSc)
#  LIST(APPEND FETCH_LIST fetch_PETSc)
#endif()

#FetchContent_Declare(
#  fetch_boost_stacktrace
#  GIT_REPOSITORY https://github.com/boostorg/stacktrace.git
#  GIT_TAG master   
#  GIT_SHALLOW TRUE
#  FIND_PACKAGE_ARGS NAMES Boost COMPONENTS stacktrace_basic
#                                           stacktrace_addr2line
#                                           stacktrace_noop
#  EXCLUDE_FROM_ALL
#)
#LIST(APPEND FETCH_LIST fetch_boost_stacktrace)


FetchContent_MakeAvailable(${FETCH_LIST})
# This is necessary to make the sanitizer work correctly. Also we do want to 
# propagate flags, in general. 
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





