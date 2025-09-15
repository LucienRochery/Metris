
# Defines Metris_FOUND
#         METRIS_LIBRARIES
#         METRIS_INCLUDE_DIRS

# Requires - METRIS_DIR as CMake or environment variable
#            e.g. /path/to/Metris/ <- contains build/, src/ etc.
#          - BIN_DIR_NAME as CMake variable (e.g. release_clang)

# Searches for installed (make install) libraries in 
# $METRIS_DIR/build/$BIN_DIR_NAME/install 

set(Metris_FOUND FALSE)

if(NOT DEFINED METRIS_DIR)
  if(DEFINED ENV{METRIS_DIR})
    set(METRIS_DIR $ENV{METRIS_DIR} CACHE PATH "Path to Metris installation")
    message(STATUS "Using Metris from environment variable METRIS_DIR: ${METRIS_DIR}")
  else()
    message(STATUS "Undefined METRIS_DIR") 
    return()
  endif()
endif()


string(TOLOWER ${BIN_DIR_NAME} BIN_DIR_NAME_LOWER)

# First try to find cmake config files
find_path(METRIS_CMAKE_LIB_DIR 
          NAMES MetrisConfig.cmake
          PATHS ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install/lib/cmake/Metris
          NO_DEFAULT_PATH)

if(METRIS_CMAKE_LIB_DIR)
  message(STATUS "Found Metris CMake config in ${METRIS_CMAKE_LIB_DIR}")
  set(Metris_DIR ${METRIS_CMAKE_LIB_DIR})
  find_package(Metris CONFIG REQUIRED)
  if(Metris_FOUND)
    set(METRIS_LIBRARIES Metris::libMetris)
    get_target_property(METRIS_INCLUDE_DIRS Metris::libMetris INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Using Metris from CMake config.")
  else()
    message(WARNING "Metris found but not configured correctly")
    return()
  endif()
else()
  # This path is not optimal, as CMake won't know about any Metris
  # dependencies. This should probably never be used.
  message(WARNING "Metris CMake config not found, falling back to library search")
  message(WARNING "CMake won't know about Metris dependencies!")
  message(STATUS "Searching in ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install/lib")
  find_library(METRIS_LIBRARIES
    NAMES libMetris libMetris
    PATHS ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install/lib
    NO_DEFAULT_PATH
  )
  find_path(METRIS_INCLUDE_DIRS
    NAMES Metris.h
    PATHS ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install/include/src
    NO_DEFAULT_PATH
  )
endif()

if(METRIS_LIBRARIES)
  set(Metris_FOUND TRUE)
  message(STATUS "Found Metris in ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install")
else()
  message("## METRIS NOT FOUND in ${METRIS_DIR}/build/${BIN_DIR_NAME_LOWER}/install")
endif()

