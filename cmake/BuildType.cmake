#==================================================
# Default compiler flags, these can be modified under
# the advanced options using ccmake
# Taken from Solution Adaptive Numerical Simulator (SANS) with permission:
# https://darmofal.mit.edu/solution-adaptive-numerical-simulator-sans/
#==================================================
IF( NOT CMAKE_BUILD_TYPE )

  #===============================
  # Set the build type to release by default, but debug if the binary directory contains the name debug
  SET( BUILD_TYPE_STRING "Options are: Debug Release RelWithDebInfo Memcheck." )

  SET( BIN_NAMES )
  IF( BIN_DIR_NAME MATCHES "DEBUG" )
    SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    LIST(APPEND BIN_NAMES ${CMAKE_BUILD_TYPE} )
  ENDIF()
  IF( BIN_DIR_NAME MATCHES "RELEASE" OR BIN_DIR_NAME MATCHES "DEPLOY" )
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    LIST(APPEND BIN_NAMES ${CMAKE_BUILD_TYPE} )
  ENDIF()
  IF( BIN_DIR_NAME MATCHES "RELWITHDEBINFO" )
    SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    LIST(APPEND BIN_NAMES ${CMAKE_BUILD_TYPE} )
  ENDIF()
  IF( BIN_DIR_NAME MATCHES "MEMCHECK" )
    SET(CMAKE_BUILD_TYPE "Memcheck" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    LIST(APPEND BIN_NAMES ${CMAKE_BUILD_TYPE} )
  ENDIF()
  LIST(LENGTH BIN_NAMES BIN_NAME_COUNT )
  IF( BIN_NAME_COUNT GREATER 1)
    MESSAGE( "" )
    MESSAGE( "=================================" )
    MESSAGE( "Build directory name should contain only one of: ${BIN_NAMES}")
    MESSAGE( "${BUILD_TYPE_STRING}" )
    MESSAGE( "=================================" )
    MESSAGE( "" )
    MESSAGE(FATAL_ERROR "" )
  ENDIF()

  IF( BIN_NAME_COUNT LESS 1)
    SET(CMAKE_BUILD_TYPE "Release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
  ENDIF()

ENDIF()