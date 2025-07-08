#==================================================
# Fetches the compiler from directory name
# Taken from Solution Adaptive Numerical Simulator (SANS) with permission:
# https://darmofal.mit.edu/solution-adaptive-numerical-simulator-sans/
#==================================================

SET( GNU_VERSIONS "4.8" "4.9" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" )
SET( CLANG_VERSIONS "3.4" "3.5" "3.6" "3.7" "3.8" "3.9" "4.0" "5.0" "6.0" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" )

SET( BIN_NAMES )

IF( BIN_DIR_NAME MATCHES "CLANG" )

  #Set to default version if no version is given
  SET( CLANG_C_COMPILER clang )
  SET( CLANG_CXX_COMPILER clang++ )

  #Choose a specific version if specified
  FOREACH( VERSION ${CLANG_VERSIONS} )
    STRING( REPLACE "." "" CLANGVERSION "CLANG${VERSION}" )
    IF( BIN_DIR_NAME MATCHES ${CLANGVERSION} )
      SET( CLANG_C_COMPILER clang-${VERSION} )
      SET( CLANG_CXX_COMPILER clang++-${VERSION} )
    ENDIF()
  ENDFOREACH()

  FIND_PROGRAM( CLANG_C_COMPILER_FOUND ${CLANG_C_COMPILER} )
  FIND_PROGRAM( CLANG_CXX_COMPILER_FOUND ${CLANG_CXX_COMPILER} )


  #If the specific version was not found, look for generic version
  IF( NOT CLANG_C_COMPILER_FOUND OR NOT CLANG_CXX_COMPILER_FOUND )
    message(WARNING "CLANG not found")
    UNSET( CLANG_C_COMPILER )
    UNSET( CLANG_CXX_COMPILER )
    FIND_PROGRAM( CLANG_C_COMPILER clang )
    FIND_PROGRAM( CLANG_CXX_COMPILER clang++ )
  else()
    SET(CLANG_C_COMPILER ${CLANG_C_COMPILER_FOUND})
    SET(CLANG_CXX_COMPILER ${CLANG_CXX_COMPILER_FOUND})
  ENDIF()

  UNSET( CLANG_C_COMPILER_FOUND CACHE )
  UNSET( CLANG_CXX_COMPILER_FOUND CACHE )

  SET( CMAKE_C_COMPILER ${CLANG_C_COMPILER}  )
  SET( CMAKE_CXX_COMPILER ${CLANG_CXX_COMPILER} )
  LIST(APPEND BIN_NAMES "Clang" )

ENDIF()

# Also check for generators
IF( BIN_DIR_NAME MATCHES "XCODE" )
  LIST(APPEND BIN_NAMES "Xcode" )
ENDIF()


IF( BIN_DIR_NAME MATCHES "INTEL" )
  # First look for the Intel(R) oneAPI DPC++/C++ Compiler (ICX)
  FIND_PROGRAM( INTEL_DPC_C_COMPILER icx )
  FIND_PROGRAM( INTEL_DPC_CXX_COMPILER icpx )
  IF( INTEL_DPC_C_COMPILER AND INTEL_DPC_CXX_COMPILER )
    message("Found IntelLLVM compilers C = ${INTEL_DPC_C_COMPILER}, C++ = ${INTEL_DPC_CXX_COMPILER}")
    SET( CMAKE_C_COMPILER ${INTEL_DPC_C_COMPILER} )
    SET( CMAKE_CXX_COMPILER ${INTEL_DPC_CXX_COMPILER} )
  ELSE()
    # Use classic intel compilers
    FIND_PROGRAM( INTEL_C_COMPILER icc )
    FIND_PROGRAM( INTEL_CXX_COMPILER icpc )
    message("Found classic Intel compilers C = ${INTEL_DPC_C_COMPILER}, C++ = ${INTEL_DPC_CXX_COMPILER}")
    SET( CMAKE_C_COMPILER ${INTEL_C_COMPILER} )
    SET( CMAKE_CXX_COMPILER ${INTEL_CXX_COMPILER}  )
  ENDIF()
  LIST(APPEND BIN_NAMES "Intel" )

  if(NOT INTEL_DPC_C_COMPILER OR NOT INTEL_DPC_CXX_COMPILER)
    message(FATAL_ERROR "Failed to find either IntelLLVM or Classic Intel compilers")
  endif()

  
  UNSET( INTEL_DPC_C_COMPILER CACHE )
  UNSET( INTEL_DPC_CXX_COMPILER CACHE )

ENDIF()


IF( BIN_DIR_NAME MATCHES "GNU" OR BIN_DIR_NAME MATCHES "GCC")
  #Set to default version if no version is given
  SET( GNU_C_COMPILER gcc )
  SET( GNU_CXX_COMPILER g++ )

  #Choose a specific version if specified
  FOREACH( VERSION ${GNU_VERSIONS} )
    STRING( REPLACE "." "" GNUVERSION "GNU${VERSION}" )
    IF( BIN_DIR_NAME MATCHES ${GNUVERSION} )
      SET( GNU_C_COMPILER gcc-${VERSION} )
      SET( GNU_CXX_COMPILER g++-${VERSION} )
      SET( GCOV_COMMAND gcov-${VERSION} )
    ENDIF()
  ENDFOREACH()

  FIND_PROGRAM( GNU_C_COMPILER_FOUND ${GNU_C_COMPILER} )
  FIND_PROGRAM( GNU_CXX_COMPILER_FOUND ${GNU_CXX_COMPILER} )

  #If the specific version was not find, look for generic version
  IF( NOT GNU_C_COMPILER_FOUND OR NOT GNU_CXX_COMPILER_FOUND )
    MESSAGE("## Could not set requested GCC version as it was not found in the system: ${GNU_C_COMPILER} or ${GNU_CXX_COMPILER}")
    UNSET( GNU_C_COMPILER )
    UNSET( GNU_CXX_COMPILER )
    FIND_PROGRAM( GNU_C_COMPILER gnu )
    FIND_PROGRAM( GNU_CXX_COMPILER g++ )
  else()
    SET(GNU_C_COMPILER ${GNU_C_COMPILER_FOUND})
    SET(GNU_CXX_COMPILER ${GNU_CXX_COMPILER_FOUND})
  ENDIF()

  UNSET( GNU_C_COMPILER_FOUND CACHE )
  UNSET( GNU_CXX_COMPILER_FOUND CACHE )

  SET( CMAKE_C_COMPILER ${GNU_C_COMPILER}  )
  SET( CMAKE_CXX_COMPILER ${GNU_CXX_COMPILER} )
  
  LIST(APPEND BIN_NAMES "gnu" )

ENDIF()


LIST(LENGTH BIN_NAMES BIN_NAME_COUNT )
IF( BIN_NAME_COUNT GREATER 1)
  MESSAGE( "" )  
  MESSAGE( "=================================" )
  MESSAGE( "Build directory name should contain only one of: ${BIN_NAMES}")
  MESSAGE( "=================================" )
  MESSAGE( "" )  
  MESSAGE(FATAL_ERROR "" )  
ENDIF()

IF( BIN_NAME_COUNT EQUAL 0)
  message("")
  MESSAGE( "=================================" )
  MESSAGE( "Build directory name containts no recognized compiler name, using CMake's system default")
  MESSAGE( "=================================" )
  MESSAGE( "" )  
ENDIF()
