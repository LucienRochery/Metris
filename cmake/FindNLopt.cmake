#
#Find the NLopt includes and libraries
#
#Code taken from SANS https://darmofal.mit.edu/solution-adaptive-numerical-simulator-sans/

IF(NOT DEFINED NLOPT_DIR)
  PKG_CHECK_MODULES(NLOPT QUIET libnlopt)

  IF(NLOPT_FOUND)
    SET(NLOPT_DEFINITIONS ${NLOPT_CFLAGS_OTHER})
  ELSE()
    FIND_PATH(NLOPT_INCLUDE_DIRS nlopt.h
                  HINTS $ENV{NLOPT_DIR} /
              include
                  PATH_SUFFIXES nlopt)

    FIND_LIBRARY(NLOPT_LIBRARIES NAMES nlopt nlopt_cxx
                    HINTS $ENV{NLOPT_DIR} /
                lib)
  ENDIF()
ELSE()
  FIND_PATH(NLOPT_INCLUDE_DIRS nlopt.h
                HINTS ${NLOPT_DIR} /
            include
                PATH_SUFFIXES nlopt NO_DEFAULT_PATH)

  FIND_LIBRARY(NLOPT_LIBRARIES NAMES nlopt nlopt_cxx
                  HINTS ${NLOPT_DIR} /
              lib NO_DEFAULT_PATH)
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(NLOPT DEFAULT_MSG
                                      NLOPT_LIBRARIES NLOPT_INCLUDE_DIRS)

MARK_AS_ADVANCED(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)
