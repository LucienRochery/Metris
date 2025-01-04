# Credit : https://github.com/libigl/CoMISo/blob/master/cmake/FindCLP.cmake
# - Try to find CLP
# Once done this will define
#  CLP_FOUND - System has CLP
#  CLP_INCLUDE_DIRS - The CLP include directories
#  CLP_LIBRARIES - The libraries needed to use CLP

if(NOT(DEFINED ENV{CLP_DIR}))
  message(WARNING "Environment variable CLP_DIR not set. Searching common directories. On MacOS with brew: /opt/homebrew/opt/clp/")
endif()

if(NOT(DEFINED ENV{COINUTILS_DIR}))
  message(WARNING "If fails: set the environment variable COINUTILS_DIR. On MacOS with brew: /opt/homebrew/opt/coinutils/")
endif()

## I8 Search paths for windows libraries
#if ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*Win64" )
#  SET(VS_SEARCH_PATH "c:/libs/vs2012/x64/")
#elseif ( CMAKE_GENERATOR MATCHES "^Visual Studio 11.*" )
#  SET(VS_SEARCH_PATH "c:/libs/vs2012/x32/")
#elseif ( CMAKE_GENERATOR MATCHES "^Visual Studio 12.*Win64" )
#  SET(VS_SEARCH_PATH "c:/libs/vs2013/x64/")
#elseif ( CMAKE_GENERATOR MATCHES "^Visual Studio 12.*" )
#  SET(VS_SEARCH_PATH "c:/libs/vs2013/x32/")
#endif()


if (CLP_INCLUDE_DIR AND COINUTILS_INCLUDE_DIR)
  # in cache already
  set(CLP_FOUND TRUE)
  set(CLP_INCLUDE_DIRS ${CLP_INCLUDE_DIR} )
  list(APPEND CLP_INCLUDE_DIRS ${COINUTILS_INCLUDE_DIR})
  set(CLP_LIBRARIES "${CLP_LIBRARY}" )
  list(APPEND CLP_LIBRARIES ${COINUTILS_LIBRARY})
else()

  find_path(CLP_INCLUDE_DIR 
            NAMES ClpConfig.h
            PATHS "$ENV{CLP_DIR}/include/clp/coin"
                   "/usr/include/clp/coin"
				   "/opt/homebrew/opt/clp/include/clp/coin"
           #        "C:\\libs\\clp\\include"
  				 #"${VS_SEARCH_PATH}CBC-2.9.4/Clp/include"
           )
  find_path(COINUTILS_INCLUDE_DIR
            NAMES CoinPragma.hpp
            PATHS "$ENV{COINUTILS_DIR}/include/coinutils/coin"
                   "/usr/include/coinutils/coin"
				   "/opt/homebrew/opt/coinutils/include/coinutils/coin"
            )

  find_library( CLP_LIBRARY 
                NAMES Clp libClp
                PATHS "$ENV{CLP_DIR}/lib"
                      "/usr/lib"
                      "/usr/lib/coin"
					  "/opt/homebrew/opt/clp/lib"
             #          "C:\\libs\\clp\\lib"
  					#"${VS_SEARCH_PATH}CBC-2.9.4/Clp/lib"
            )

  find_library( COINUTILS_LIBRARY 
                NAMES CoinUtils
                PATHS "$ENV{COINUTILS_DIR}/lib"
                      "/usr/lib"
                      "/usr/lib/coin"
					  "/opt/homebrew/opt/coinutils/lib"
             #          "C:\\libs\\clp\\lib"
            #"${VS_SEARCH_PATH}CBC-2.9.4/Clp/lib"
            )
			

  set(CLP_INCLUDE_DIRS "${CLP_INCLUDE_DIR}" )
  list(APPEND CLP_INCLUDE_DIRS ${COINUTILS_INCLUDE_DIR})
  set(CLP_LIBRARIES "${CLP_LIBRARY}" )
  list(APPEND CLP_LIBRARIES ${COINUTILS_LIBRARY})

  message("Found CLP_INCLUDE_DIRS = ${CLP_INCLUDE_DIRS}" )
  message("Found CLP_LIBRARIES = ${CLP_LIBRARIES}" )


  include(FindPackageHandleStandardArgs)
  # handle the QUIETLY and REQUIRED arguments and set CLP_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(CLP  DEFAULT_MSG
                                    CLP_LIBRARY CLP_INCLUDE_DIR)

  mark_as_advanced(CLP_INCLUDE_DIR CLP_LIBRARY)

endif()

