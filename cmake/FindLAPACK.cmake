
# Sets the variables LAPACK_INCLUDE_DIRS containing lapacke.h
# and LAPACK_LIBRARIES containing the lapack or lapacke libraries.
# If these are both set, nothing is done. 
# Otherwise, LAPACK_DIR can be provided if it contains both include/ and lib/
# Otherwise, we look for LAPACK using:
#  - C_INCLUDE_PATH
#  - CPLUS_INCLUDE_PATH
#  - default system directories

if(DEFINED LAPACK_INCLUDE_DIRS AND DEFINED LAPACK_LIBRARIES)

  if(NOT EXISTS ${LAPACK_INCLUDE_DIRS} OR NOT EXISTS ${LAPACK_INCLUDE_DIRS}/lapacke.h)
    message(FATAL_ERROR "Provided LAPACK_INCLUDE_DIRS = ${LAPACK_INCLUDE_DIRS} does not exist or does not contain lapacke.h")
  endif()

  if(NOT EXISTS ${LAPACK_LIBRARIES})
    message(FATAL_ERROR "Provided LAPACK_LIBRARIES = ${LAPACK_LIBRARIES} does not exist")
  endif()

  set(LAPACK_FOUND ON)
  message("Using provided LAPACK_INCLUDE_DIRS = ${LAPACK_INCLUDE_DIRS} " 
          "and LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")

elseif(DEFINED ENV{LAPACK_DIR} OR DEFINED LAPACK_DIR)

  if(NOT DEFINED LAPACK_DIR)
    set(LAPACK_DIR $ENV{LAPACK_DIR})
    message("Using environment variable LAPACK_DIR = $ENV{LAPACK_DIR}")
  else()
    message("Using CMake variable LAPACK_DIR = ${LAPACK_DIR}")
  endif()


  set(LAPACK_INCLUDE_DIRS "${LAPACK_DIR}/include/")
  set(LAPACK_LIBRARIES "${LAPACK_DIR}/lib/")

  if(NOT EXISTS ${LAPACK_INCLUDE_DIRS} OR NOT EXISTS ${LAPACK_INCLUDE_DIRS}/lapacke.h)
    message(FATAL_ERROR "LAPACK_INCLUDE_DIRS = ${LAPACK_INCLUDE_DIRS} does not exist or does not contain lapacke.h")
  endif()

  if(NOT EXISTS ${LAPACK_LIBRARIES})
    message(FATAL_ERROR "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES} does not exist")
  endif()

else()

  include(GNUInstallDirs) # Also handles MacOS

  string(REPLACE ":" ";" CPLUS_INCLUDE_PATH_LIST $ENV{CPLUS_INCLUDE_PATH})
  string(REPLACE ":" ";" C_INCLUDE_PATH_LIST $ENV{C_INCLUDE_PATH})

  find_path(LAPACK_INCLUDE_DIRS NAMES lapacke.h
            HINTS ${CPLUS_INCLUDE_PATH_LIST} ${C_INCLUDE_PATH_LIST} 
                  ${CMAKE_INSTALL_INCLUDEDIR} ${CMAKE_INSTALL_FULL_INCLUDEDIR}# defined by GNUInstallDirs
                  ${CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES} ${CMAKE_C_STANDARD_INCLUDE_DIRECTORIES}
            PATH_SUFFIXES include)

  if(NOT LAPACK_INCLUDE_DIRS OR NOT EXISTS ${LAPACK_INCLUDE_DIRS})
    message(FATAL_ERROR "find_path failed to find LAPACK_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR}"
      " using hint 1 ${CMAKE_INSTALL_INCLUDEDIR} 2 ${CMAKE_INSTALL_FULL_INCLUDEDIR}" 
      " 3 ${CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES}"
      " 4 ${CMAKE_C_STANDARD_INCLUDE_DIRECTORIES} 5 ${C_INCLUDE_PATH_LIST}"
      " 6 ${CPLUS_INCLUDE_PATH_LIST}")
  endif()

  if(NOT EXISTS ${LAPACK_INCLUDE_DIRS}/lapacke.h)
    message(FATAL_ERROR "LAPACK_INCLUDE_DIRS = ${LAPACK_INCLUDE_DIRS} does not not contain lapacke.h")
  endif()

  # Try first the directory in which include/lapacke.h was found
  get_filename_component(LAPACK_INCLUDE_DIR_PARENT ${LAPACK_INCLUDE_DIRS} DIRECTORY)

  find_library(LAPACK_LIBRARIES NAMES lapacke lapack 
               HINTS ${LAPACK_INCLUDE_DIR_PARENT} 
                     ${CMAKE_INSTALL_LIBDIR} # defined by GNUInstallDirs
               PATH_SUFFIXES lib)

  if(NOT LAPACK_LIBRARIES OR NOT EXISTS ${LAPACK_LIBRARIES})
    message(FATAL_ERROR "find_path failed to find lapacke.h with hint ${CMAKE_INSTALL_INCLUDEDIR}")
  endif()

endif()

set(LAPACK_FOUND ON)
set(LAPACKE_FOUND ON)