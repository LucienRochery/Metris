
include(FindPackageHandleStandardArgs)

# Try pkg-config first
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(PC_NLOPT QUIET nlopt)
    if(PC_NLOPT_FOUND)
      message(STATUS "Found NLopt via pkg-config: ${PC_NLOPT_VERSION}")
      # pkg-config succeeded, use its results directly
      set(NLOPT_INCLUDE_DIRS ${PC_NLOPT_INCLUDE_DIRS})
      set(NLOPT_LIBRARIES ${PC_NLOPT_LIBRARIES})
      set(NLOPT_LIBRARY_DIRS ${PC_NLOPT_LIBRARY_DIRS})
      set(NLOPT_VERSION ${PC_NLOPT_VERSION})

      # Create the imported target
      if(NOT TARGET NLopt::nlopt)
        add_library(NLopt::nlopt INTERFACE IMPORTED)
        set_target_properties(NLopt::nlopt PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PC_NLOPT_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${PC_NLOPT_LIBRARIES}"
            INTERFACE_LINK_DIRECTORIES "${PC_NLOPT_LIBRARY_DIRS}")
        
        if(PC_NLOPT_CFLAGS_OTHER)
            set_target_properties(NLopt::nlopt PROPERTIES
                INTERFACE_COMPILE_OPTIONS "${PC_NLOPT_CFLAGS_OTHER}")
        endif()
        # Mark as found and return early
        find_package_handle_standard_args(NLopt
            REQUIRED_VARS PC_NLOPT_LIBRARIES PC_NLOPT_INCLUDE_DIRS
            VERSION_VAR PC_NLOPT_VERSION)
        mark_as_advanced(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)
        return()
      endif()
    endif()
endif()

## Case where PkgConfig failed

# Find include directory
find_path(NLOPT_INCLUDE_DIRS
    NAMES nlopt.h
    HINTS ${NLOPT_DIR} $ENV{NLOPT_DIR}
    PATH_SUFFIXES include
)

# Find library
find_library(NLOPT_LIBRARIES
    NAMES nlopt nlopt_cxx
    HINTS ${NLOPT_DIR} $ENV{NLOPT_DIR}
    PATH_SUFFIXES lib lib64
)

# Handle standard arguments
find_package_handle_standard_args(NLopt
    REQUIRED_VARS NLOPT_LIBRARIES NLOPT_INCLUDE_DIRS
)

if(NLopt_FOUND)
  # Create modern imported target
  if(NOT TARGET NLopt::nlopt)
    add_library(NLopt::nlopt UNKNOWN IMPORTED)
    set_target_properties(NLopt::nlopt PROPERTIES
        IMPORTED_LOCATION "${NLOPT_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${NLOPT_INCLUDE_DIRS}"
    )
    
    # Add any compile definitions from pkg-config
    if(PC_NLOPT_CFLAGS_OTHER)
      set_target_properties(NLopt::nlopt PROPERTIES
          INTERFACE_COMPILE_OPTIONS "${PC_NLOPT_CFLAGS_OTHER}"
      )
    endif()
  endif()
endif()

mark_as_advanced(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)