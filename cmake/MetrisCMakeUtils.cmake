# CMake utility functions and settings for Metris
# Centralizes common CMake configuration

# ============================================================================
# RPATH Configuration
# ============================================================================
# Centralized RPATH handling for better cross-platform support

function(metris_setup_rpath)
  # Enable RPATH on macOS
  if(APPLE)
    set(CMAKE_MACOSX_RPATH 1 PARENT_SCOPE)
  endif()
  
  # Use install RPATH (modern approach)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE PARENT_SCOPE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE PARENT_SCOPE)
  
  # Set install RPATH to include library directory
  if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}" PARENT_SCOPE)
  else()
    set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}" PARENT_SCOPE)
  endif()
  
  # Add any additional RPATH entries
  if(DEFINED METRIS_ADDITIONAL_RPATH)
    list(APPEND CMAKE_INSTALL_RPATH ${METRIS_ADDITIONAL_RPATH})
    set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH} PARENT_SCOPE)
  endif()
endfunction()

# ============================================================================
# Version Requirements
# ============================================================================

# Set minimum version requirements for dependencies
# These can be overridden by find_package calls with explicit versions
set(METRIS_DEPENDENCY_VERSIONS
  "Boost:1.70"  # Minimum Boost version
  "CMake:3.14"  # Already enforced in main CMakeLists.txt
  CACHE INTERNAL "Minimum dependency versions"
)

# ============================================================================
# Boost Modernization Helpers
# ============================================================================

# Helper to get Boost libraries as targets
# Prefers modern target-based approach over old-style variables
function(metris_get_boost_libraries output_var)
  set(boost_libs "")
  
  # Use modern target-based approach if available
  if(TARGET Boost::program_options)
    list(APPEND boost_libs Boost::program_options)
  endif()
  
  # Fallback to old-style variables for compatibility
  if(boost_libs STREQUAL "" AND DEFINED Boost_LIBRARIES)
    set(boost_libs ${Boost_LIBRARIES})
  endif()
  
  # Always include headers target if available
  if(TARGET Boost::headers)
    list(APPEND boost_libs Boost::headers)
  endif()
  
  set(${output_var} ${boost_libs} PARENT_SCOPE)
endfunction()

# ============================================================================
# Install Configuration Helpers
# ============================================================================

# Set up proper install RPATH for a target
function(metris_setup_target_rpath target_name)
  set_target_properties(${target_name} PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
    BUILD_WITH_INSTALL_RPATH TRUE
  )
  
  # Set install RPATH
  if(APPLE)
    set_target_properties(${target_name} PROPERTIES
      INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}"
    )
  else()
    set_target_properties(${target_name} PROPERTIES
      INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    )
  endif()
endfunction()

# ============================================================================
# Dependency Version Checking
# ============================================================================

# Check if a dependency version meets requirements
function(metris_check_dependency_version dep_name found_version min_version)
  if(found_version AND min_version)
    if(found_version VERSION_LESS min_version)
      message(WARNING 
        "${dep_name} version ${found_version} is less than required ${min_version}. "
        "This may cause compatibility issues."
      )
      return()
    endif()
    message(STATUS "${dep_name} version ${found_version} meets requirement (>= ${min_version})")
  endif()
endfunction()

