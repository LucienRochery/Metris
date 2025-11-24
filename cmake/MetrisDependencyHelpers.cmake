# Helper functions for managing dependencies in Metris
# These functions provide a consistent interface for finding or fetching dependencies
# and properly registering them for export in MetrisConfig.cmake

include(FetchContent)

# Improved dependency registration using structured format
# This replaces the string-based format for better maintainability
# Note: This must be set in the parent scope, not as a cache variable
# The parent scope should initialize: set(METRIS_CONFIG_DEPENDENCIES "")

function(metris_register_dependency IMPORT_TYPE DEP_NAME DEP_COMPONENTS)
  # Store using pipe separator to avoid CMake list splitting issues
  # Format: "TYPE|NAME|COMPONENTS"
  # The parsing code in CMakeLists.txt handles this format
  # Use a cache-backed internal variable so that entries added from nested
  # function calls (e.g., metris_find_eigen3() calling this function) are
  # persisted and visible in the caller scopes. Using the cache avoids
  # tricky scoping when functions call other functions which then register
  # dependencies.
  if(NOT DEFINED METRIS_CONFIG_DEPENDENCIES)
    set(METRIS_CONFIG_DEPENDENCIES "" CACHE INTERNAL "Metris dependencies list (internal)")
  endif()
  # Read current cached value, append our new entry, and write it back to cache.
  set(_metris_deps_cached "${METRIS_CONFIG_DEPENDENCIES}")
  # Build the canonical entry string
  set(_new_entry "${IMPORT_TYPE}|${DEP_NAME}|${DEP_COMPONENTS}")
  if(_metris_deps_cached STREQUAL "")
    set(_metris_deps_cached "${_new_entry}")
  else()
    # Deduplicate: only append if not already present
    # Convert cached value into a list and check membership
    list(APPEND _tmp_list ${_metris_deps_cached})
    list(FIND _tmp_list "${_new_entry}" _found_index)
    if(_found_index EQUAL -1)
      list(APPEND _metris_deps_cached "${_new_entry}")
    endif()
    unset(_tmp_list)
  endif()
  # Store back into cache so top-level code can read it after includes
  set(METRIS_CONFIG_DEPENDENCIES "${_metris_deps_cached}" CACHE INTERNAL "Metris dependencies list (internal)")
  unset(_metris_deps_cached)
endfunction()

# Helper function to find or fetch a dependency
# Usage:
#   metris_find_or_fetch_dependency(
#     NAME <name>
#     [VERSION <version>]
#     [COMPONENTS <components>]
#     [GIT_REPOSITORY <url>]
#     [GIT_TAG <tag>]
#     [URL <url>]
#     [TARGET <target_name>]  # Expected target name (e.g., "fmt::fmt")
#     [INCLUDE_VAR <var_name>]  # Variable to set with include directory
#     [REQUIRED]  # Make dependency required
#   )
function(metris_find_or_fetch_dependency)
  set(options REQUIRED)
  set(oneValueArgs NAME VERSION GIT_REPOSITORY GIT_TAG URL TARGET INCLUDE_VAR)
  set(multiValueArgs COMPONENTS)
  cmake_parse_arguments(METRIS_DEP "${options}" "${oneValue_args}" "${multiValue_args}" ${ARGN})

  if(NOT METRIS_DEP_NAME)
    message(FATAL_ERROR "metris_find_or_fetch_dependency: NAME is required")
  endif()

  # First check if target already exists (from parent project)
  if(METRIS_DEP_TARGET AND TARGET ${METRIS_DEP_TARGET})
    message(STATUS "Using ${METRIS_DEP_TARGET} target provided by parent project")
    target_link_libraries(metris_deps INTERFACE ${METRIS_DEP_TARGET})
    set(components_str "")
    if(METRIS_DEP_COMPONENTS)
      string(REPLACE ";" " " components_str "${METRIS_DEP_COMPONENTS}")
    endif()
    metris_register_dependency("find_package" "${METRIS_DEP_NAME}${METRIS_DEP_REQUIRED}" "${components_str}")
    return()
  endif()

  # Try find_package
  if(METRIS_DEP_VERSION)
    find_package(${METRIS_DEP_NAME} ${METRIS_DEP_VERSION} QUIET ${METRIS_DEP_COMPONENTS})
  else()
    find_package(${METRIS_DEP_NAME} QUIET ${METRIS_DEP_COMPONENTS})
  endif()

  # Check again for target after find_package (it might have been created)
  if(${METRIS_DEP_NAME}_FOUND OR (METRIS_DEP_TARGET AND TARGET ${METRIS_DEP_TARGET}))
    # Found via find_package
    message(STATUS "Found ${METRIS_DEP_NAME} via find_package")
    
    # Register for export
    set(components_str "")
    if(METRIS_DEP_COMPONENTS)
      string(REPLACE ";" " " components_str "${METRIS_DEP_COMPONENTS}")
    endif()
    metris_register_dependency("find_package" "${METRIS_DEP_NAME}${METRIS_DEP_REQUIRED}" "${components_str}")
    
    # Link to metris_deps if target exists
    if(METRIS_DEP_TARGET AND TARGET ${METRIS_DEP_TARGET})
      target_link_libraries(metris_deps INTERFACE ${METRIS_DEP_TARGET})
    endif()
    
    return()
  endif()

  # Not found, try FetchContent
  if(METRIS_DEP_GIT_REPOSITORY)
    message(STATUS "${METRIS_DEP_NAME} not found via find_package, fetching from ${METRIS_DEP_GIT_REPOSITORY}")
    FetchContent_Declare(
      ${METRIS_DEP_NAME}_fetch
      GIT_REPOSITORY ${METRIS_DEP_GIT_REPOSITORY}
      GIT_TAG ${METRIS_DEP_GIT_TAG}
      EXCLUDE_FROM_ALL
    )
  elseif(METRIS_DEP_URL)
    message(STATUS "${METRIS_DEP_NAME} not found via find_package, fetching from ${METRIS_DEP_URL}")
    FetchContent_Declare(
      ${METRIS_DEP_NAME}_fetch
      URL ${METRIS_DEP_URL}
      EXCLUDE_FROM_ALL
    )
  else()
    if(METRIS_DEP_REQUIRED)
      message(FATAL_ERROR "${METRIS_DEP_NAME} not found and no FetchContent source provided")
    else()
      message(WARNING "${METRIS_DEP_NAME} not found and no FetchContent source provided")
      return()
    endif()
  endif()

  FetchContent_MakeAvailable(${METRIS_DEP_NAME}_fetch)
  
  # Register for export
  metris_register_dependency("FetchContent" "${METRIS_DEP_NAME}" "")
  
  # Link to metris_deps if target exists
  if(METRIS_DEP_TARGET AND TARGET ${METRIS_DEP_TARGET})
    target_link_libraries(metris_deps INTERFACE ${METRIS_DEP_TARGET})
  endif()
  
  # Set include directory variable if requested
  if(METRIS_DEP_INCLUDE_VAR)
    FetchContent_GetProperties(${METRIS_DEP_NAME}_fetch)
    if(${METRIS_DEP_NAME}_fetch_POPULATED)
      set(${METRIS_DEP_INCLUDE_VAR} "${${METRIS_DEP_NAME}_fetch_SOURCE_DIR}" PARENT_SCOPE)
    endif()
  endif()
endfunction()

# Helper for header-only dependencies (like Eigen3)
# Checks for target first, then environment variable, then find_package, then FetchContent
function(metris_find_eigen3)
  # Check if Eigen3 target already exists (from parent project)
  if(TARGET Eigen3::Eigen)
    message(STATUS "Using Eigen3::Eigen target provided by parent project")
    target_link_libraries(metris_deps INTERFACE Eigen3::Eigen)
    metris_register_dependency("find_package" "Eigen3 REQUIRED" "")
    return()
  endif()

  # Check environment variable (backward compatibility)
  if(DEFINED ENV{EIGEN_DIR})
    set(EIGEN3_INCLUDE_DIRS "$ENV{EIGEN_DIR}")
    message(STATUS "Using EIGEN_DIR from environment: ${EIGEN3_INCLUDE_DIRS}")
    if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIRS})
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
      set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
    endif()
    metris_register_dependency("provided" "Eigen3" "${EIGEN3_INCLUDE_DIRS}")
    return()
  endif()

  # Try find_package
  find_package(Eigen3 QUIET)
  if(Eigen3_FOUND)
    message(STATUS "Found Eigen3 via find_package: ${EIGEN3_INCLUDE_DIRS}")
    target_link_libraries(metris_deps INTERFACE Eigen3::Eigen)
    if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>)
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
      set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
    endif()
    metris_register_dependency("find_package" "Eigen3 REQUIRED" "")
    return()
  endif()

  # Fallback to FetchContent
  message(STATUS "Eigen3 not found, fetching from GitLab")
  FetchContent_Declare(
    Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG bcce88c99ed687b756b7a537554cb7c1780b816e
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(Eigen3)
  
  set(EIGEN3_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/_deps/eigen3-src/")
  if(METRIS_INSTALL)
    install(DIRECTORY ${EIGEN3_INCLUDE_DIRS}/Eigen ${EIGEN3_INCLUDE_DIRS}/unsupported
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  
  if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
    set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
  endif()
  
  metris_register_dependency("FetchContent" "Eigen3" "")
endfunction()

# Helper for fmt library
function(metris_find_fmt)
  # Check if fmt target already exists (from parent project)
  if(TARGET fmt::fmt)
    message(STATUS "Using fmt::fmt target provided by parent project")
    target_link_libraries(metris_deps INTERFACE fmt::fmt)
    metris_register_dependency("find_package" "fmt REQUIRED" "")
    return()
  endif()

  # Try find_package first
  find_package(fmt QUIET)
  if(fmt_FOUND)
    message(STATUS "fmt lib found: ${fmt_VERSION}")
    target_link_libraries(metris_deps INTERFACE fmt::fmt)
    metris_register_dependency("find_package" "fmt REQUIRED" "")
    return()
  endif()

  # Fetch via FetchContent
  message(STATUS "fmt lib not found, fetching")
  FetchContent_Declare(
    fmt_fetch
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG e69e5f977d458f2650bb346dadf2ad30c5320281
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(fmt_fetch)

  FetchContent_GetProperties(fmt_fetch)
  if(NOT fmt_fetch_POPULATED)
    message(FATAL_ERROR "fmt was not fetched correctly.")
  endif()

  if(METRIS_INSTALL)
    install(DIRECTORY ${fmt_fetch_SOURCE_DIR}/include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h*")
    set_target_properties(fmt PROPERTIES PUBLIC_HEADER "")
    install(TARGETS fmt
            EXPORT libMetrisTargets
            LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()

  if(NOT TARGET fmt::fmt)
    add_library(fmt::fmt ALIAS fmt)
  endif()

  target_link_libraries(metris_deps INTERFACE fmt::fmt)
  
  if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${fmt_fetch_SOURCE_DIR}/include>)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
    set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
  endif()
  
  metris_register_dependency("FetchContent" "fmt" "")
endfunction()

# Helper for nlohmann/json library
function(metris_find_json)
  # Check if json target already exists (from parent project)
  if(TARGET nlohmann_json::nlohmann_json)
    message(STATUS "Using nlohmann_json::nlohmann_json target provided by parent project")
    target_link_libraries(metris_deps INTERFACE nlohmann_json::nlohmann_json)
    metris_register_dependency("find_package" "nlohmann_json REQUIRED" "")
    return()
  endif()

  # Try find_package first
  find_package(nlohmann_json QUIET)
  if(nlohmann_json_FOUND)
    message(STATUS "nlohmann_json found: ${nlohmann_json_VERSION}")
    target_link_libraries(metris_deps INTERFACE nlohmann_json::nlohmann_json)
    metris_register_dependency("find_package" "nlohmann_json REQUIRED" "")
    return()
  endif()

  # Fetch via FetchContent
  message(STATUS "nlohmann_json not found, fetching")
  FetchContent_Declare(
    json_fetch 
    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
    EXCLUDE_FROM_ALL
  )
  FetchContent_MakeAvailable(json_fetch)
  
  FetchContent_GetProperties(json_fetch)
  if(NOT json_fetch_POPULATED)
    message(FATAL_ERROR "nlohmann_json was not fetched correctly.")
  endif()
  
  # nlohmann/json is header-only, so we need to create an INTERFACE library
  # It doesn't automatically create a target when fetched
  if(NOT TARGET nlohmann_json)
    add_library(nlohmann_json INTERFACE)
    target_include_directories(nlohmann_json INTERFACE
      $<BUILD_INTERFACE:${json_fetch_SOURCE_DIR}/single_include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
    # Ensure include path is available to all consumers
    # Always propagate to metris_deps INTERFACE_INCLUDE_DIRECTORIES
    target_include_directories(metris_deps INTERFACE $<BUILD_INTERFACE:${json_fetch_SOURCE_DIR}/single_include> $<INSTALL_INTERFACE:include/nlohmann>)
  endif()
  
  # Create the namespace alias
  if(NOT TARGET nlohmann_json::nlohmann_json)
    add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
  endif()
  
  target_link_libraries(metris_deps INTERFACE nlohmann_json::nlohmann_json)
  
  if(METRIS_INSTALL)
    # Install headers
    install(FILES ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json.hpp
                  ${json_fetch_SOURCE_DIR}/single_include/nlohmann/json_fwd.hpp
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann/)
    
    # Install the target in the export set (required for metris_deps export)
    # nlohmann_json is an INTERFACE library (header-only)
    install(TARGETS nlohmann_json
            EXPORT libMetrisTargets
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  
  # Include directories are handled via the target's INTERFACE_INCLUDE_DIRECTORIES
  # No need to add to METRIS_EXTERNAL_INCLUDE_DIRS since we're using the target
  
  metris_register_dependency("FetchContent" "json" "")
endfunction()

# Improved NLopt handling for shared dependencies
# This function checks if NLopt is already available from a parent project
function(metris_find_nlopt)
  # First, check if NLopt target already exists (from parent project)
  if(TARGET NLopt::nlopt)
    message(STATUS "Using NLopt::nlopt target provided by parent project")
    target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
    metris_register_dependency("find_package" "NLopt REQUIRED" "")
    return()
  endif()

  # Check if explicitly provided
  if(DEFINED NLOPT_LIBRARIES AND DEFINED NLOPT_INCLUDE_DIRS)
    message(STATUS "Using provided NLOPT_LIBRARIES and NLOPT_INCLUDE_DIRS")
    message(STATUS "NLOPT_LIBRARIES = ${NLOPT_LIBRARIES}")
    message(STATUS "NLOPT_INCLUDE_DIRS = ${NLOPT_INCLUDE_DIRS}")
    target_link_libraries(metris_deps INTERFACE ${NLOPT_LIBRARIES})
    if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
      list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS})
      set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
    endif()
    metris_register_dependency("provided" "NLopt" "${NLOPT_LIBRARIES} ${NLOPT_INCLUDE_DIRS}")
    return()
  endif()

  # Try NLOPT_DIR hint
  if(DEFINED NLOPT_DIR OR DEFINED ENV{NLOPT_DIR})
    if(NOT DEFINED NLOPT_DIR)
      set(NLOPT_DIR $ENV{NLOPT_DIR})
      message(STATUS "Using NLOPT_DIR from environment: ${NLOPT_DIR}")
    endif()
    message(STATUS "Attempting to find external NLopt using hint NLOPT_DIR=${NLOPT_DIR}...")
    find_package(NLopt REQUIRED HINTS "${NLOPT_DIR}")
    
    if(NLopt_FOUND)
      message(STATUS "Found external NLopt: ${NLOPT_LIBRARIES}")
      target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
      metris_register_dependency("find_package" "NLopt REQUIRED HINTS ${NLOPT_DIR}" "")
      return()
    endif()
  endif()

  # Try system-wide find_package
  find_package(NLopt QUIET)
  if(NLopt_FOUND)
    message(STATUS "Found NLopt libraries: ${NLOPT_LIBRARIES}")
    message(STATUS "Found NLopt include directories: ${NLOPT_INCLUDE_DIRS}")
    target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
    metris_register_dependency("find_package" "NLopt REQUIRED" "")
    return()
  endif()

  # Last resort: fetch and build
  message(STATUS "find_package(NLopt) failed, cloning NLopt.")
  FetchContent_Declare(
    nlopt_fetch
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG 019f61ac7253a537760d9cdd9febd927ec97320c
  )
  
  set(NLOPT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(nlopt_fetch)
  
  FetchContent_GetProperties(nlopt_fetch)
  if(NOT nlopt_fetch_POPULATED)
    message(FATAL_ERROR "NLopt was not fetched correctly.")
  endif()
  
  if(NOT TARGET NLopt::nlopt)
    add_library(NLopt::nlopt ALIAS nlopt)
  endif()
  
  target_link_libraries(metris_deps INTERFACE NLopt::nlopt)
  metris_register_dependency("FetchContent" "NLopt" "")
  
  # Set include directories for installation (must be set in parent scope)
  # Note: METRIS_EXTERNAL_INCLUDE_DIRS should be defined in the parent scope
  if(DEFINED METRIS_EXTERNAL_INCLUDE_DIRS)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<BUILD_INTERFACE:${nlopt_fetch_BINARY_DIR}>)
    list(APPEND METRIS_EXTERNAL_INCLUDE_DIRS $<INSTALL_INTERFACE:include>)
    set(METRIS_EXTERNAL_INCLUDE_DIRS ${METRIS_EXTERNAL_INCLUDE_DIRS} PARENT_SCOPE)
  endif()
  
  # Install NLopt if using FetchContent
  if(METRIS_INSTALL)
    install(DIRECTORY ${nlopt_fetch_BINARY_DIR}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            FILES_MATCHING PATTERN "*.h*")
    install(TARGETS nlopt
            EXPORT libMetrisTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endfunction()

