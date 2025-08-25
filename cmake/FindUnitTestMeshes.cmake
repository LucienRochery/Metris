# Ensure that unit tests reference meshes that have generation 
# scripts, and build up generation targets. 
# Some of these can be simple copies.

# This script parses bunit/* files for mesh names
# It then looks for scripts names gen_<meshname>.sh in 
# examples/ with path matching as follows:
# If bunit/* file contains METRIS_CASES_DIR "unit/somepath/meshname"
# then seek examples/somepath/gen_meshname.sh

# If no such script is found, an error is thrown.
# Otherwise, a target is created for the mesh / script pair.


# -- Get the required meshes by reading bunit/test_ files
set(METRIS_UNIT_TEST_MESHES "")
foreach(test_file ${TEST_SOURCES})
  file(READ ${test_file} CMAKE_FILE_CONTENT)

  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+)[\" ]" MATCHES_NOEXT ${CMAKE_FILE_CONTENT})
  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+\\.mesh[\" ])" MATCHES_PMESH ${CMAKE_FILE_CONTENT})
  string(REGEX MATCHALL "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+\\.meshb)" MATCHES_PMESHB ${CMAKE_FILE_CONTENT})

  #message(STATUS "MATCHES_NOEXT = ${MATCHES_NOEXT}")

  foreach(match IN LISTS MATCHES_NOEXT)
    #message(STATUS "Match raw (1): ${match}")
    string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_/]+)[\" ]" "\\1" meshname "${match}")
    #message(STATUS "meshname0 (1): ${meshname}")
    get_filename_component(fext ${meshname} LAST_EXT)
    if(fext STREQUAL ".mesh" OR fext STREQUAL ".meshb")
      #message(STATUS "Skipping, already has extension: ${meshname}")
      continue()
    endif()
    set(meshname "${meshname}.meshb")
    set(meshname "${METRIS_CASES_DIR}/${meshname}")
    #message(STATUS "meshname2 (1): ${meshname}")
    if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
      list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
    endif()
  endforeach()

  foreach(match IN LISTS MATCHES_PMESH)
    #message(STATUS "Match raw (2): ${match}")
    string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.mesh)" "\\1" meshname "${match}")
    get_filename_component(fext ${meshname} LAST_EXT)
    if(fext STREQUAL ".meshb")
      #message(STATUS "Skipping, already has extension: ${meshname}")
      continue()
    endif()
    #message(STATUS "PMESH case base ${match} fext = ${fext}")
    set(meshname "${METRIS_CASES_DIR}/${meshname}")
    #message(STATUS "meshname (2): ${meshname}")
    if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
      list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
    endif()
    #message(STATUS "Test ${test_file} added mesh from PMESH: ${meshname}")
  endforeach()

  foreach(match IN LISTS MATCHES_PMESHB)
    #message(STATUS "Match raw (3): ${match}")
    string(REGEX REPLACE "METRIS_CASES_DIR +\"([a-zA-Z0-9\\.-_]+\\.meshb)" "\\1" meshname "${match}")
    set(meshname "${METRIS_CASES_DIR}/${meshname}")
    #message(STATUS "meshname (3): ${meshname}")
    if(NOT (${meshname} IN_LIST METRIS_UNIT_TEST_MESHES))
      list(APPEND METRIS_UNIT_TEST_MESHES ${meshname})
    endif()
  endforeach()
endforeach()

# -- Check there's a gen_<meshname>.sh for each <meshname>.meshb
foreach(meshname ${METRIS_UNIT_TEST_MESHES})
  get_filename_component(meshname_base ${meshname} NAME_WLE)
  get_filename_component(meshdir ${meshname} DIRECTORY)
  message(STATUS "Mesh dir: ${meshdir}")
  # From meshdir, extract everything after unit
  # then append that to examples/
  string(REGEX REPLACE "^.*unit/" "" meshdir_relative ${meshdir})
  message(STATUS "Mesh dir relative: ${meshdir_relative}")
  set(gen_script "${CMAKE_CURRENT_SOURCE_DIR}/examples/${meshdir_relative}/gen_${meshname_base}.sh")
  message(STATUS "Unit test mesh: ${meshname} : ${meshname_base}")
  if(NOT EXISTS ${gen_script})
    message(FATAL_ERROR "No generation script found for unit test mesh ${meshname}.\nExpected: ${gen_script}")
  endif()
  
  add_custom_command(OUTPUT ${meshname}
                     DEPENDS metris
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/${meshdir_relative}
                     COMMAND ${CMAKE_COMMAND} -E env 
                             METRIS_CASES_DIR=${METRIS_CASES_DIR}
                             ${gen_script}
                     COMMENT "Generating ${meshname}")
  message(STATUS "New command outputs ${meshname}")
endforeach()

add_custom_target(deploy_cases DEPENDS ${METRIS_UNIT_TEST_MESHES})
