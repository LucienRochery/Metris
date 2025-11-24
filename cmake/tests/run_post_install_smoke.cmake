# Post-install smoke test for Metris package
# This script configures and builds two small consumer projects against the
# installed Metris package in CMAKE_BINARY_DIR (passed by the test invocation).

if(DEFINED INSTALL_PREFIX)
  set(_install_prefix "${INSTALL_PREFIX}")
else()
  # Fallback to common convention: install prefix inside the top-level build dir
  if(DEFINED CMAKE_BINARY_DIR)
    set(_install_prefix "${CMAKE_BINARY_DIR}/install")
  else()
    message(FATAL_ERROR "INSTALL_PREFIX or CMAKE_BINARY_DIR must be set for post-install smoke test")
  endif()
endif()
message(STATUS "Post-install smoke test using install prefix: ${_install_prefix}")

# Configure the no_deps consumer
execute_process(
  COMMAND ${CMAKE_COMMAND} -S "${CMAKE_CURRENT_LIST_DIR}/no_deps" -B "${CMAKE_CURRENT_BINARY_DIR}/no_deps_build" -DMETRIS_DIR=${_install_prefix} -DCMAKE_PREFIX_PATH=${_install_prefix}
  RESULT_VARIABLE _r1
  OUTPUT_VARIABLE _o1
  ERROR_VARIABLE _e1
)
if(NOT _r1 EQUAL 0)
  message(FATAL_ERROR "Configuring no_deps consumer failed:\n${_o1}\n${_e1}")
endif()

# Configure the with_deps consumer
execute_process(
  COMMAND ${CMAKE_COMMAND} -S "${CMAKE_CURRENT_LIST_DIR}/with_deps" -B "${CMAKE_CURRENT_BINARY_DIR}/with_deps_build" -DMETRIS_DIR=${_install_prefix} -DCMAKE_PREFIX_PATH=${_install_prefix}
  RESULT_VARIABLE _r2
  OUTPUT_VARIABLE _o2
  ERROR_VARIABLE _e2
)
if(NOT _r2 EQUAL 0)
  message(FATAL_ERROR "Configuring with_deps consumer failed:\n${_o2}\n${_e2}")
endif()

message(STATUS "Post-install consumers configured successfully")

# Build both consumers to ensure targets can be linked
execute_process(
  COMMAND ${CMAKE_COMMAND} --build "${CMAKE_CURRENT_BINARY_DIR}/no_deps_build" --target all
  RESULT_VARIABLE _b1
  OUTPUT_VARIABLE _bo1
  ERROR_VARIABLE _be1
)
if(NOT _b1 EQUAL 0)
  message(FATAL_ERROR "Building no_deps consumer failed:\n${_bo1}\n${_be1}")
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} --build "${CMAKE_CURRENT_BINARY_DIR}/with_deps_build" --target all
  RESULT_VARIABLE _b2
  OUTPUT_VARIABLE _bo2
  ERROR_VARIABLE _be2
)
if(NOT _b2 EQUAL 0)
  message(FATAL_ERROR "Building with_deps consumer failed:\n${_bo2}\n${_be2}")
endif()

message(STATUS "Post-install consumers built successfully")
