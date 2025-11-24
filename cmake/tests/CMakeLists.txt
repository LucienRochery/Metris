cmake_minimum_required(VERSION 3.14)
project(MetrisConsumerTests NONE)

include(CTest)

option(METRIS_ENABLE_CONSUMER_TESTS "Enable consumer CMake tests" ON)
if(NOT METRIS_ENABLE_CONSUMER_TESTS)
  return()
endif()

# This test folder will configure two small consumer projects into their own build directories
# and then run CMake configure to exercise find_package(Metris) resolution.

add_custom_target(consumer_tests
  COMMENT "Configure and build consumer test projects"
)

# Helper function to add a consumer
function(add_consumer name src_dir)
  set(out_build_dir "${CMAKE_CURRENT_BINARY_DIR}/${name}_build")
  file(MAKE_DIRECTORY "${out_build_dir}")

  add_custom_command(TARGET consumer_tests PRE_BUILD
    # Configure the consumer against the build tree so it uses the generated
    # MetrisConfig.cmake in the top-level build directory instead of any
    # previously installed package in the install prefix. This avoids
    # consumers reading stale installed files during the build phase.
    COMMAND ${CMAKE_COMMAND} -S "${src_dir}" -B "${out_build_dir}" -DMETRIS_DIR=${CMAKE_BINARY_DIR} -DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}
    COMMENT "Configuring consumer ${name}"
  )

  add_custom_command(TARGET consumer_tests PRE_BUILD
    COMMAND ${CMAKE_COMMAND} --build "${out_build_dir}" --target all
    COMMENT "Building consumer ${name}"
  )

endfunction()

# Add consumers
add_consumer(no_deps "${CMAKE_CURRENT_LIST_DIR}/no_deps")
add_consumer(with_deps "${CMAKE_CURRENT_LIST_DIR}/with_deps")

# Make this target part of test fixtures so it runs during `ctest` if desired
add_test(NAME consumer_tests COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target consumer_tests)
set_tests_properties(consumer_tests PROPERTIES TIMEOUT 300)

# A post-install smoke test that configures and builds the consumers against
# the install prefix. This test is not part of the install process; CI or a
# developer should run it after running `make install` (or the equivalent).
add_test(NAME post_install_smoke
  COMMAND ${CMAKE_COMMAND} -DINSTALL_PREFIX=${CMAKE_BINARY_DIR}/install -P "${CMAKE_CURRENT_LIST_DIR}/run_post_install_smoke.cmake"
)
set_tests_properties(post_install_smoke PROPERTIES TIMEOUT 300)
