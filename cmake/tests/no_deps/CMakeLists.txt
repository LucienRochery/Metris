cmake_minimum_required(VERSION 3.14)
project(metris_consumer_no_deps)

# Consumer that only calls find_package(Metris)
find_package(Metris CONFIG REQUIRED)

add_executable(consumer_no_deps main.cpp)
target_link_libraries(consumer_no_deps PRIVATE Metris::libMetris)

message(STATUS "Consumer no_deps configured. Metris include dirs: ${METRIS_INCLUDE_DIRS}")
