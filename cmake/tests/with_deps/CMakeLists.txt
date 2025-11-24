cmake_minimum_required(VERSION 3.14)
project(metris_consumer_with_deps)

# Consumer that already finds NLopt and Eigen before finding Metris
find_package(Eigen3 REQUIRED)
find_package(NLopt REQUIRED)

find_package(Metris CONFIG REQUIRED)

add_executable(consumer_with_deps main.cpp)
# Link to Metris library and to NLopt/Eigen to ensure co-existence
target_link_libraries(consumer_with_deps PRIVATE Metris::libMetris NLopt::nlopt Eigen3::Eigen)

message(STATUS "Consumer with_deps configured. Eigen include dirs: ${EIGEN3_INCLUDE_DIRS}")
