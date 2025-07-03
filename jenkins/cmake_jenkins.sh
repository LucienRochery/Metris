#!/bin/bash

WORKSPACE=$(git rev-parse --show-toplevel)

echo "Running CMake in directory $(pwd) with CMakeLists.txt in $WORKSPACE"

cmake $CMAKEARGS \
      $@ \
      $WORKSPACE
