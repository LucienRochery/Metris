#!/bin/bash

shopt -s nocasematch

dir=`pwd`
builddir=`basename $dir`

CAS_VERSION=7.8.1
export CASREV=$CAS_VERSION
if [[ `hostname` == *"reynolds"* ]]; then

  if [[ $builddir == *"intel"* || $builddir == *"coverage"* ]]; then
    source /opt/intel/oneapi/setvars.sh > intelvars.out 2>&1
  fi
  # This makes sure cmake finds open-mpi rather than intel mpi
  export PATH=/usr/bin:$PATH
  export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH

  export LAPACK_DIR=/home/jenkins/util/lapack/lapack-3.6.1/install
  export ESP_DIR=/home/jenkins/util/ESP/EngSketchPad
  export CAS_DIR=/home/jenkins/util/ESP/OpenCASCADE-$CAS_VERSION
  if [ -d /home/jenkins/util/NLOPT ]; then
    export NLOPT_DIR=`ls -d /home/jenkins/util/NLOPT/nlopt-*/install`
  fi
  if [ -d /home/jenkins/util/ccache ]; then
    export PATH=`ls -d /home/jenkins/util/ccache/ccache-*/install/bin`:$PATH
  fi

elif [[ `hostname` == *"macys"* ]]; then

  export PATH=/usr/local/bin:/usr/local/sbin:/usr/local/opt/python/libexec/bin:$PATH
  export LAPACK_DIR=/usr/local/opt/lapack/
  export ESP_DIR=/Users/jenkins/util/ESP/EngSketchPad
  export CAS_DIR=/Users/jenkins/util/ESP/OpenCASCADE-$CAS_VERSION
  export PATH=/Users/jenkins/util/fefloa/bin:$PATH

elif [[ `hostname` == *"viggen"* ]]; then

  export PATH=/opt/homebrew/bin:/opt/homebrew/sbin:/opt/homebrew/opt/python/libexec/bin:$PATH
  export LAPACK_DIR=/opt/homebrew/opt/lapack/
  export ESP_DIR=/Users/jenkins/util/ESP/EngSketchPad
  export CAS_DIR=/Users/jenkins/util/ESP/OpenCASCADE-$CAS_VERSION
else
  if [ -z "$LAPACK_DIR" ]; then
    echo "Please set LAPACK_DIR in your environment."
    exit 0
  fi
fi

#Use the system version by default
export Boost_NO_SYSTEM_PATHS=OFF

export ASAN_OPTIONS=detect_odr_violation=0

env > env.log 2>&1
