#!/bin/bash

# Nightly script: run regression tests for merging onto main

cd $WORKSPACE

export METRIS_CASES_DIR=$WORKSPACE/examples/

#Files might linger if a build was aborted
echo "Removing any lingering untracked files"
for f in `git ls-files --others --exclude=build --exclude=external --exclude=bunit/regression`; do
  echo "rm -rf $f"
done

cmakedir=$WORKSPACE/build/$builddir

#Create the build directory if it does not exist
mkdir -p $cmakedir
cd $cmakedir

source $WORKSPACE/jenkins/jenkins_env.sh

# Number of processors used in compile
nproc=1
CMAKEARGS="-DMETRIS_BUILD_CORES=$nproc"


time source $WORKSPACE/jenkins/cmake_jenkins.sh

# Copy over the makefile that pipes parallel execution to files
cp $WORKSPACE/jenkins/Makefile.parallel .


#Build basic Metris targets
#time make -j $nproc -f Makefile.parallel metris
#time make -j $nproc -f Makefile.parallel unit_build
#time make -f Makefile.parallel install

# This should handle all the necessary compilations
ctest --output-on-failure -V -L regression 2>&1 |tee regression.log

#Fail the build if any files were generated in source.
#Count the number of files
cd $WORKSPACE
tmpfiles=`git ls-files --others --exclude=build --exclude=external --exclude=bunit/regression | wc -l`

if [ $tmpfiles -ne 0 ]; then
  set -x
  echo "error: Files should not be generated in code pushed to the repository."
  echo "error: Files found :"
  for f in `git ls-files --others --exclude=build --exclude=external --exclude=bunit/regression`; do
    echo "error: $f"
  done
  exit 1
fi
