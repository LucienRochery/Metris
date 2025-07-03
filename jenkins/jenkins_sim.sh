#!/bin/bash

if [ $# -eq 1 ]; then
  echo "Using workspace $1"
  JENKINS_TEST_WORKSPACE=$1
fi

if [ -z ${JENKINS_TEST_WORKSPACE:+x} ]; then
  echo "Workspace unset. Define JENKINS_TEST_WORKSPACE or pass in to script"
  exit 1
fi

WORKSPACE="${JENKINS_TEST_WORKSPACE}/Metris"

CURDIR=$(pwd)

cd $JENKINS_TEST_WORKSPACE


if [ -d "Metris" ]; then
  echo "$WORKSPACE exists, pulling."
  cd Metris
  git pull
else
  echo "$WORKSPACE does not exist, cloning."
  git clone git@github.com:LucienRochery/Metris.git
  if [ ! -d "Metris" ]; then
    echo "Directory Metris was not created in $JENKINS_TEST_WORKSPACE."
    exit 1
  fi
fi

cd Metris
git checkout LucienRochery/develop
cd ..

cp -r $METRIS_DIR/jenkins $WORKSPACE/

read  -n 1 -p "About to run Metris_Commit.sh" dummyvariable

cd $WORKSPACE/jenkins
builddir="release_clang"
source Metris_Commit.sh


cd $CURDIR