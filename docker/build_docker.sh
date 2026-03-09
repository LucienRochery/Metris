#!/bin/bash

# Run this from within Metris/docker/
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
if [ "$PWD" != "$SCRIPT_DIR" ]; then
  echo "Script must be run from within <METRIS_ROOT>/docker"
  echo "SCRIPT_DIR = $SCRIPT_DIR"
  exit 1
fi

ARCH=$(uname -m)
[[ "$ARCH" == "arm64" ]] && ARCH="aarch64"
OCC_TARBALL="OCC781-linux-${ARCH}.tgz"

cd downloads/
if [ ! -f "$OCC_TARBALL" ]; then
  wget "https://acdl.mit.edu/ESP/${OCC_TARBALL}"
fi
if [ ! -d "OCC" ]; then
  mkdir OCC
  tar --strip-components=1 -xf "$OCC_TARBALL" -C OCC/
fi

if [ ! -d "ESP" ] ; then
  wget https://acdl.mit.edu/ESP/ESP.tgz
  mkdir ESP
  tar --strip-components=1 -xvf ESP.tgz -C ESP/
fi
cd ..
# Back in docker/


# Run dockerfile from root directory
cd ..
echo ""
echo "-----------------------------------------------"
echo "-- Building Metris with ALL systemwide libraries"
echo "-----------------------------------------------"
docker build --build-arg INSTALL_NLOPT=true \
             --build-arg INSTALL_EIGEN3=true \
             --build-arg INSTALL_LAPACK=true \
             --build-arg INSTALL_FMT=true \
             --target metris_setup \
             -t metris_alldeps_dev -f docker/Dockerfile .
# -v mounts to the Docker filesystem: this allows changing files 
# without having to rebuild because of COPY
docker run -it --rm -v "$(pwd)":/Metris metris_alldeps_dev /bin/bash

docker build --build-arg INSTALL_NLOPT=true \
             --build-arg INSTALL_EIGEN3=true \
             --build-arg INSTALL_LAPACK=true \
             --build-arg INSTALL_FMT=true \
             -t metris_alldeps -f docker/Dockerfile .


echo ""
echo "-----------------------------------------------"
echo "-- Building Metris with NO systemwide libraries"
echo "-----------------------------------------------"
docker build --build-arg INSTALL_NLOPT=false \
             --build-arg INSTALL_EIGEN3=false \
             --build-arg INSTALL_LAPACK=false \
             --build-arg INSTALL_FMT=false \
             -t metris_nodeps -f docker/Dockerfile .



# --rm cleans up after
docker run --rm metris_nodeps
docker run --rm metris_alldeps