# This should be run from project root !

# CLI arg and default value.
ARG BASE_IMAGE=ubuntu:22.04
# metris_setup is a base image that just prepares the environment
# This way we can then have a development image (doesn't build)
# and a unit test image (builds and runs tests by default)
FROM ${BASE_IMAGE} AS metris_setup


# Prevent prompts
ENV DEBIAN_FRONTEND=noninteractive

# Nothing comes by default, we need 
# build-essential (gcc, make), cmake and git.
# Each RUN creates and can only affect a layer of the image 
# so cleanup has to be done per RUN. 
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  git \
  libboost-all-dev \
  wget \
  && rm -rf /var/lib/apt/lists/*

# -- Install EGADS
# This assumes called by run_docker.sh which downloads these directories
COPY docker/downloads/OCC /progs/OCC/
ENV CASROOT=/progs/OCC/
COPY docker/downloads/ESP /progs/ESP/

COPY docker/install_EGADS.sh /tmp/install_EGADS.sh
RUN /bin/sh /tmp/install_EGADS.sh


# THESE MUST GO AFTER FROM
ARG INSTALL_NLOPT=false
ARG INSTALL_EIGEN3=false
ARG INSTALL_LAPACK=false
ARG INSTALL_FMT=false

# Move args to env so they are available at runtime
# e.g. in entrypoint.sh for _deps checks
ENV INSTALL_NLOPT=${INSTALL_NLOPT}
ENV INSTALL_EIGEN3=${INSTALL_EIGEN3}
ENV INSTALL_LAPACK=${INSTALL_LAPACK}
ENV INSTALL_FMT=${INSTALL_FMT}

COPY docker/optional_deps_helper.sh /tmp/optional_deps_helper.sh
# set -e: exit on error
RUN set -e; \ 
    OPTIONAL_DEPS=""; \
    # Docker uses sh which doesn't support source (. is equivalent)
    . /tmp/optional_deps_helper.sh; \
    # Defined in optional_deps_helper.sh
    add_if_enabled ${INSTALL_NLOPT} "libnlopt-dev"; \
    add_if_enabled ${INSTALL_EIGEN3} "libeigen3-dev"; \
    add_if_enabled ${INSTALL_LAPACK} "liblapack-dev"; \
    add_if_enabled ${INSTALL_FMT} "libfmt-dev"; \
    if [ -n "$OPTIONAL_DEPS" ]; then \
      echo "Installing optional deps: $OPTIONAL_DEPS"; \
      apt-get update && \
      apt-get install -y $OPTIONAL_DEPS && rm -rf /var/lib/apt/lists/*; \
    else \
      echo "No optional deps to install"; \
    fi; 

# Copy entrypoint (sources ESPenv.sh) as it'll be used for all images
# starting from this
COPY docker/entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["entrypoint.sh"]


#----------------------------------------------------------------
# --- Building Metris
FROM metris_setup AS metris_build

# Set the working directory inside the container.
WORKDIR /Metris

# Copy all your project files from your local machine into the container's /app directory.
COPY . .

# Simplest approach is to use cmake facilities for building within
# directories, rather than creating manually and moving there, 
# given docker constraints.
RUN set -e; \
    . /progs/ESP/ESPenv.sh && \
    cmake -B build/release_clang -S . -DMETRIS_BUILD_CORES=4
RUN set -e; \
    . /progs/ESP/ESPenv.sh && \
    cmake --build build/release_clang --parallel 4

# Defined in stage metris_setup
#ENTRYPOINT ["entrypoint.sh"]
# Default command to run when a container is started from this image.
CMD ctest --test-dir build/release_clang --output-on-failure \
          -V -L unit