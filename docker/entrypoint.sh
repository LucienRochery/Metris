#!/bin/sh
set -e

echo "Sourcing EGADS_env.sh at /progs/ESP/ESPenv.sh"
. /progs/ESP/ESPenv.sh


#verify_dep_build() {
#  local system_install="$1"
#  local dep_name="$2"
#  local deps_path="/Metris/build/release_clang/_deps/${3}"
#
#  echo "Verifying configuration for: $dep_name"
#
#  if [ "$system_install" = "true" ]; then
#    # If we installed it with apt, the _deps directory should NOT exist.
#    if [ -d "$deps_path" ]; then
#      echo "Error: $dep_name was built from source but should have been found on the system." >&2
#      exit 1
#    else
#      echo "OK: $dep_name was correctly found on the system."
#    fi
#  else
#    # If we did NOT install it with apt, the _deps directory SHOULD exist.
#    if [ ! -d "$deps_path" ]; then
#      echo "Error: $dep_name was not built from source but should have been." >&2
#      exit 1
#    else
#      echo "OK: $dep_name was correctly built from source."
#    fi
#  fi
#}

# Docker concatenates ENTRYPOINT and CMD
# With exec, we are basically chaining this script into calling the first
# argument of CMD "as if sourced" (from the same shell)
# "$@" just pipes the arguments as-is to that command.
exec "$@"