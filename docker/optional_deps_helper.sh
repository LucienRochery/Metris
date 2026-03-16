#!/bin/sh

add_if_enabled() {
  if [ "$1" = "true" ]; then
    echo "Adding $2 to optional deps"
    OPTIONAL_DEPS="$OPTIONAL_DEPS $2"
  fi
}

