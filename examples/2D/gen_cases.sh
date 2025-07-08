#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

DIR0=$(pwd)
cd square
for scr in ./gen*sh
do
	echo "Running $scr"
	./$scr
done
cd $DIR0

cp -r misc $METRIS_CASES_DIR/unit/2D/
