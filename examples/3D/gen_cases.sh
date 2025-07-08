#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

DIR0=$(pwd)
cd cube
for scr in ./gen*sh
do
	echo "Running $scr"
	./$scr
done
cd $DIR0
