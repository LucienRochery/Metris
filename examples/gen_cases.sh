#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

DIR0=$(pwd)

cd $METRIS_CASES_DIR

# For directories that don't contain a gen*sh
shopt -s nullglob

for dir in */
do
	cd $dir
	for scr in ./gen*sh
	do
		echo "Running $scr"
		./$scr
	done
	cd $METRIS_CASES_DIR
done

cd $DIR0

shopt -u nullglob