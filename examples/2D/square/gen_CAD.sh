#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/2D/square/
OUTFIL=$OUTDIR/square.egads

if [ -f $OUTFIL ]; then
	echo "File $OUTFIL already exists, skipping."
	exit 0
fi

mkdir -p $OUTDIR
mkdir -p tmp/

cp square.egads $OUTFIL