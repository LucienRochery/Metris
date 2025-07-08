#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/2D/square/
OUTFIL=$OUTDIR/circmet.p2.500.meshb

if [ -f $OUTFIL ]; then
	echo "File $OUTFIL already exists, skipping."
	exit 0
fi

mkdir -p $OUTDIR
mkdir -p tmp/

rmetris -in square -anamet 2 -sclmet 0.5 -out $OUTFIL -prefix tmp/ -adapt 10 -opt-niter 20 -t 2 -adp-opt-niter -1