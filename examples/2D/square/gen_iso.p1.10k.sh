#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/2D/square/
mkdir -p $OUTDIR
OUTFIL=$OUTDIR/iso.p1.10k.meshb

#if [ -f $OUTFIL ]; then
#	echo "File $OUTFIL already exists, skipping."
#	exit 0
#fi

mkdir -p $OUTDIR
mkdir -p tmp/
metris -in square -cad square -anamet 1 -sclmet 0.21 -prefix tmp/ -out $OUTFIL -adapt 20