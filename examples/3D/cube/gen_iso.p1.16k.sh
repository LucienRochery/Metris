#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/3D/cube/
OUTFIL=$OUTDIR/iso.p1.16k.meshb

#if [ -f $OUTFIL ]; then
#	echo "File $OUTFIL already exists, skipping."
#	exit 0
#fi

mkdir -p $OUTDIR
mkdir -p tmp/
metris -in cube -cad cube -anamet 1 -sclmet 1 -prefix tmp/ -out $OUTFIL -adapt 10 