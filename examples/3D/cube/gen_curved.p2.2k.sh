#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/3D/cube/
OUTFIL=$OUTDIR/curved.p2.2k.meshb

#if [ -f $OUTFIL ]; then
#	echo "File $OUTFIL already exists, skipping."
#	exit 0
#fi

mkdir -p $OUTDIR
mkdir -p tmp/

metris -in cube -cad cube.egads -prefix tmp/ -opt-niter 20 -anamet 2 -out $OUTFIL -t 2