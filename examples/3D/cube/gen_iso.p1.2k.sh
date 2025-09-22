#!/bin/bash

if [[ -z "$METRIS_CASES_DIR" ]]; then
    echo "Define METRIS_CASES_DIR environment variable" 1>&2
    exit 1
fi

OUTDIR=$METRIS_CASES_DIR/unit/3D/cube/
mkdir -p $OUTDIR

cp cube.meshb $OUTDIR/iso.p1.2k.meshb