#!/bin/bash

if ! [ -x "$(command -v metris)" ]; then
  echo 'Error: metris is not in path.' >&2
  exit 1
else
  echo 'Using metris = $(which metris)'
fi


rm -r unifP2
mkdir unifP2

for sclmet in 0.1 0.05 0.02
do
	metris -in square -anamet 1 -sclmet $sclmet -prefix unifP2/ -main-in-prefix -out outP2_$sclmet -adapt 10 -qopt-niter 5 -t 2 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
	metris -in square -anamet 1 -sclmet $sclmet -prefix unifP2/ -main-in-prefix -out outCADP2_$sclmet -adapt 10 -qopt-niter 5 -cad square.egads  -t 2 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
done

rm -r anamet2P2
mkdir anamet2P2
for sclmet in 0.5 0.2 0.1
do
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2P2/ -main-in-prefix -out outP2_$sclmet -adapt 10 -qopt-niter 5 -t 2 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2P2/ -main-in-prefix -out outCADP2_$sclmet -adapt 10 -qopt-niter 5 -cad square.egads -t 2 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
done