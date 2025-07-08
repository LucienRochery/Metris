#!/bin/bash

if ! [ -x "$(command -v metris)" ]; then
  echo 'Error: metris is not in path.' >&2
  exit 1
else
  echo "Using metris = $(which metris)"
fi


rm -r unif
mkdir unif

for sclmet in 0.5 0.2 0.1
do
	metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -main-in-prefix -out outP1_$sclmet -adapt 10 -qopt-niter 5 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
	metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -main-in-prefix -out outCADP1_$sclmet -adapt 10 -qopt-niter 5 -cad square.egads -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
done

rm -r anamet2 
mkdir anamet2
for sclmet in 0.5 0.2 0.1 0.05
do
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -main-in-prefix -out outP1_$sclmet -adapt 10 -qopt-niter 5 -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -main-in-prefix -out outCADP1_$sclmet -adapt 10 -qopt-niter 5 -cad square.egads -adp-opt-niter 0 -adp-unit-stop 99 -opt-smoo-tol 0.05
done