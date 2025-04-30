#!/bin/bash

rm -r unif
mkdir unif

for sclmet in 0.1 0.05 0.02 0.01
do
	echo "running metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -o out$sclmet -adapt 10 -qopt-niter 5"
	metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -out out$sclmet -adapt 10 -qopt-niter 5
	echo "running metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -o outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads"
	metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -out outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads
done

rm -r anamet2 
mkdir anamet2
for sclmet in 0.5 0.2 0.1 0.05
do
	echo "running metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -o out$sclmet -adapt 10 -qopt-niter 5"
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -out out$sclmet -adapt 10 -qopt-niter 5
	echo "running metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -o outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads"
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -out outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads
done