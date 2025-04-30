#!/bin/bash

rm -r unifP2
mkdir unifP2

for sclmet in 0.1 0.05 0.02
do
	echo "running metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -o out$sclmet -adapt 10 -qopt-niter 5"
	metris -in square -anamet 1 -sclmet $sclmet -prefix unifP2/ -out out$sclmet -adapt 10 -qopt-niter 5 -t 2
	echo "running metris -in square -anamet 1 -sclmet $sclmet -prefix unif/ -o outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads"
	metris -in square -anamet 1 -sclmet $sclmet -prefix unifP2/ -out outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads  -t 2
done

rm -r anamet2P2
mkdir anamet2P2
for sclmet in 0.5 0.2 0.1
do
	echo "running metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -o out$sclmet -adapt 10 -qopt-niter 5"
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2P2/ -out out$sclmet -adapt 10 -qopt-niter 5 -t 2
	echo "running metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2/ -o outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads"
	metris -in square -anamet 2 -sclmet $sclmet -prefix anamet2P2/ -out outCAD$sclmet -adapt 10 -qopt-niter 5 -cad square.egads -t 2
done