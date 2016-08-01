#!/bin/bash

for f in `ls *.rxyz`
do
	echo -n $f
	
	molecule.graph $f 1.2 | \
		grep -E "(Randic|Wiener|InverseWiener|Balaban|MolecularTopological|Kirchhoff|KirchhoffSum|WienerSum|JOmega)" | \
		awk '{printf "%10.3f", $3}END{print ""}'
		
done > connectivityIndices.dat