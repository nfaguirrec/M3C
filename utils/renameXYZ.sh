#!/bin/bash

charge=$1

if [ -z "$charge" ]
then
	charge="0"
fi

declare -A labelsMap

for f in `ls *.xyz`
do
	formula=`molecule.chemicalFormula $f`
	mult=`molecule.minMult $f $charge | awk '{print $1}'`
	labelBase="$formula.q${charge}.m${mult}"
	
	if [ -z "${labelsMap[$labelBase]}" ]
	then
		labelsMap[$labelBase]="1"
	fi
	
	label="$labelBase-${labelsMap[$labelBase]}"
	
	echo -n "Moving $f to $label.xyz ... "
	mv $f $label.xyz
	echo "OK"
	
	labelsMap[$labelBase]=$(( ${labelsMap[$labelBase]} + 1 ))
done

