#!/bin/bash

iFile=$1

awk '
BEGIN{
	id=1
}
{
	if($0==$1){
		n=$1
		fname=sprintf("mol-%d.xyz",id)
		
		printf( "Generating ... %s ", fname )
		
		print n >> fname
		for(i=1;i<=n+1;i++){
			getline
			print $0 >> fname
		}
		
		print "OK"
	}
	
	id++
}
' $iFile

echo ""
echo ""

declare -A labelsMap

for f in `ls mol-*.xyz`
do
	formula=`molecule.chemicalFormula $f`
	mult=`molecule.minMult $f | awk '{print $1}'`
	labelBase="$formula.q0.m${mult}"
	
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

