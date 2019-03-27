#!/bin/bash

# for f in `ls *.out`; do echo ${f%*.com.out}; mv $f ${f%*.com.out}.out; done

# for f in `ls *.out`; do echo $f; oGaussian2rxyz.sh $f > ${f%*.out}.rxyz; done

for f in `ls *.out`
do
	echo -n "$f --> "
	oGaussian2rxyz.sh $f > mol.xyz
	formula=`molecule.chemicalFormula mol.xyz`
	
	charge=`grep Multiplicity $f | tail -n1 | awk '{print $3}'`
	mult=`grep Multiplicity $f | tail -n1 | awk '{print $6}'`
	
	id=1
	while true
	do
		if [ ! -f ../"${formula}.q${charge}.m${mult}-$id.out" ]
		then
			break
		fi
		
		id=$(( $id+1 ))
	done
	
	echo "${formula}.q${charge}.m${mult}-$id.out"
	
	cp ${f%.*}.com ../"${formula}.q${charge}.m${mult}-$id.com"
	cp $f ../"${formula}.q${charge}.m${mult}-$id.out"
	
done

rm mol.xyz
