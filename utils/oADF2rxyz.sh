#!/bin/bash

iFile=$1

if [ -f "$M3C_HOME/bin/adfDriver.sh" ]
then
	source $M3C_HOME/bin/adfDriver.sh
elif [ -f "$M3C_HOME/src/adfDriver.sh" ]
then
	source $M3C_HOME/src/adfDriver.sh
else
	echo "### ERROR ### adfDriver.sh: No such file"
	echo "              Check \$M3C_HOME variable"
	exit
fi

nAtoms=`grep "NAtoms=" $iFile | awk '{print $2; exit}'`
( grep -A$(( $nAtoms+4 )) "Input orientation:" $iFile | tail -n$nAtoms | \
while read a1 a2 a3 a4 a5 a6
do
	printf "%5s%12.6f%12.6f%12.6f\n" ${ATOMIC_SYMBOL[$a2]} $a4 $a5 $a6
done ) > .geom

if [ "$nAtoms" -gt 0  ]
then
	energy=`grep -E "^[[:blank:]]+CCSD\(T\)= " $iFile | sed 's/D/E/g' | gawk '{printf "%.10f\n", $2}'`
	
	if [ -z "$energy" ]
	then
		energy=`grep "SCF Done" $iFile | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
	fi
	
	cat /dev/null > .freqs
	
	grep "Frequencies" $iFile | while read a1 a2 freq1 freq2 freq3
	do
		if [ -n "$freq1" ]
		then
			echo $freq1 >> .freqs
		fi
		
		if [ -n "$freq2" ]
		then
			echo $freq2 >> .freqs
		fi
		
		if [ -n "$freq3" ]
		then
			echo $freq3 >> .freqs
		fi
	done
	fv=`cat .freqs | wc -l`
	
	echo $nAtoms
	echo "Energy = $energy"
	cat .geom
	echo ""
	
	echo "FREQUENCIES $fv"
	cat .freqs
	echo ""
	
	if [ "$nAtoms" -eq 1  ]
	then
		echo "SYMMETRY R3"
		echo "ELECTRONIC_STATE ??"
	else
		group=`grep "Full point group" $iFile | tail -n1 | gawk '{print $4}'`
		if [ "$group" = "Nop" -o "$group" = "NOp" ]
		then
			echo "SYMMETRY ??"
		else
			echo "SYMMETRY $group"
		fi
		
		state=`grep "The electronic state is" $iFile | tail -n1 | gawk '{print $5}' | sed 's/\.//'`
		if [ -n "$state" ]
		then
			echo "ELECTRONIC_STATE $state"
		else
			echo "ELECTRONIC_STATE ??"
		fi
	fi
fi

rm .geom
rm .freqs

