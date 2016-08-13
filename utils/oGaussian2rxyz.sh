iFile=$1
source $M3C_HOME/bin/gaussianDriver.sh

nAtoms=`grep "NAtoms=" $iFile | awk '{print $2; exit}'`
( grep -A$(( $nAtoms+4 )) "Input orientation:" $iFile | tail -n$nAtoms | \
while read a1 a2 a3 a4 a5 a6
do
	printf "%5s%12.6f%12.6f%12.6f\n" ${ATOMIC_SYMBOL[$a2]} $a4 $a5 $a6
done ) > .geom

if [ "$nAtoms" -gt 0  ]
then
	energy=`grep "SCF Done" $iFile | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
	
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
		echo "SYMMETRY SO3"
		echo "ELECTRONIC_STATE ??"
	else
		group=`grep "Full point group" $iFile | gawk '{print $4}'`
		if [ "$group" = "Nop" -o "$group" = "NOp" ]
		then
			echo "SYMMETRY ??"
		else
			echo "SYMMETRY $group"
		fi
		
		state=`grep "The electronic state is" $iFile | gawk '{print $5}' | sed 's/\.//'`
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
