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

awk '
BEGIN{
	loc=0
}
{
	if($0~/^ </) loc=0
	if(loc==1 &&$1~/^[[:digit:]]+./){
		gsub("^[[:digit:]]+.","",$1)
		print $0
	}
	if($0~/>>>> FRAGM/) loc=1
}
' $iFile > .geom
nAtoms=`cat .geom | wc -l`

if [ "$nAtoms" -eq 0  ]
then
	awk '
	BEGIN{
		angs=1.88972612456506
		loc=0
	}
	{
		if($0~/^[[:blank:]]*$/) loc=0
		if(loc==1 && $1~/^[[:digit:]]+$/){
			printf("%5s%15.8f%15.8f%15.8f\n",$2,$3/angs,$4/angs,$5/angs)
		}
		if($0~/Geometry/) loc=1
	}' $iFile > .geom
	nAtoms=`cat .geom | wc -l`
fi

if [ "$nAtoms" -gt 0  ]
then
	energy=`grep "<.*Bond Energy.*a.u." $iFile | tail -n1 | gawk '{print $5}'`  # <<< When only a frequencies calculation was carried out.
	
	[ -z "$energy" ] && energy=`grep "<.*current energy" $iFile | tail -n1 | gawk '{print $5}'` # <<< I think this is just for geometry optimization
	[ -z "$energy" ] && energy=`grep "Total Energy (hartree)" $iFile | tail -n1 | awk '{print $NF}'` # <<< This is for DFTB
	
	awk '
		BEGIN{ loc=0 }
		(loc==2&&$0~/^[[:blank:]]*$/){ loc=0 }
		(loc==2){ print $1 }
		($0~/List of All Frequencies:/){loc=1}
		(loc==1&&$1=="----------"){loc=2}
	' $iFile > .freqs
	fv=`cat .freqs | wc -l`
	
	if [ "$fv" -eq 0 ]
	then
		awk '
			($0~/Index:.*Frequency \(cm-1\)/){print $5}
		' $iFile > .freqs
		fv=`cat .freqs | wc -l`
	fi
	
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
		group=`grep "Symmetry:" $iFile | tail -n1 | gawk '{print $2}'`
		[ -z "$group" ] && group="NOSYM" # < Just to fix a warning message in dftb. dftb does not use symmetry
		
		group=${SYMMETRY_GROUP_MAP[$group]}
		if [ "$group" = "Nop" -o "$group" = "NOp" ]
		then
			echo "SYMMETRY ??"
		else
			echo "SYMMETRY $group"
		fi
		
		# @todo This is not available in ADF
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

