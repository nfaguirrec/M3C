#!/bin/bash

ATOMIC_SYMBOL[1]="H" 
ATOMIC_SYMBOL[2]="He"
ATOMIC_SYMBOL[3]="Li"
ATOMIC_SYMBOL[4]="Be"
ATOMIC_SYMBOL[5]="B" 
ATOMIC_SYMBOL[6]="C" 
ATOMIC_SYMBOL[7]="N" 
ATOMIC_SYMBOL[8]="O" 
ATOMIC_SYMBOL[9]="F" 
ATOMIC_SYMBOL[10]="Ne"
ATOMIC_SYMBOL[11]="Na"
ATOMIC_SYMBOL[12]="Mg"
ATOMIC_SYMBOL[13]="Al"
ATOMIC_SYMBOL[14]="Si"
ATOMIC_SYMBOL[15]="P" 
ATOMIC_SYMBOL[16]="S" 
ATOMIC_SYMBOL[17]="Cl"
ATOMIC_SYMBOL[18]="Ar"

##
# @brief
##
function runGAUSSIAN()
{
	local iFile=$1
	
	g09  < $iFile
}

##
# @brief
##
function xyz2geom()
{
	local iFile=$1
	
	gawk 'BEGIN{i=0}( NR>2 && $0!~/^[[:blank:]]*$/ ){ print $0 }' $iFile
}

##
# @brief
##
function geom2xyz()
{
	local iFile=$1
	
	local nAtoms=""
	
	nAtoms=`cat $iFile | wc -l`
	
	echo $nAtoms
	echo "Geometry from GAUSSIAN"
	
	cat $iFile | while read a1 a2 a3 a4 a5 a6
	do
		echo ${ATOMIC_SYMBOL[$a2]} $a4 $a5 $a6
	done
}

##
# @brief
##
function fillTemplate()
{
	local template=$1
	local xyzFile=$2
	local charge=$3
	local mult=$4
	
	local SID="-$xyzFile$RANDOM"
	
	if [ -z "$charge" ]
	then
		charge="0"
	fi
	
	if [ -z "$mult" ]
	then
		mult="1"
	fi
	
	xyz2geom $xyzFile > .geom$SID
	
	gawk '{
		if( $0~/@GEOMETRY/ ){
			while( ( getline line < "'.geom$SID'" ) > 0 ) {
				print line
			}
		}else{
			print $0
		}
	}
	' $template > .tmp$SID
	
	sed -i 's/@CHARGE/'$charge'/g' .tmp$SID
	sed -i 's/@MULT/'$mult'/g' .tmp$SID
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID
}

##
# @brief
##
function optgGAUSSIANTemplate()
{
	local template=$1
	local xyzFile=$2
	local charge=$3
	local mult=$4
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.com
			
		runGAUSSIAN input$SID.com > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		
		if grep "Normal termination" input$SID.out > /dev/null
		then
			grep -A$(( $nAtoms+4 )) "Standard orientation:" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			geom2xyz .finalGeom$SID
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.com input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsGAUSSIANTemplate()
{
	local template=$1
	local xyzFile=$2
	local charge=$3
	local mult=$4
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	local fv="" # Vibrational degrees of freedom
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.com
		
		runGAUSSIAN input$SID.com > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		
		energy=`grep "SCF Done" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
		
		cat /dev/null > .freqs$SID
		
		grep "Frequencies" input$SID.out | while read a1 a2 freq1 freq2 freq3
		do
			if [ -n "$freq1" ]
			then
				echo $freq1 >> .freqs$SID
			fi
			
			if [ -n "$freq2" ]
			then
				echo $freq2 >> .freqs$SID
			fi
			
			if [ -n "$freq3" ]
			then
				echo $freq3 >> .freqs$SID
			fi
		done
		fv=`cat .freqs$SID | wc -l`
		
		echo $nAtoms
		echo "Energy = $energy"
		cat $xyzFile | gawk '(NR>2){print $0}'
		echo ""
		echo "FREQUENCIES   $fv"
		cat .freqs$SID
	fi
	
	rm -rf .freqs$SID input$SID.com input$SID.out
}

##
# @brief
##
function GAMESSTemplate_fixEnergyInRXYZ()
{
	local template=$1
	local rxyzFile=$2
	local charge=$3
	local mult=$4
	
# 	local SID="-$rxyzFile$RANDOM"
# 	
# 	local nAtoms=""
# 	local energy=""
# 	
# 	if [ -z "$charge" ]
# 	then
# 		charge="0"
# 	fi
# 	
# 	if [ -z "$mult" ]
# 	then
# 		mult="1"
# 	fi
# 	
# 	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
# 	
# 	if [ "$nAtoms" -gt 0  ]
# 	then
# 		rm -rf input$SID.dat
# 		
# 		rxyz2xyz $rxyzFile > .xyzFile$SID
# 		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.com
# 		rm .xyzFile$SID
# 		
# 		rungms input$SID.com > input$SID.out 2>&1
# 		
# 		energy=`grep "FINAL .* ENERGY IS" input$SID.out | gawk '{print $5}'`
# 		
# 		sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
# 	fi
# 	
# 	if [ "$debug" == "debug" ]
# 	then
# 		mv input$SID.out ${rxyzFile%.*}.out
# 		rm -rf .finalGeom$SID input$SID.dat input$SID.com input$SID.rst .finalGeom$SID
# 	else
# 		rm -rf .finalGeom$SID input$SID.dat input$SID.com input$SID.out input$SID.rst .finalGeom$SID
# 	fi
}
