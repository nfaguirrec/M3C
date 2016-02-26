#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2012-2016)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2012-2016)
#          nestor.aguirre@uam.es
#    * Dr. Sergio Díaz-Tendero (2012-2015)
#          sergio.diaztendero@uam.es
#    * Prof. M. Paul-Antoine Hervieux (2012-2015)
#          Paul-Antoine.Hervieux@ipcms.unistra.fr
#    * Prof. Fernando Martín (2012-2015)
#          fernando.martin@uam.es
#    * Prof. Manuel Alcamí (2012-2015)
#          manuel.alcami@uam.es
#  
#  Redistribution and use in source and binary forms, with or
#  without modification, are permitted provided that the
#  following conditions are met:
#  
#   * Redistributions of binary or source code must retain
#     the above copyright notice and this list of conditions
#     and/or other materials provided with the distribution.
#   * All advertising materials mentioning features or use of
#     this software must display the following acknowledgement:
#     
#     This product includes software from M3C project.
#
##################################################################

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
		printf "%5s%12.6f%12.6f%12.6f\n" ${ATOMIC_SYMBOL[$a2]} $a4 $a5 $a6
	done
}

##
# @brief
##
function rxyz2xyz()
{
	local iFile=$1
	
	gawk 'BEGIN{ n=1 }(NR==1){nAtoms=$1}(NR==2){title=$0; print nAtoms; print title;}( NR>2 && n<=nAtoms ){ print $0; n++ }' $iFile
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
		cp input$SID.com ${xyzFile%.*}.com 2> /dev/null
		
		if grep "Normal termination" input$SID.out > /dev/null
		then
# 			grep -A$(( $nAtoms+4 )) "Standard orientation:" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			grep -A$(( $nAtoms+4 )) "Input orientation:" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			geom2xyz .finalGeom$SID
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
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
		cp input$SID.com ${xyzFile%.*}.com 2> /dev/null
		
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
function ienerGAUSSIANTemplate()
{
	local template=$1
	local rxyzFile=$2
	local charge=$3
	local mult=$4
	
	local SID="-$rxyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		rxyz2xyz $rxyzFile > .xyzFile$SID
		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.com
		rm .xyzFile$SID
			
		runGAUSSIAN input$SID.com > input$SID.out 2>&1
		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
		cp input$SID.com ${rxyzFile%.*}.com 2> /dev/null
		
		if grep "Normal termination" input$SID.out > /dev/null
		then
			energy=`grep -E "^[[:blank:]]+CCSD\(T\)= " input$SID.out | sed 's/D/E/g' | gawk '{printf "%.10f\n", $2}'`
			
			sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
		else
			echo "***** FAILURE CONVERGE *****"
		fi
	fi
	
	rm -rf .finalGeom$SID input$SID.com input$SID.out .finalGeom$SID
}
