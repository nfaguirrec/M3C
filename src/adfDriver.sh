#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2017-2017)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2017-2017)
#          nestor.aguirre@uam.es
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

##
# @brief
##
function runADF()
{
	local iFile=$1
	local nProcShared=$2
	
	export ADFHOME=$M3C_ADF_HOME
	export ADFBIN=$ADFHOME/bin
	export ADFRESOURCES=$ADFHOME/atomicdata
	export SCMLICENSE=$ADFHOME/license.txt
	export SCM_TMPDIR=$M3C_ADF_SCRATCH
	
# 	if [ ! -d  "$SCM_TMPDIR" ]
# 	then
		mkdir -p $SCM_TMPDIR/${iFile%.*}
# 	fi
	
	cp $iFile $SCM_TMPDIR/${iFile%.*}/
	
	pushd . &> /dev/null
	cd $SCM_TMPDIR/${iFile%.*}/
	
	if [ -n "$nProcShared" ]
	then
		$M3C_ADF_HOME/bin/adf -n $nProcShared < $iFile
	else
		$M3C_ADF_HOME/bin/adf < $iFile
	fi
	
	rm -rf $SCM_TMPDIR/${iFile%.*}/
	popd &> /dev/null
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
	local energy=$2
	
	local header=""
	local nAtoms=""
	
	if [ -z "$energy" ]
	then
		header="Geometry from ADF"
	else
		header="Energy = $energy"
	fi
	
	nAtoms=`gawk 'BEGIN{i=0}($0!~/^[[:blank:]]*$/){i++}END{print i}' $iFile`
	
	echo "$nAtoms"
	echo "$header"
	cat $iFile
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
# 	sed -i 's/@MULT/'$mult'/g' .tmp$SID
	sed -i 's/@NSPIN/'$(($mult-1))'/g' .tmp$SID
	
	if [ "$mult" -ne "1" ]
	then
		sed -ri 's/(^[[:blank:]]*Charge .*)/\1\nUNRESTRICTED/gI' .tmp$SID
	fi
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID
}

##
# @brief
##
function optgADFTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.adf
			
		runADF input$SID.adf $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.adf ${xyzFile%.*}.adf 2> /dev/null
		
		if grep "NORMAL TERMINATION" input$SID.out > /dev/null  # @todo No esta claro cual es la linea que asegura la convergencia apropiada
		then
			energy=`grep "<.*Total energy" input$SID.out | tail -n1 | gawk '{print $5}'`
			grep -A$(( $nAtoms+1 )) "Coordinates in Geometry Cycle" input$SID.out | tail -n$nAtoms | sed -r 's/^[[:blank:]]*[[:digit:]]+\.//g' > .finalGeom$SID
			
			geom2xyz .finalGeom$SID $energy
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.adf input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsADFTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	echo "### Error ### freqsADFTemplate is not implemented yet"
	exit
	
# 	local SID="-$xyzFile$RANDOM"
# 	
# 	local nAtoms=""
# 	local energy=""
# 	local fv="" # Vibrational degrees of freedom
# 	
# 	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
# 	
# 	if [ "$nAtoms" -gt 0  ]
# 	then
# 		fillTemplate $template $xyzFile $charge $mult > input$SID.adf
# 		
# 		runADF input$SID.adf $nProcShared > input$SID.out 2>&1
# 		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
# 		cp input$SID.adf ${xyzFile%.*}.adf 2> /dev/null
# 		
# 		if grep "Normal termination" input$SID.out > /dev/null
# 		then
# 			energy=`grep "SCF Done" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
# 			
# 			cat /dev/null > .freqs$SID
# 			
# 			grep "Frequencies" input$SID.out | while read a1 a2 freq1 freq2 freq3
# 			do
# 				if [ -n "$freq1" ]
# 				then
# 					echo $freq1 >> .freqs$SID
# 				fi
# 				
# 				if [ -n "$freq2" ]
# 				then
# 					echo $freq2 >> .freqs$SID
# 				fi
# 				
# 				if [ -n "$freq3" ]
# 				then
# 					echo $freq3 >> .freqs$SID
# 				fi
# 			done
# 			fv=`cat .freqs$SID | wc -l`
# 			
# 			echo $nAtoms
# 			echo "Energy = $energy"
# 			cat $xyzFile | gawk '(NR>2){print $0}'
# 			echo ""
# 			
# 			echo "FREQUENCIES $fv"
# 			cat .freqs$SID
# 			echo ""
# 			
# 			if [ "$nAtoms" -eq 1  ]
# 			then
# 				echo "SYMMETRY SO3"
# 				echo "ELECTRONIC_STATE ??"
# 			else
# 				group=`grep "Full point group" input$SID.out | gawk '{print $4}'`
# 				if [ "$group" = "Nop" -o "$group" = "NOp" ]
# 				then
# 					echo "SYMMETRY ??"
# 				else
# 					echo "SYMMETRY $group"
# 				fi
# 				
# 				state=`grep "The electronic state is" input$SID.out | gawk '{print $5}' | sed 's/\.//'`
# 				if [ -n "$state" ]
# 				then
# 					echo "ELECTRONIC_STATE $state"
# 				else
# 					echo "ELECTRONIC_STATE ??"
# 				fi
# 			fi
# 		else
# 			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
# 		fi
# 	fi
# 	
# 	rm -rf .freqs$SID input$SID.adf input$SID.out
}

##
# @brief
##
function ienerADFTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	echo "### Error ### ienerADFTemplate is not implemented yet"
	exit
	
# 	local SID="-$rxyzFile$RANDOM"
# 	
# 	local nAtoms=""
# 	
# 	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
# 	
# 	if [ "$nAtoms" -gt 0  ]
# 	then
# 		rxyz2xyz $rxyzFile > .xyzFile$SID
# 		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.adf
# 		rm .xyzFile$SID
# 			
# 		runADF input$SID.adf $nProcShared > input$SID.out 2>&1
# 		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
# 		cp input$SID.adf ${rxyzFile%.*}.adf 2> /dev/null
# 		
# 		if grep "Normal termination" input$SID.out > /dev/null
# 		then
# 			energy=`grep -E "^[[:blank:]]+CCSD\(T\)= " input$SID.out | sed 's/D/E/g' | gawk '{printf "%.10f\n", $2}'`
# 			
# 			sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
# 		else
# 			echo "***** FAILURE CONVERGE *****"
# 		fi
# 	fi
# 	
# 	rm -rf .finalGeom$SID input$SID.adf input$SID.out .finalGeom$SID
}
