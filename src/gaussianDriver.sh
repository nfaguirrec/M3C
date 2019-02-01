#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) 2013-2016 Departamento de Química                                   #
#                         Universidad Autónoma de Madrid                            #
#                         All rights reserved.                                      #
#                                                                                   #
#                         * Néstor F. Aguirre (2013-2016)                           #
#                           nestor.aguirre@uam.es                                   #
#                         * Sergio Díaz-Tendero (2013-2016)                         #
#                           sergio.diaztendero@uam.es                               #
#                         * M. Paul-Antoine Hervieux (2013-2015)                    #
#                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  #
#                         * Manuel Alcamí (2013-2016)                               #
#                           manuel.alcami@uam.es                                    #
#                         * Fernando Martín (2013-2016)                             #
#                           fernando.martin@uam.es                                  #
#                                                                                   #
#  Redistribution and use in source and binary forms, with or without               #
#  modification, are permitted provided that the following conditions are met:      #
#                                                                                   #
#  1. Redistributions of source code must retain the above copyright notice, this   #
#     list of conditions and the following disclaimer.                              #
#  2. Redistributions in binary form must reproduce the above copyright notice,     #
#     this list of conditions and the following disclaimer in the documentation     #
#     and/or other materials provided with the distribution.                        #
#  3. Neither the name of the copyright holders nor the names of its contributors   #
#     may be used to endorse or promote products derived from this software         #
#     without specific prior written permission.                                    #
#                                                                                   #
#  The copyright holders provide no reassurances that the source code provided      #
#  does not infringe any patent, copyright, or any other intellectual property      #
#  rights of third parties.  The copyright holders disclaim any liability to any    #
#  recipient for claims brought against recipient by any third party for            #
#  infringement of that parties intellectual property rights.                       #
#                                                                                   #
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  #
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    #
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           #
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR  #
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     #
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      #
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       #
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    #
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     #
#                                                                                   #
#####################################################################################

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
	local nProcShared=$2
	
	export GAUSS_EXEDIR=$M3C_GAUSSIAN_HOME
	export GAUSS_SCRDIR=$M3C_GAUSSIAN_SCRATCH
	if [ ! -d  "$GAUSS_SCRDIR" ]
	then
		mkdir -p $GAUSS_SCRDIR
	fi
	
	if [ -n "$nProcShared" ]
	then
		awk '
		BEGIN{
			line=1
		}
		($1!~/NProcShared/){
			map[line]=$0
			line++
		}
		END{
			print "%NProcShared='$nProcShared'"
			for(i=1;i<line;i+=1)
				print map[i]
		}
		' $iFile > .$$tmp-$iFile
		mv .$$tmp-$iFile $iFile
	fi
	
	$M3C_GAUSSIAN_HOME/g09  < $iFile
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
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.com
			
		runGAUSSIAN input$SID.com $nProcShared > input$SID.out 2>&1
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
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	local fv="" # Vibrational degrees of freedom
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.com
		
		runGAUSSIAN input$SID.com $nProcShared > input$SID.out 2>&1
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
		
		echo "FREQUENCIES $fv"
		cat .freqs$SID
		echo ""
		
		if [ "$nAtoms" -eq 1  ]
		then
			echo "SYMMETRY SO3"
			echo "ELECTRONIC_STATE ??"
		else
			#group=`grep "Full point group" input$SID.out | gawk '{print $4}'`
			group=`grep "Framework group" input$SID.out | sed 's/\[/ /g' | gawk '{print $3}'`
			if [ "$group" = "Nop" -o "$group" = "NOp" ]
			then
				echo "SYMMETRY ??"
			else
				echo "SYMMETRY $group"
			fi
			
			state=`grep "The electronic state is" input$SID.out | gawk '{print $5}' | sed 's/\.//'`
			if [ -n "$state" ]
			then
				echo "ELECTRONIC_STATE $state"
			else
				echo "ELECTRONIC_STATE ??"
			fi
		fi
	fi
	
	rm -rf .freqs$SID input$SID.com input$SID.out
}

##
# @brief
##
function ienerGAUSSIANTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$rxyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		rxyz2xyz $rxyzFile > .xyzFile$SID
		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.com
		rm .xyzFile$SID
			
		runGAUSSIAN input$SID.com $nProcShared > input$SID.out 2>&1
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
