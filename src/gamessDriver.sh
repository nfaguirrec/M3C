#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) 2012-2016 Departamento de Química                                   #
#                         Universidad Autónoma de Madrid                            #
#                         All rights reserved.                                      #
#                                                                                   #
#                         * Néstor F. Aguirre (2012-2016)                           #
#                           nestor.aguirre@uam.es                                   #
#                         * Sergio Díaz-Tendero (2012-2016)                         #
#                           sergio.diaztendero@uam.es                               #
#                         * M. Paul-Antoine Hervieux (2012-2015)                    #
#                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  #
#                         * Manuel Alcamí (2012-2016)                               #
#                           manuel.alcami@uam.es                                    #
#                         * Fernando Martín (2012-2016)                             #
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

declare -A ATOMIC_NUMBER

ATOMIC_NUMBER["H"]=1
ATOMIC_NUMBER["He"]=2
ATOMIC_NUMBER["Li"]=3
ATOMIC_NUMBER["Be"]=4
ATOMIC_NUMBER["B"]=5
ATOMIC_NUMBER["C"]=6
ATOMIC_NUMBER["N"]=7
ATOMIC_NUMBER["O"]=8
ATOMIC_NUMBER["F"]=9
ATOMIC_NUMBER["Ne"]=10
ATOMIC_NUMBER["Na"]=11
ATOMIC_NUMBER["Mg"]=12
ATOMIC_NUMBER["Al"]=13
ATOMIC_NUMBER["Si"]=14
ATOMIC_NUMBER["P"]=15
ATOMIC_NUMBER["S"]=16
ATOMIC_NUMBER["Cl"]=17
ATOMIC_NUMBER["Ar"]=18

ATOMIC_NUMBER["H"]=1
ATOMIC_NUMBER["HE"]=2
ATOMIC_NUMBER["LI"]=3
ATOMIC_NUMBER["BE"]=4
ATOMIC_NUMBER["B"]=5
ATOMIC_NUMBER["C"]=6
ATOMIC_NUMBER["N"]=7
ATOMIC_NUMBER["O"]=8
ATOMIC_NUMBER["F"]=9
ATOMIC_NUMBER["NE"]=10
ATOMIC_NUMBER["NA"]=11
ATOMIC_NUMBER["MG"]=12
ATOMIC_NUMBER["AL"]=13
ATOMIC_NUMBER["SI"]=14
ATOMIC_NUMBER["P"]=15
ATOMIC_NUMBER["S"]=16
ATOMIC_NUMBER["CL"]=17
ATOMIC_NUMBER["AR"]=18

##
# @brief
##
function runGAMESS()
{
	local iFile=$1
	local nProcShared=$2
	
	$M3C_GAMESS_HOME/rungms $iFile 01 $nProcShared
}

##
# @brief
##
function xyz2geom()
{
	local iFile=$1
	
	local SID="-$iFile$RANDOM"
	
	local label=""
	local id=""
	local nAtoms=""
	local i=""
	
	cp $iFile .tmp$SID
	
	for i in "${!ATOMIC_NUMBER[@]}"
	do
		label=$i
		id=${ATOMIC_NUMBER[$i]}
		
		sed -i 's/^[[:blank:]]*'$label' / '$label'   '$id'./g' .tmp$SID
	done
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' .tmp$SID`
	tail -n $nAtoms .tmp$SID
	
	rm .tmp$SID
}

##
# @brief
##
function geom2xyz()
{
	local iFile=$1
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}($0!~/^[[:blank:]]*$/){i++}END{print i}' $iFile`
	
	echo "$nAtoms"
	echo ""
	gawk '{print $1"   "$3"   "$4"   "$5}' $iFile
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
	local oshell=$5
	
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
	
	if [ "$mult" -eq 1 ]
	then
		sed -i 's/@MULT/'$mult' scftyp=rhf/g' .tmp$SID
	else
		if [ -z "$oshell" ]
		then
			sed -i 's/@MULT/'$mult' scftyp=uhf/g' .tmp$SID
		else
			sed -i 's/@MULT/'$mult' scftyp='$oshell'/g' .tmp$SID
		fi
	fi
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID
}

##
# @brief
##
function optgGAMESSTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	local oshell=$6
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		rm -rf input$SID.dat 2> /dev/null
		fillTemplate $template $xyzFile $charge $mult $oshell > input$SID.inp
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${xyzFile%.*}.inp 2> /dev/null
		
		if grep "EXECUTION OF GAMESS TERMINATED NORMALLY" input$SID.out > /dev/null
		then
			grep -A $(( $nAtoms + 3 )) "***** EQUILIBRIUM GEOMETRY LOCATED *****" input$SID.out \
				| tail -n $nAtoms > .finalGeom$SID
			
			geom2xyz .finalGeom$SID
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.inp input$SID.dat input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsGAMESSTemplate()
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
		rm -rf input$SID.dat
		fillTemplate $template $xyzFile $charge $mult > input$SID.inp
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${xyzFile%.*}.inp 2> /dev/null
		
		energy=`grep "FINAL .* ENERGY IS" input$SID.out | head -n1 | gawk '{print $5}'`
		fv=`molecule.fv $xyzFile`
		
		echo $nAtoms
		echo "Energy = $energy"
		cat $xyzFile | gawk '(NR>2){print $0}'
		echo ""
		grep "       FREQUENCY:" input$SID.out \
			| gawk '
				BEGIN{n=0}{for(i=2;i<=NF;i++) if($i>100.0){ val[n]=$i; n++ }}
				END{
					asort(val)
					fv='`echo $fv`'
					print "FREQUENCIES   ", fv
					for (i = n; i > sqrt( (fv-n)**2 ); i--) print val[i]
				}'
	fi
	
	rm -rf input$SID.dat input$SID.inp input$SID.out input$SID.rst
}

##
# @brief
##
function ienerGAMESSTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$rxyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		rm -rf input$SID.dat
		
		rxyz2xyz $rxyzFile > .xyzFile$SID
		fillTemplate $template .xyzFile$SID $charge $mult rohf > input$SID.inp
		rm .xyzFile$SID
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${rxyzFile%.*}.inp 2> /dev/null
		
		if grep "EXECUTION OF GAMESS TERMINATED NORMALLY" input$SID.out > /dev/null
		then
			energy=`grep "COUPLED-CLUSTER ENERGY" input$SID.out | gawk '{print $5}'`  # RHF
			
			if [ -z "$energy" ]
			then
				energy=`grep "CCSD ENERGY" input$SID.out | gawk '{print $3}'`  # ROHF
			fi
			
			sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
		else
			echo "***** FAILURE CONVERGE *****"
		fi
	fi
	
	rm -rf input$SID.dat input$SID.inp input$SID.out
}
