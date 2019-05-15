#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) by authors                                                          #
#                                                                                   #
#  Authors:                                                                         #
#                         * Néstor F. Aguirre (2017-2017)                           #
#                           nfaguirrec@gmail.com                                    #
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

declare -A SYMMETRY_GROUP_MAP

SYMMETRY_GROUP_MAP["ATOM"]="R3"
SYMMETRY_GROUP_MAP["C1"]="C1"
SYMMETRY_GROUP_MAP["Cs"]="CS"
SYMMETRY_GROUP_MAP["Ci"]="CI"
SYMMETRY_GROUP_MAP["C2"]="C2"
SYMMETRY_GROUP_MAP["C3"]="C3"
SYMMETRY_GROUP_MAP["C4"]="C4"
SYMMETRY_GROUP_MAP["C5"]="C5"
SYMMETRY_GROUP_MAP["C6"]="C6"
SYMMETRY_GROUP_MAP["C7"]="C7"
SYMMETRY_GROUP_MAP["C8"]="C8"
SYMMETRY_GROUP_MAP["D2"]="D2"
SYMMETRY_GROUP_MAP["D3"]="D3"
SYMMETRY_GROUP_MAP["D4"]="D4"
SYMMETRY_GROUP_MAP["D5"]="D5"
SYMMETRY_GROUP_MAP["D6"]="D6"
SYMMETRY_GROUP_MAP["C2v"]="C2V"
SYMMETRY_GROUP_MAP["C3v"]="C3V"
SYMMETRY_GROUP_MAP["C4v"]="C4V"
SYMMETRY_GROUP_MAP["C5v"]="C5V"
SYMMETRY_GROUP_MAP["C6v"]="C6V"
SYMMETRY_GROUP_MAP["C2h"]="C2H"
SYMMETRY_GROUP_MAP["C3h"]="C3H"
SYMMETRY_GROUP_MAP["C4h"]="C4H"
SYMMETRY_GROUP_MAP["C5h"]="C5H"
SYMMETRY_GROUP_MAP["C6h"]="C6H"
SYMMETRY_GROUP_MAP["D2h"]="D2H"
SYMMETRY_GROUP_MAP["D3h"]="D3H"
SYMMETRY_GROUP_MAP["D4h"]="D4H"
SYMMETRY_GROUP_MAP["D5h"]="D5H"
SYMMETRY_GROUP_MAP["D6h"]="D6H"
SYMMETRY_GROUP_MAP["D8h"]="D8H"
SYMMETRY_GROUP_MAP["D2d"]="D2D"
SYMMETRY_GROUP_MAP["D3d"]="D3D"
SYMMETRY_GROUP_MAP["D4d"]="D4D"
SYMMETRY_GROUP_MAP["D5d"]="D5D"
SYMMETRY_GROUP_MAP["D6d"]="D6D"
SYMMETRY_GROUP_MAP["S4"]="S4"
SYMMETRY_GROUP_MAP["S6"]="S6"
SYMMETRY_GROUP_MAP["S8"]="S8"
SYMMETRY_GROUP_MAP["T"]="T"
SYMMETRY_GROUP_MAP["Th"]="TH"
SYMMETRY_GROUP_MAP["Td"]="TD"
SYMMETRY_GROUP_MAP["O"]="O"
SYMMETRY_GROUP_MAP["Oh"]="OH"
SYMMETRY_GROUP_MAP["I"]="I"
SYMMETRY_GROUP_MAP["Ih"]="IH"

# Nwchem 6.6 do not use non-Abelian groups
# SYMMETRY_GROUP_MAP["C(LIN)"]="C*V"
# SYMMETRY_GROUP_MAP["D(LIN)"]="D*H"

##
# @brief
##
function runNWCHEM()
{
	local iFile=$1
	local nProcShared=$2
	
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$M3C_NWCHEM_HOME/lib/LINUX64
	local SCRATCH_DIR="$M3C_NWCHEM_SCRATCH/${iFile%.*}"
	
# 	if [ ! -d  "$SCRATCH_DIR" ]
# 	then
		mkdir -p $SCRATCH_DIR
# 	fi
	
	awk '
	BEGIN{
		line=1
	}
	( $1!~/SCRATCH_DIR/ || $1!~/scratch_dir/ || $1!~/PERMANENT_DIR/ || $1!~/permanent_dir/ ){
		map[line]=$0
		line++
	}
	END{
		print "SCRATCH_DIR '$SCRATCH_DIR'"
		print "PERMANENT_DIR '$SCRATCH_DIR'"
		for(i=1;i<line;i+=1)
			print map[i]
	}
	' $iFile > .$$tmp-$iFile
	mv .$$tmp-$iFile $iFile
	
	if [ "$nProcShared" -eq 1 ]
	then
		#$M3C_NWCHEM_HOME/bin/LINUX64/nwchem $iFile
		nwchem $iFile
	else
		# Using MPI
#		mpirun -n $nProcShared $M3C_NWCHEM_HOME/bin/LINUX64/nwchem $iFile
		# If you have all nodes connected via shared memory and you have installed the ch_shmem version of MPICH
 		$M3C_NWCHEM_HOME/bin/LINUX64/nwchem -np $nProcShared $iFile
# 		nwchem -np $nProcShared $iFile
	fi
	
	rm -rf $SCRATCH_DIR
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
		header="Geometry from NWCHEM"
	else
		header="Energy = $energy"
	fi
	
	nAtoms=`gawk 'BEGIN{i=0}($0!~/^[[:blank:]]*$/){i++}END{print i}' $iFile`

	echo "$nAtoms"
	echo "$header"
	gawk '{print $2"   "$4"   "$5"   "$6}' $iFile
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
function optgNWCHEMTemplate()
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
		fillTemplate $template $xyzFile $charge $mult > input$SID.nw
			
		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.nw ${xyzFile%.*}.nw 2> /dev/null
		
# 		if grep "Optimization converged" input$SID.out > /dev/null
		if grep "Please cite the following reference when publishing" input$SID.out > /dev/null
		then
			energy=`grep "Total DFT energy =" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
			grep -A$(( $nAtoms+3 )) "Output coordinates in angstroms" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			
			geom2xyz .finalGeom$SID $energy
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.nw input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsNWCHEMTemplate()
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
		fillTemplate $template $xyzFile $charge $mult > input$SID.nw
		
		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.nw ${xyzFile%.*}.nw 2> /dev/null
		
		if grep "Please cite the following reference when publishing" input$SID.out > /dev/null
		then
			energy=`grep "Total DFT energy =" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
			
			cat /dev/null > .freqs$SID
			
			grep "P.Frequency" input$SID.out | gawk '{for(i=2;i<=NF;i++) if(sqrt(($i-0.0)**2)>1.0) print $i}' > .freqs$SID
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
				echo "SYMMETRY R3"
				echo "ELECTRONIC_STATE ??"
			else
				group=`grep "Group name" input$SID.out | gawk '{print $3}'`
				if [ -z "$group" ]
				then
					echo "SYMMETRY ??"
				else
					group=${SYMMETRY_GROUP_MAP[$group]}
					echo "SYMMETRY $group"
				fi
				
				# @todo This is not available in NWCHEM
				state=`grep "The electronic state is" input$SID.out | gawk '{print $5}' | sed 's/\.//'`
				if [ -n "$state" ]
				then
					echo "ELECTRONIC_STATE $state"
				else
					echo "ELECTRONIC_STATE ??"
				fi
			fi
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	fi
	
	rm -rf .freqs$SID input$SID.nw input$SID.out
}

##
# @brief
##
function ienerNWCHEMTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
        echo "### Error ### ienerNWCHEMTemplate is not implemented yet"
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
# 		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.nw
# 		rm .xyzFile$SID
# 			
# 		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
# 		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
# 		cp input$SID.nw ${rxyzFile%.*}.nw 2> /dev/null
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
# 	rm -rf .finalGeom$SID input$SID.nw input$SID.out .finalGeom$SID
}
