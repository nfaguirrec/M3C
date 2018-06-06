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
SYMMETRY_GROUP_MAP["NOSYM"]="C1"
SYMMETRY_GROUP_MAP["C(I)"]="CI"
SYMMETRY_GROUP_MAP["C(S)"]="CS"
SYMMETRY_GROUP_MAP["C(LIN)"]="C*V"
SYMMETRY_GROUP_MAP["C(N)"]="C2"
SYMMETRY_GROUP_MAP["C(2V)"]="C2V"
SYMMETRY_GROUP_MAP["C(2H)"]="C2H"
SYMMETRY_GROUP_MAP["D(LIN)"]="D*H"
SYMMETRY_GROUP_MAP["C(3)"]="C3"
SYMMETRY_GROUP_MAP["C(3V)"]="C3V"
SYMMETRY_GROUP_MAP["C(3H)"]="C3H"
SYMMETRY_GROUP_MAP["S(6)"]="S6"  # <<< It looks is no implemented in ADF
SYMMETRY_GROUP_MAP["C(4)"]="C4"
SYMMETRY_GROUP_MAP["C(4V)"]="C4V"
SYMMETRY_GROUP_MAP["C(4H)"]="C4H"
SYMMETRY_GROUP_MAP["D(2)"]="D2"
SYMMETRY_GROUP_MAP["D(2D)"]="D2D"
SYMMETRY_GROUP_MAP["D(2H)"]="D2H"
SYMMETRY_GROUP_MAP["C(6)"]="C6"
SYMMETRY_GROUP_MAP["C(6V)"]="C6V"
SYMMETRY_GROUP_MAP["C(6H)"]="C6H"
SYMMETRY_GROUP_MAP["D(3)"]="D3"
SYMMETRY_GROUP_MAP["D(3D)"]="D3D"
SYMMETRY_GROUP_MAP["D(3H)"]="D3H"
SYMMETRY_GROUP_MAP["D(4)"]="D4"
SYMMETRY_GROUP_MAP["D(4D)"]="D4D"
SYMMETRY_GROUP_MAP["D(4H)"]="D4H"
SYMMETRY_GROUP_MAP["D(6)"]="D6"
SYMMETRY_GROUP_MAP["D(6D)"]="D6D"
SYMMETRY_GROUP_MAP["D(6H)"]="D6H"
SYMMETRY_GROUP_MAP["T"]="T"  # <<< It looks is no implemented in ADF
SYMMETRY_GROUP_MAP["T(D)"]="TD"
SYMMETRY_GROUP_MAP["O(H)"]="OH"

# ELECTRONIC_STATE_MAP[""]="BG"
# ELECTRONIC_STATE_MAP[""]="A1''"
# ELECTRONIC_STATE_MAP[""]="A"
# ELECTRONIC_STATE_MAP[""]="T"
# ELECTRONIC_STATE_MAP[""]="E''"
# ELECTRONIC_STATE_MAP[""]="EU"
# ELECTRONIC_STATE_MAP[""]="B2U"
# ELECTRONIC_STATE_MAP[""]="PIG"
# ELECTRONIC_STATE_MAP[""]="E1"
# ELECTRONIC_STATE_MAP[""]="A1U"
# ELECTRONIC_STATE_MAP[""]="E3"
# ELECTRONIC_STATE_MAP[""]="??"
# ELECTRONIC_STATE_MAP[""]="E"
# ELECTRONIC_STATE_MAP[""]="E1G"
# ELECTRONIC_STATE_MAP[""]="PI"
# ELECTRONIC_STATE_MAP[""]="E'"
# ELECTRONIC_STATE_MAP[""]="A1'"
# ELECTRONIC_STATE_MAP[""]="AU"
# ELECTRONIC_STATE_MAP[""]="A1"
# ELECTRONIC_STATE_MAP[""]="SGG"
# ELECTRONIC_STATE_MAP[""]="T1G"
# ELECTRONIC_STATE_MAP[""]="A'"
# ELECTRONIC_STATE_MAP[""]="B3U"
# ELECTRONIC_STATE_MAP[""]="DLTU"
# ELECTRONIC_STATE_MAP[""]="A2U"
# ELECTRONIC_STATE_MAP[""]="E2G"
# ELECTRONIC_STATE_MAP[""]="T1"
# ELECTRONIC_STATE_MAP[""]="A2'"
# ELECTRONIC_STATE_MAP[""]="BU"
# ELECTRONIC_STATE_MAP[""]="PHIG"
# ELECTRONIC_STATE_MAP[""]="B1"
# ELECTRONIC_STATE_MAP[""]="E2''"
# ELECTRONIC_STATE_MAP[""]="PIU"
# ELECTRONIC_STATE_MAP[""]="T2G"
# ELECTRONIC_STATE_MAP[""]="B1G"
# ELECTRONIC_STATE_MAP[""]="B3"
# ELECTRONIC_STATE_MAP[""]="A1G
# ELECTRONIC_STATE_MAP[""]="E1''"
# ELECTRONIC_STATE_MAP[""]="E1U"
# ELECTRONIC_STATE_MAP[""]="B"
# ELECTRONIC_STATE_MAP[""]="SGU"
# ELECTRONIC_STATE_MAP[""]="A''"
# ELECTRONIC_STATE_MAP[""]="E2"
# ELECTRONIC_STATE_MAP[""]="T1U"
# ELECTRONIC_STATE_MAP[""]="E1'"
# ELECTRONIC_STATE_MAP[""]="EG"
# ELECTRONIC_STATE_MAP[""]="B2G"
# ELECTRONIC_STATE_MAP[""]="A1G"
# ELECTRONIC_STATE_MAP[""]="A2"
# ELECTRONIC_STATE_MAP[""]="DLTA"
# ELECTRONIC_STATE_MAP[""]="E2U"
# ELECTRONIC_STATE_MAP[""]="SG"
# ELECTRONIC_STATE_MAP[""]="PHIU"
# ELECTRONIC_STATE_MAP[""]="AG"
# ELECTRONIC_STATE_MAP[""]="T2U"
# ELECTRONIC_STATE_MAP[""]="B1U"
# ELECTRONIC_STATE_MAP[""]="E2'"
# ELECTRONIC_STATE_MAP[""]="PHI"
# ELECTRONIC_STATE_MAP[""]="B3G"
# ELECTRONIC_STATE_MAP[""]="T2"
# ELECTRONIC_STATE_MAP[""]="DLTG"
# ELECTRONIC_STATE_MAP[""]="A2G"
# ELECTRONIC_STATE_MAP[""]="A2''"
# ELECTRONIC_STATE_MAP[""]="B2"

# OH  
# # "A1G"
# # "A1U"
# # "A2G"
# # "A2U"
# # "EG"
# # "EU"
# # "T1G"
# # "T1U"
# # "T2G"
# # "T2U"
# # "??"
# # 
# O   
# # "A1"
# # "A2"
# # "E"
# # "T1"
# # "T2"
# # "??"
# # 
# TD  
# # "A1"
# # "A2"
# # "E"
# # "T1"
# # "T2"
# # "??"
# # 
# T   
# # "A"
# # "E"
# # "T"
# # "??"
# # 
# # 
# D*  D*H 
# # "SGG"
# # "SGU"
# # "PIG"
# # "PIU"
# # "DLTG"
# # "DLTU"
# # "PHIG"
# # "PHIU"
# # "??"
# # 
# D2H
# # "AG"
# # "AU"
# # "B1G"
# # "B1U"
# # "B2G"
# # "B2U"
# # "B3G"
# # "B3U"
# # "??"
# # 
# D3H
# # "A1''"
# # "A1'"
# # "A2''"
# # "A2'"
# # "E''"
# # "E'"
# # "??"
# # 
# D4H
# # "A1G"
# # "A1U"
# # "A2G"
# # "A2U"
# # "B1G"
# # "B1U"
# # "B2G"
# # "B2U"
# # "EG"
# # "EU"
# # "??"
# # 
# D5H
# # "A1''"
# # "A1'"
# # "A2''"
# # "A2'"
# # "E1''"
# # "E1'"
# # "E2''"
# # "E2'"
# # "??"
# # 
# D6H
# # "A1G"
# # "A1U"
# # "A2G"
# # "A2U"
# # "B1G"
# # "B1U"
# # "B2G"
# # "B2U"
# # "E1G"
# # "E1U"
# # "E2G"
# # "E2U"
# # "??"
# # 
# D2D
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E"
# # "??"
# # 
# D3D
# # "A1G "
# # "A1U"
# # "A2G"
# # "A2U"
# # "EG"
# # "EU"
# # "??"
# # 
# D4D
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E1"
# # "E2"
# # "E3"
# # "??"
# # 
# D2 
# # "A"
# # "B1"
# # "B2"
# # "B3"
# # "??"
# # 
# D3 
# # "A1"
# # "A2"
# # "E"
# # "??"
# # 
# D4 
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E"
# # "??"
# # 
# D6
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E1"
# # "E2"
# # "??"
# # 
# # 
# C*V C*  
# # "SG"
# # "PI"
# # "DLTA"
# # "PHI"
# # "??"
# # 
# C2H
# # "AG"
# # "AU"
# # "BG"
# # "BU"
# # "??"
# # 
# C3H
# # "A''"
# # "A'"
# # "E''"
# # "E'"
# # "??"
# # 
# C4H
# # "AG"
# # "AU"
# # "BG"
# # "BU"
# # "EG"
# # "EU"
# # "??"
# # 
# C6H
# # "AG"
# # "AU"
# # "BG"
# # "BU"
# # "E1G"
# # "E1U"
# # "E2G"
# # "E2U"
# # "??"
# # 
# C2V
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "??"
# # 
# C3V
# # "A1"
# # "A2"
# # "E"
# # "??"
# # 
# C4V
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E"
# # "??"
# # 
# C5V
# # "A1"
# # "A2"
# # "E1"
# # "E2"
# # "??"
# # 
# C6V
# # "A1"
# # "A2"
# # "B1"
# # "B2"
# # "E1"
# # "E2"
# # "??"
# # 
# CS  
# # "A''"
# # "A'"
# # "??"
# # 
# CI  
# # "AG"
# # "AU"
# # "??"
# # 
# C2 
# # "A"
# # "B"
# # "??"
# # 
# C3 
# # "A"
# # "E"
# # "??"
# # 
# C4 
# # "A"
# # "B"
# # "E"
# # "??"
# # 
# C6 
# # "A"
# # "B"
# # "E1"
# # "E2"
# # "??"
# # 
# OTHER
# # "A1"
# # "??"

##
# @brief
##
function runADF()
{
	local program=$1
	local iFile=$2
	local nProcShared=$3
	local copyT21to=$4
	
	source /home/nestor/Software/adf2017.106/adfrc.sh
	export M3C_ADF_HOME=$HOME/Software/adf2017.106
	M3C_ADF_SCRATCH=/scratch/$USER/ADF
	
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
		$M3C_ADF_HOME/bin/$program -n $nProcShared < $iFile
	else
		$M3C_ADF_HOME/bin/$program < $iFile
	fi
	
	if [ "$program" = "adf" -a -n "$copyT21to" ]
	then
		cp TAPE21 $copyT21to 2> /dev/null
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
	local program=$1
	local template=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
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
	sed -i 's/@NSPIN/'$(($mult-1))'/g' .tmp$SID
	
	if [ "$program" = "adf" -a "$mult" -ne "1" ]
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
	
	local program="Unknown"
	if [ -z "`grep -m1 -E "^[[:blank:]]*DFTB[[:blank:]]*$" $template`" ]
	then
		program="adf"
	else
		program="dftb"
	fi
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $program $template $xyzFile $charge $mult > input$SID.adf
			
		runADF $program input$SID.adf $nProcShared $PWD/${xyzFile%.*}.t21 &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.adf ${xyzFile%.*}.adf 2> /dev/null
		
# 		if grep "ERROR:" input$SID.out > /dev/null
		if [ -n "`grep "ERROR" input$SID.out`" ]
		then
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		else
# 			energy=`grep "<.*Total energy" input$SID.out | tail -n1 | gawk '{print $5}'`  # <<< This not available in all cases, I think is just atoms
			
			if [ "$program" = "adf" ]
			then
			
				energy=`grep "<.*current energy" input$SID.out | tail -n1 | gawk '{print $5}'`
				grep -A$(( $nAtoms+1 )) "Coordinates in Geometry Cycle" input$SID.out | tail -n$nAtoms | sed -r 's/^[[:blank:]]*[[:digit:]]+\.//g' > .finalGeom$SID
				geom2xyz .finalGeom$SID $energy
				
			# We can afford geometry convergence errors (looking for ERROR). It is not necessary to recover any good geometry candidate
			elif [ "$program" = "dftb" ]
			then
			
				energy=`grep "Total Energy (hartree)" input$SID.out | tail -n1 | awk '{print $NF}'`
				grep -A$nAtoms "Index Symbol   x (angstrom)   y (angstrom)   z (angstrom)" input$SID.out | tail -n$nAtoms | gawk '{$1=""; print $0}' > .finalGeom$SID
				geom2xyz .finalGeom$SID $energy
			fi
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
	
	local SID="-$xyzFile$RANDOM"
	
	local program="Unknown"
	if [ -z "`grep -m1 -E "^[[:blank:]]*DFTB[[:blank:]]*$" $template`" ]
	then
		program="adf"
	else
		program="dftb"
	fi
	
	local nAtoms=""
	local energy=""
	local fv="" # Vibrational degrees of freedom
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		fillTemplate $program $template $xyzFile $charge $mult > input$SID.adf
		
		runADF $program input$SID.adf $nProcShared $PWD/${xyzFile%.*}.t21 &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.adf ${xyzFile%.*}.adf 2> /dev/null
		
# 		if grep "ERROR:" input$SID.out > /dev/null
		if [ -n "`grep "ERROR" input$SID.out`" ]
		then
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		else
# 			energy=`grep "<.*Total energy" input$SID.out | tail -n1 | gawk '{print $5}'`  # <<< This not available in all cases, I think is just atoms
			
			if [ "$program" = "adf" ]
			then
				energy=`grep "<.*current energy" input$SID.out | tail -n1 | gawk '{print $5}'`
				
			# We can afford geometry convergence errors (looking for ERROR). It is not necessary to recover any good geometry candidate
# 			elif [ "$program" = "dftb" -a -z "`grep "ERROR DETECTED" input$SID.out`" ]
			elif [ "$program" = "dftb" ]
			then
				energy=`grep "Total Energy (hartree)" input$SID.out | tail -n1 | awk '{print $NF}'`
			fi
			
			cat /dev/null > .freqs$SID
			
			awk '
				($0~/Index:.*Frequency \(cm-1\)/){print $5}
			' input$SID.out >> .freqs$SID
			
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
				group=`grep "Symmetry:" input$SID.out | tail -n1 | gawk '{print $2}'`
				
				[ -z "$group" ] && group="NOSYM" # < Just to fix a warning message in dftb. dftb does not use symmetry
				
				group=${SYMMETRY_GROUP_MAP[$group]}
				if [ "$group" = "Nop" -o "$group" = "NOp" ]
				then
					echo "SYMMETRY ??"
				else
					echo "SYMMETRY $group"
				fi
				
				# @todo This is not available in ADF
				state=`grep "The electronic state is" input$SID.out | gawk '{print $5}' | sed 's/\.//'`
				if [ -n "$state" ]
				then
					echo "ELECTRONIC_STATE $state"
				else
					echo "ELECTRONIC_STATE ??"
				fi
			fi
		fi
	fi
	
	rm -rf .freqs$SID input$SID.adf input$SID.out
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
}
