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

##
# @brief
##
function formatLATTEoutput()
{
	local oFileRaw=$1
	
	cat $oFileRaw
	
	nAtoms="NONE"
	if [ -f "monitorrelax.xyz" ]
	then
		nAtoms=`head -n1 monitorrelax.xyz`
		echo "NATOMS = $nAtoms"
		echo ""
		echo "LAST GEOMETRY"
		tail -n$(( $nAtoms )) monitorrelax.xyz
	elif [ -f "restart_singlepoint.dat" ]
	then
		nAtoms=`head -n1 restart_singlepoint.dat`
		echo "NATOMS = $nAtoms"
		echo ""
		echo "LAST GEOMETRY"
		head -n$(( $nAtoms+4 )) restart_singlepoint.dat | sed '{1,4d}'
	fi
	
	if [ "$nAtoms" = "NONE" ]
	then
		echo "### ERROR ### Parsing final geometry from monitorrelax.xyz or restart_singlepoint"
		exit 1
	fi

	echo ""
}

##
# @brief
##
function runLATTE()
{
	local iFile=$1
	local nProcShared=$2
	
	mkdir -p $M3C_LATTE_SCRATCH/${iFile%.*}
	
	awk '{ if( $0~/<<<<<<<<<<<<</ ) exit; else print $0 }' $iFile > $M3C_LATTE_SCRATCH/${iFile%.*}/latte.in
	awk 'BEGIN{loc=0}{ if(loc==1) print $0; if( $0~/<<<<<<<<<<<<</ ) loc=1 }' $iFile > $M3C_LATTE_SCRATCH/${iFile%.*}/geometry.dat
	
	pushd . &> /dev/null
	cd $M3C_LATTE_SCRATCH/${iFile%.*}/
	
	export OMP_NUM_THREADS=1
	if [ -n "$nProcShared" ]
	then
		export OMP_NUM_THREADS=$nProcShared
	fi
	
	$M3C_LATTE_HOME/LATTE_DOUBLE &> rawOutput
	formatLATTEoutput rawOutput
	
	rm -rf $M3C_LATTE_SCRATCH/${iFile%.*}/
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
		header="Geometry from LATTE"
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
	
	TEMPLATE_REALPATH=`realpath $template`
	TEMPLATE_REALPATH=${TEMPLATE_REALPATH%/*}
	
	PARAMPATH=`grep "PARAMPATH=" $template | awk '{print $2}' | sed 's/\"//g'`
	COORDSFILE=`grep "COORDSFILE=" $template | awk '{print $2}' | sed 's/\"//g'`
	
	cat $template > template$SID
	echo "<<<<<<<<<<<<<"  >> template$SID
	cat $TEMPLATE_REALPATH/$COORDSFILE >> template$SID
	
	xyz2geom $xyzFile > .geom$SID
	nAtoms=`cat .geom$SID | wc -l`
	
	gawk '{
		if( $0~/@GEOMETRY/ ){
			while( ( getline line < "'.geom$SID'" ) > 0 ) {
				print line
			}
		}else{
			print $0
		}
	}
	' template$SID > .tmp$SID
	
	PARAMPATH=`echo "$TEMPLATE_REALPATH/$PARAMPATH" | sed 's/\//@@/g'`
	sed -i 's/PARAMPATH= .*$/PARAMPATH= "'$PARAMPATH'"/g' .tmp$SID
	sed -i 's/@@/\//g' .tmp$SID
	
	sed -i 's/COORDSFILE= .*$/COORDSFILE= geometry.dat/g' .tmp$SID
	
	sed -i 's/@NATOMS/'$nAtoms'/g' .tmp$SID
	sed -i 's/@CHARGE/'$charge'/g' .tmp$SID
	sed -i 's/@MULT/'$mult'/g' .tmp$SID
	sed -i 's/@NSPIN/'$(($mult-1))'/g' .tmp$SID
	
	if [ "$mult" -ne "1" ]
	then
		sed -ri 's/^[[:blank:]]*SPINON= [01]+[[:blank:]]+/ SPINON= 1 /gI' .tmp$SID
	fi
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID template$SID
}

##
# @brief
##
function optgLATTETemplate()
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
		fillTemplate $template $xyzFile $charge $mult > input$SID.latte
			
		runLATTE input$SID.latte $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.latte ${xyzFile%.*}.latte 2> /dev/null
		
# 		if grep "ERROR:" input$SID.out > /dev/null
		if [ -n "`grep "### ERROR ###" input$SID.out`" ]
		then
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		else
			energy=`grep "FREE ENERGY" input$SID.out | awk '{print $NF}'`
			awk 'BEGIN{loc=0}($0!~/^[[:blank:]]*$/){if(loc==1) print $0; if($0~/LAST GEOMETRY/) loc=1}' input$SID.out > .finalGeom$SID
			geom2xyz .finalGeom$SID $energy
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.latte input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsLATTETemplate()
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
		fillTemplate $template $xyzFile $charge $mult > input$SID.latte
		
		runLATTE input$SID.latte $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.latte ${xyzFile%.*}.latte 2> /dev/null
		
# 		if grep "ERROR:" input$SID.out > /dev/null
		if [ -n "`grep "ERROR:" input$SID.out`" ]
		then
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		else
			cat /dev/null > .freqs$SID
			
			energy=`grep "<.*Bond Energy.*a.u." input$SID.out | tail -n1 | gawk '{print $5}'`  # <<< When only a frequencies calculation was carried out.
			
			awk '
				BEGIN{ loc=0 }
				(loc==2&&$0~/^[[:blank:]]*$/){ loc=0 }
				(loc==2){ print $1 }
				($0~/List of All Frequencies:/){loc=1}
				(loc==1&&$1=="----------"){loc=2}
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
				
				# @todo This is not available in LATTE
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
	
	rm -rf .freqs$SID input$SID.latte input$SID.out
}

##
# @brief
##
function ienerLATTETemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	echo "### Error ### ienerLATTETemplate is not implemented yet"
	exit
}