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

iFile=$1

if [ -f "$M3C_HOME/bin/gaussianDriver.sh" ]
then
	source $M3C_HOME/bin/gaussianDriver.sh
elif [ -f "$M3C_HOME/src/gaussianDriver.sh" ]
then
	source $M3C_HOME/src/gaussianDriver.sh
else
	echo "### ERROR ### gaussianDriver.sh: No such file"
	echo "              Check \$M3C_HOME variable"
	exit
fi

nAtoms=`grep "NAtoms=" $iFile | awk '{print $2; exit}'`
( grep -A$(( $nAtoms+4 )) "Input orientation:" $iFile | tail -n$nAtoms | \
while read a1 a2 a3 a4 a5 a6
do
	printf "%5s%12.6f%12.6f%12.6f\n" ${ATOMIC_SYMBOL[$a2]} $a4 $a5 $a6
done ) > .geom

if [ "$nAtoms" -gt 0  ]
then
	energy=`grep -E "^[[:blank:]]+CCSD\(T\)= " $iFile | sed 's/D/E/g' | gawk '{printf "%.10f\n", $2}'`
	
	if [ -z "$energy" ]
	then
		energy=`grep "SCF Done" $iFile | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
	fi
	
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
		#group=`grep "Full point group" $iFile | tail -n1 | gawk '{print $4}'`
		group=`grep "Framework group" $iFile | sed 's/\[/ /g' | gawk '{print $3}' | tail -n1`
		if [ "$group" = "Nop" -o "$group" = "NOp" ]
		then
			echo "SYMMETRY ??"
		else
			echo "SYMMETRY $group"
		fi
		
		state=`grep "The electronic state is" $iFile | tail -n1 | gawk '{print $5}' | sed 's/\.//'`
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

