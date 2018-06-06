#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) 2016-2016 Departamento de Química                                   #
#                         Universidad Autónoma de Madrid                            #
#                         All rights reserved.                                      #
#                                                                                   #
#                         * Néstor F. Aguirre (2016-2016)                           #
#                           nestor.aguirre@uam.es                                   #
#                         * Sergio Díaz-Tendero (2016-2016)                         #
#                           sergio.diaztendero@uam.es                               #
#                         * Manuel Alcamí (2016-2016)                               #
#                           manuel.alcami@uam.es                                    #
#                         * Fernando Martín (2016-2016)                             #
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

function parallel()
{
	iFile=$1
	nThreads=$2
	
	nCommands=`cat $iFile | wc -l`
	
	if [ "$nCommands" -lt 1 ]
	then
		echo "### ERROR ### parallel()"
		echo "              There is no any command to execute in '$iFile' file"
		exit
	fi
	
	startTime=`date "+%s"`
	
	ij=0
	for (( i=1; i<=$(( $nCommands/$nThreads )); i++ ))
	do
		iStartTime=`date "+%s"`
		
		echo "Running: "
		
		for (( j=1; j<=$nThreads; j++ ))
		do
			ij=$(( $j+$nThreads*$(($i-1)) ))
			
			echo -n "   $ij) "
			command=`sed ''$ij'q;d' $iFile`
			
			echo "$command"
			`echo "$command"` &
		done
		
		wait
		
		iEndTime=`date "+%s"`
		elapsedTime=$(( $iEndTime-$iStartTime ))
		echo "      `echo -n  "Time elapsed:"` $(( $elapsedTime / 3600 ))h $(( ( $elapsedTime / 60 ) % 60 ))m $(( $elapsedTime % 60 ))s"
	done
	
	if (( $ij < $nCommands ))
	then
		iStartTime=`date "+%s"`
		
		echo "Running: "
		
		k=1
		for (( i=$(( $ij+1 )); i<=$nCommands; i++ ))
		do
			echo -n "   $i) "
			command=`sed ''${i}'q;d' $iFile`
			
			echo "$command"
			`echo "$command"` &
			
			k=$(( $k+1 ))
		done
		
		for (( j=$k; j<=$nThreads; j++ ))
		do
			echo -n "   $i) "
			echo "None"
			i=$(( $i+1 ))
		done
		
		wait
		
		iEndTime=`date "+%s"`
		elapsedTime=$(( $iEndTime-$iStartTime ))
		echo "      `echo -n  "Time elapsed:"` $(( $elapsedTime / 3600 ))h $(( ( $elapsedTime / 60 ) % 60 ))m $(( $elapsedTime % 60 ))s"
	fi
	
	endTime=`date "+%s"`
	elapsedTime=$(( $endTime-$startTime ))
	
	echo ""
	echo "      ` echo -n  "Total"`: $(( $elapsedTime / 3600 ))h $(( ( $elapsedTime / 60 ) % 60 ))m $(( $elapsedTime % 60 ))s"
	
	#-----------------------------------------------------------------------
}

# cat > commands << EOF
# sleep 1s
# sleep 2s
# sleep 1s
# sleep 3s
# sleep 1s
# EOF
# 
# parallel commands 2

