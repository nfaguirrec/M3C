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

REFRESH_TIME="2s"

if [ -f "$M3C_HOME/bin/stack.sh" ]
then
        source $M3C_HOME/bin/stack.sh
elif [ -f "$M3C_HOME/src/stack.sh" ]
then
        source $M3C_HOME/src/stack.sh
else
        echo "### ERROR ### stack.sh: No such file"
        echo "              Check \$M3C_HOME variable"
        exit
fi

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
	
	stack_new jobs
	commands=()
	pidJobs=()
	iStartTime=()
	
	for (( i=1; i<=$nCommands; i++ ))
	do
		command=`sed ''$i'q;d' $iFile`
		stack_push jobs "$command"
	done
	
	firstTime=1
	while true
	do
		stack_size jobs sizeJobs
		
		if [ "$sizeJobs" -eq 0 ]
		then
			break
		fi
		
		if [ "${#pidJobs[@]}" -le "$nThreads" -a "$sizeJobs" -ne 0 ]
		then
			stack_pop jobs command
			
			if [ -n "$command" ]
			then
				`echo "$command"` &
				pid="$!"
				
				pidJobs[$pid]=$pid
				commands[$pid]=$command
				echo "Running: $command ($pid)"
				iStartTime[$pid]=`date "+%s"`
			fi
		fi
		
		if [ "$firstTime" -eq 1 ]
		then
			if [ "${#pidJobs[@]}" -eq "$nThreads" ]
			then
				firstTime=0
			else
				continue
			fi
		fi
		
		if [ "${#pidJobs[@]}" -le "$nThreads" ]
		then
			while true
			do
				delete=()
				for pid in ${pidJobs[@]}
				do
					loc=0
					for cpid in `ps -U $USER | gawk '($1~/^[[:digit:]]+$/){print $1}'`
					do
						if [ "$pid" -eq "$cpid" ]
						then
							loc=1
							break
						fi
					done
					
					if [ $loc -eq 0 ]
					then
						delete[$pid]=$pid
					fi
				done
				
				if [ "${#delete[@]}" -ne 0 ]
				then
					break
				fi
				
				sleep $REFRESH_TIME
			done
		fi
		
		if [ "${#delete[@]}" -ne 0 ]
		then
			for pid in ${delete[@]}
			do
				echo -n "   Finished: ${commands[$pid]} ($pid)  "
				unset pidJobs[$pid]
				
				iEndTime=`date "+%s"`
				elapsedTime=$(( $iEndTime-${iStartTime[$pid]} ))
				echo "      `echo -n  "Time elapsed:"` $(( $elapsedTime / 3600 ))h $(( ( $elapsedTime / 60 ) % 60 ))m $(( $elapsedTime % 60 ))s"
				unset iStartTime[$pid]
			done
		fi
	done
	
	endTime=`date "+%s"`
	elapsedTime=$(( $endTime-$startTime ))
	
	echo ""
	echo "` echo -n  "Total"`: $(( $elapsedTime / 3600 ))h $(( ( $elapsedTime / 60 ) % 60 ))m $(( $elapsedTime % 60 ))s"
}

# cat > commands << EOF
# sleep 1s
# sleep 2s
# sleep 3s
# sleep 4s
# sleep 5s
# sleep 4s
# sleep 3s
# sleep 2s
# sleep 1s
# EOF
# 
# parallel commands 3
# rm commands
# 
