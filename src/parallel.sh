#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2016-2017)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2016-2017)
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
