#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2016-2016)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2016-2016)
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

