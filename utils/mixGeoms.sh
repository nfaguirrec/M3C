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

##
# @brief
##
function main()
{
	local nSpinStates=$1
	
	[ -z "$nSpinStates" ] && nSpinStates=1
	
	mkdir backup-`date +%Y%m%d-%H.%M.%S` 2> /dev/null
	cp *.xyz backup-`date +%Y%m%d-%H.%M.%S`/
	
	categories=`ls *.xyz | sed 's/.*\///g;s/\..*//g' | gawk '{map[$1]=1}END{for(key in map) print key}'`
	
	for category in $categories
	do
		charges=`ls $category.* | awk 'BEGIN{FS="[.-]"}{map[$2]=1}END{ for( key in map ) print key }' | sed 's/q//g'`
		
		echo ""
		echo "---------------------"
		echo $category
		echo "---------------------"
		
		targetFiles=`ls $category.* 2> /dev/null`
		i=1
		for target in $targetFiles
		do
			echo "$target --> geom-$i.xyz"
			mv $target geom-$i.xyz
			
			i=$(( $i+1 ))
		done
		
		n=$(( $i-1 ))
		
		echo ""
		
		for i in `seq 1 $n`
		do
			echo -n "Replicating geom-$i.xyz ... "
			
			for charge in $charges
			do
				minMult=`molecule.minMult geom-$i.xyz $charge`
				maxMult=$(( $minMult+2*($nSpinStates-1) ))
				
				for multiplicity in `seq $minMult 2 $maxMult`
				do
					cp geom-$i.xyz $category.q$charge.m$multiplicity-$i.xyz
				done
			done
			
			echo "OK"
			
			rm geom-$i.xyz
		done
	done
}

main $*

