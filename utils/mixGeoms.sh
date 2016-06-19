#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2012-2016)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2012-2016)
#          nestor.aguirre@uam.es
#    * Dr. Sergio Díaz-Tendero (2012-2015)
#          sergio.diaztendero@uam.es
#    * Prof. M. Paul-Antoine Hervieux (2012-2015)
#          Paul-Antoine.Hervieux@ipcms.unistra.fr
#    * Prof. Fernando Martín (2012-2015)
#          fernando.martin@uam.es
#    * Prof. Manuel Alcamí (2012-2015)
#          manuel.alcami@uam.es
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
	mkdir backup-`date +%Y%m%d-%H.%M.%S` 2> /dev/null
	cp *.xyz backup-`date +%Y%m%d-%H.%M.%S`/
	
	categories=`ls *.xyz | sed 's/.*\///g;s/\..*//g' | gawk '{map[$1]=1}END{for(key in map) print key}'`
	
	for category in $categories
	do
		charges=`ls $category.* | awk 'BEGIN{FS="[.-]"}{map[$2]=1}END{ for( key in map ) print key }'`
		multiplicities=`ls $category.* | awk 'BEGIN{FS="[.-]"}{map[$3]=1}END{ for( key in map ) print key }'`
		
		echo ""
		echo "---------------------"
		echo $category
		echo "---------------------"
		
		targetFiles=`ls $category.* 2> /dev/null`
		i=1
		for target in $targetFiles
		do
			echo "$target --> geom-$i"
			mv $target geom-$i
			
			i=$(( $i+1 ))
		done
		
		n=$(( $i-1 ))
		
		echo ""
		
		for i in `seq 1 $n`
		do
			echo -n "Replicating geom-$i ... "
			
			for charge in $charges
			do
				for multiplicity in $multiplicities
				do
					cp geom-$i $category.$charge.$multiplicity-$i.xyz
				done
			done
			
			echo "OK"
			
			rm geom-$i
		done
	done
}

main $*

