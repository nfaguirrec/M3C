#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
#                                                                                   #
# Copyright (c) 2020-2020 by authors                                                #
# Authors:                                                                          #
#                         * Néstor F. Aguirre (2020-2020)                           #
#                           nfaguirrec@gmail.com                                    #
#                                                                                   #
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

if [ -f "$M3C_HOME/bin/parallel.sh" ]
then
	source $M3C_HOME/bin/parallel.sh
else
	source $M3C_HOME/src/parallel.sh
fi

runM3C(){
	local iFile=$1
	local excitationEnergy=$2
	local numberOfEvents=$3
	local keepOutputFiles=$4

	# @todo Si excitationEnergy no está en el input el programa se cae. Hay que verificar primero su existencia
	# @todo Hay que hacer que no solo busque los casos " = ", sino todos los "\s+=\s+"
	sed -r -i 's/(^[[:blank:]]*)excitationEnergy = .*$/\1excitationEnergy = '${excitationEnergy}' # eV  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)numberOfEvents = [0-9]+[[:blank:]]*$/\1numberOfEvents = '${numberOfEvents}' # <-- Generated by M3C.p/g' $iFile
	
	sed -r -i 's/(^[[:blank:]]*)geometryHistoryFilePrefix = .*$/\1geometryHistoryFilePrefix = E_'${excitationEnergy}'-geom #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)energyHistoryFile = .*$/\1energyHistoryFile = E_'${excitationEnergy}'.ehistory #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)KERHistoryFilePrefix = .*$/\1KERHistoryFilePrefix = E_'${excitationEnergy}'-KER #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)weightHistoryFile = .*$/\1weightHistoryFile = E_'${excitationEnergy}'.whistory #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)JHistoryFile = .*$/\1JHistoryFile = E_'${excitationEnergy}'.jhistory #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)LHistoryFile = .*$/\1LHistoryFile = E_'${excitationEnergy}'.lhistory #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)histogramFile = .*$/\1histogramFile = E_'${excitationEnergy}'.histogram #  <-- Generated by M3C.p/g' $iFile
	sed -r -i 's/(^[[:blank:]]*)(END MARKOV_CHAIN[[:blank:]]*$)/\1\tgenEbkl = TRUE  #  <-- Generated by M3C.p\n\1\2/g' $iFile
	
	if [ "$keepOutputFiles" = "1" ]
	then
		M3C -i $iFile > ${iFile%.*}.out
	else
		M3C -i $iFile > /dev/null
	fi
}

main(){
	local i=0
	local j=0
	local ij=0
	local iFile=""
	local nThreads=""
	local outputDir=""
# 	local scratch=""
# 	local work=""
	local iFileEff=""

	NPROCSHARED=1  # M3C is not parallelized
	
	if [ -z "$M3C_NTHREADS" ]
	then
		M3C_NTHREADS=`cat /proc/cpuinfo | grep processor | wc -l`
	fi
	nThreads=$(( $M3C_NTHREADS/$NPROCSHARED ))
	
	while getopts "i:n:" OPTNAME
	do
		case $OPTNAME in
				"i" )
						iFile=$OPTARG
						;;
				"n" )
						nThreads=$OPTARG
						;;
				* )
						exit
						;;
		esac
	done
	
	if [ -z "$iFile" ]
	then
			echo "Usage:"
			echo "      $ M3C.p -i iFile.inp [ -n nThreads ]"
			exit
	fi
	
	outputDir="$PWD/${iFile%.*}.data"
	
	#--------------------------------------------------------------------------------------------------------------------------------
	# La variable de salida de este bloque es el vector excitationEnergy el resto de 
	# variables definidas no deberían de utilizarse más adelante
	#--------------------------------------------------------------------------------------------------------------------------------
	if [ -z "`grep "EXCITATION_ENERGY_SCAN" $iFile`" ]
	then
		echo "### ERROR ### M3C.p. Check your input file $iFile, EXCITATION_ENERGY_SCAN block is required."
		exit 0
	fi
	gawk 'BEGIN{ loc=0 }{ if($0~/END EXCITATION_ENERGY_SCAN/) loc=0; if(loc==1) print $0; if( $0~/BEGIN EXCITATION_ENERGY_SCAN/ ) loc=1 }' $iFile > .$iFile-erange$$
	
	keepOutputFiles=`gawk '($1~/keepOutputFiles/){ print $3 }' .$iFile-erange$$`

	#-------------------------------------------------------
	# Parsing energy grid
	#   output: excitationEnergy array
	#-------------------------------------------------------
	gridType=`gawk '($1~/excitationEnergy/){ split($3,arr,":"); if( arr[1]=="file" ) print "file"; else print "range" }' .$iFile-erange$$`
	
	if [ "$gridType" = "range" ]
	then
	
		minValue=`gawk '($1~/excitationEnergy/){ split($3,arr,":"); print arr[1] }' .$iFile-erange$$`
		maxValue=`gawk '($1~/excitationEnergy/){ split($3,arr,":"); print arr[2] }' .$iFile-erange$$`
		nPoints=`gawk '($1~/excitationEnergy/){ split($3,arr,":"); print arr[3] }' .$iFile-erange$$`
		
		stepSize=`echo "($maxValue-($minValue))/($nPoints-1.0)" | bc -l`
		
		excitationEnergy=( `seq -s " " -f "%10.5f" $minValue $stepSize $maxValue` )
		
	elif [ "$gridType" = "file" ]
	then
	
		gridFile=`gawk '($1~/excitationEnergy/){ split($3,arr,":"); print arr[2] }' .$iFile-erange$$`
		
		excitationEnergy=( `gawk '($1!~/^[[:blank:]]*#/){ print $1 }' $gridFile` )
	fi

	#-------------------------------------------------------
	# Parsing numberOfEvents grid
	#   output: numberOfEvents array
	#-------------------------------------------------------
	gridType=`gawk '($1~/numberOfEvents/){ split($3,arr,":"); if( arr[1]=="file" ) print "file"; else print "range" }' .$iFile-erange$$`
	
	if [ "$gridType" = "range" ]
	then
		
		minValue=`gawk '($1~/numberOfEvents/){ split($3,arr,":"); print arr[1] }' .$iFile-erange$$`
		maxValue=`gawk '($1~/numberOfEvents/){ split($3,arr,":"); print arr[2] }' .$iFile-erange$$`
		nPoints=`gawk '($1~/numberOfEvents/){ split($3,arr,":"); print arr[3] }' .$iFile-erange$$`
		
		stepSize=`echo "($maxValue-($minValue))/($nPoints-1.0)" | bc -l`
		
		numberOfEvents=( `seq -s " " -f "%10.0f" $minValue $stepSize $maxValue` )
		
	elif [ "$gridType" = "file" ]
	then
	
		gridFile=`gawk '($1~/numberOfEvents/){ split($3,arr,":"); print arr[2] }' .$iFile-erange$$`
		
		numberOfEvents=( `gawk '($1!~/^[[:blank:]]*#/){ print $1 }' $gridFile` )
	else
		numberOfEventsBase=`gawk 'BEGIN{ loc=0 }{ if($0~/END MARKOV_CHAIN/) loc=0; if(loc==1) print $0; if( $0~/BEGIN MARKOV_CHAIN/ ) loc=1 }' $iFile | gawk '($1~/numberOfEvents/){print $3}'`
		numberOfEvents=( `echo $numberOfEventsBase | gawk '{ for(i=1;i<='${#excitationEnergy[@]}';i++) print $1 }'` )
	fi
	
	#-------------------------------------------------------
	# Checking consistency and cleaning
	#-------------------------------------------------------
	
	rm .$iFile-erange$$
	
	if [ "${#excitationEnergy[@]}" -ne "${#numberOfEvents[@]}" ]
	then
		echo "### ERROR ### EXCITATION_ENERGY_SCAN inconsistency. Size of grids excitationEnergy and numberOfEvents are not the same."
		exit
	fi
	
	#--------------------------------------------------------------------------------------------------------------------------------
	
	startTime=`date "+%s"`
	
	pushd . > /dev/null 2> /dev/null
	
	if [ -d "$outputDir" ]
	then
		echo "### ERROR ### There is already an output directory ($outputDir), please remove or rename it before to run the calculation"
		exit
	else
		mkdir $outputDir
	fi
	cp $iFile $outputDir
	cp *.xyz $outputDir &> /dev/null
	cp *.rxyz $outputDir &> /dev/null
	cp *.molden $outputDir &> /dev/null
	cd $outputDir

	cat /dev/null > .commands$$
	for (( i=0; i<${#excitationEnergy[@]}; i++ ))
	do
		iFileEff="E_${excitationEnergy[$i]}.inp"
		cp $iFile $iFileEff
		echo "runM3C $iFileEff ${excitationEnergy[$i]} ${numberOfEvents[$i]} $keepOutputFiles" >> .commands$$
	done
	
	parallel .commands$$ $nThreads
	rm .commands$$
	
	#-----------------------------------------------------------------------
}

main $*
