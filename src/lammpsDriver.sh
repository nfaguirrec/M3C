#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) by authors                                                          #
#                                                                                   #
#  Authors:                                                                         #
#                         * Néstor F. Aguirre (2018-2019)                           #
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

TEMPLATE_PATH=""

declare -A ATOMIC_NUMBER

# awk '($2~/Number/){number=$4}($2~/Symbol/){map[number]="ATOMIC_NUMBER[\""$4"\"]="number}END{for(k in map) print map[k]}' NIST.dat
ATOMIC_NUMBER["H"]=1
ATOMIC_NUMBER["He"]=2
ATOMIC_NUMBER["Li"]=3
ATOMIC_NUMBER["Be"]=4
ATOMIC_NUMBER["B"]=5
ATOMIC_NUMBER["C"]=6
ATOMIC_NUMBER["N"]=7
ATOMIC_NUMBER["O"]=8
ATOMIC_NUMBER["F"]=9
ATOMIC_NUMBER["Ne"]=10
ATOMIC_NUMBER["Na"]=11
ATOMIC_NUMBER["Mg"]=12
ATOMIC_NUMBER["Al"]=13
ATOMIC_NUMBER["Si"]=14
ATOMIC_NUMBER["P"]=15
ATOMIC_NUMBER["S"]=16
ATOMIC_NUMBER["Cl"]=17
ATOMIC_NUMBER["Ar"]=18
ATOMIC_NUMBER["K"]=19
ATOMIC_NUMBER["Ca"]=20
ATOMIC_NUMBER["Sc"]=21
ATOMIC_NUMBER["Ti"]=22
ATOMIC_NUMBER["V"]=23
ATOMIC_NUMBER["Cr"]=24
ATOMIC_NUMBER["Mn"]=25
ATOMIC_NUMBER["Fe"]=26
ATOMIC_NUMBER["Co"]=27
ATOMIC_NUMBER["Ni"]=28
ATOMIC_NUMBER["Cu"]=29
ATOMIC_NUMBER["Zn"]=30
ATOMIC_NUMBER["Ga"]=31
ATOMIC_NUMBER["Ge"]=32
ATOMIC_NUMBER["As"]=33
ATOMIC_NUMBER["Se"]=34
ATOMIC_NUMBER["Br"]=35
ATOMIC_NUMBER["Kr"]=36
ATOMIC_NUMBER["Rb"]=37
ATOMIC_NUMBER["Sr"]=38
ATOMIC_NUMBER["Y"]=39
ATOMIC_NUMBER["Zr"]=40
ATOMIC_NUMBER["Nb"]=41
ATOMIC_NUMBER["Mo"]=42
ATOMIC_NUMBER["Tc"]=43
ATOMIC_NUMBER["Ru"]=44
ATOMIC_NUMBER["Rh"]=45
ATOMIC_NUMBER["Pd"]=46
ATOMIC_NUMBER["Ag"]=47
ATOMIC_NUMBER["Cd"]=48
ATOMIC_NUMBER["In"]=49
ATOMIC_NUMBER["Sn"]=50
ATOMIC_NUMBER["Sb"]=51
ATOMIC_NUMBER["Te"]=52
ATOMIC_NUMBER["I"]=53
ATOMIC_NUMBER["Xe"]=54
ATOMIC_NUMBER["Cs"]=55
ATOMIC_NUMBER["Ba"]=56
ATOMIC_NUMBER["La"]=57
ATOMIC_NUMBER["Ce"]=58
ATOMIC_NUMBER["Pr"]=59
ATOMIC_NUMBER["Nd"]=60
ATOMIC_NUMBER["Pm"]=61
ATOMIC_NUMBER["Sm"]=62
ATOMIC_NUMBER["Eu"]=63
ATOMIC_NUMBER["Gd"]=64
ATOMIC_NUMBER["Tb"]=65
ATOMIC_NUMBER["Dy"]=66
ATOMIC_NUMBER["Ho"]=67
ATOMIC_NUMBER["Er"]=68
ATOMIC_NUMBER["Tm"]=69
ATOMIC_NUMBER["Yb"]=70
ATOMIC_NUMBER["Lu"]=71
ATOMIC_NUMBER["Hf"]=72
ATOMIC_NUMBER["Ta"]=73
ATOMIC_NUMBER["W"]=74
ATOMIC_NUMBER["Re"]=75
ATOMIC_NUMBER["Os"]=76
ATOMIC_NUMBER["Ir"]=77
ATOMIC_NUMBER["Pt"]=78
ATOMIC_NUMBER["Au"]=79
ATOMIC_NUMBER["Hg"]=80
ATOMIC_NUMBER["Tl"]=81
ATOMIC_NUMBER["Pb"]=82
ATOMIC_NUMBER["Bi"]=83
ATOMIC_NUMBER["Po"]=84
ATOMIC_NUMBER["At"]=85
ATOMIC_NUMBER["Rn"]=86
ATOMIC_NUMBER["Fr"]=87
ATOMIC_NUMBER["Ra"]=88
ATOMIC_NUMBER["Ac"]=89
ATOMIC_NUMBER["Th"]=90
ATOMIC_NUMBER["Pa"]=91
ATOMIC_NUMBER["U"]=92
ATOMIC_NUMBER["Np"]=93
ATOMIC_NUMBER["Pu"]=94
ATOMIC_NUMBER["Am"]=95
ATOMIC_NUMBER["Cm"]=96
ATOMIC_NUMBER["Bk"]=97
ATOMIC_NUMBER["Cf"]=98
ATOMIC_NUMBER["Es"]=99
ATOMIC_NUMBER["Fm"]=100
ATOMIC_NUMBER["Md"]=101
ATOMIC_NUMBER["No"]=102
ATOMIC_NUMBER["Lr"]=103
ATOMIC_NUMBER["Rf"]=104
ATOMIC_NUMBER["Db"]=105
ATOMIC_NUMBER["Sg"]=106
ATOMIC_NUMBER["Bh"]=107
ATOMIC_NUMBER["Hs"]=108
ATOMIC_NUMBER["Mt"]=109
ATOMIC_NUMBER["Ds"]=110
ATOMIC_NUMBER["Rg"]=111
ATOMIC_NUMBER["Cn"]=112
ATOMIC_NUMBER["Uut"]=113
ATOMIC_NUMBER["Fl"]=114
ATOMIC_NUMBER["Uup"]=115
ATOMIC_NUMBER["Lv"]=116
ATOMIC_NUMBER["Uus"]=117
ATOMIC_NUMBER["Uuo"]=118

ATOMIC_NUMBER["H"]=1
ATOMIC_NUMBER["HE"]=2
ATOMIC_NUMBER["LI"]=3
ATOMIC_NUMBER["BE"]=4
ATOMIC_NUMBER["B"]=5
ATOMIC_NUMBER["C"]=6
ATOMIC_NUMBER["N"]=7
ATOMIC_NUMBER["O"]=8
ATOMIC_NUMBER["F"]=9
ATOMIC_NUMBER["NE"]=10
ATOMIC_NUMBER["NA"]=11
ATOMIC_NUMBER["MG"]=12
ATOMIC_NUMBER["AL"]=13
ATOMIC_NUMBER["SI"]=14
ATOMIC_NUMBER["P"]=15
ATOMIC_NUMBER["S"]=16
ATOMIC_NUMBER["CL"]=17
ATOMIC_NUMBER["AR"]=18
ATOMIC_NUMBER["K"]=19
ATOMIC_NUMBER["CA"]=20
ATOMIC_NUMBER["SC"]=21
ATOMIC_NUMBER["TI"]=22
ATOMIC_NUMBER["V"]=23
ATOMIC_NUMBER["CR"]=24
ATOMIC_NUMBER["MN"]=25
ATOMIC_NUMBER["FE"]=26
ATOMIC_NUMBER["CO"]=27
ATOMIC_NUMBER["NI"]=28
ATOMIC_NUMBER["CU"]=29
ATOMIC_NUMBER["ZN"]=30
ATOMIC_NUMBER["GA"]=31
ATOMIC_NUMBER["GE"]=32
ATOMIC_NUMBER["AS"]=33
ATOMIC_NUMBER["SE"]=34
ATOMIC_NUMBER["BR"]=35
ATOMIC_NUMBER["KR"]=36
ATOMIC_NUMBER["RB"]=37
ATOMIC_NUMBER["SR"]=38
ATOMIC_NUMBER["Y"]=39
ATOMIC_NUMBER["ZR"]=40
ATOMIC_NUMBER["NB"]=41
ATOMIC_NUMBER["MO"]=42
ATOMIC_NUMBER["TC"]=43
ATOMIC_NUMBER["RU"]=44
ATOMIC_NUMBER["RH"]=45
ATOMIC_NUMBER["PD"]=46
ATOMIC_NUMBER["AG"]=47
ATOMIC_NUMBER["CD"]=48
ATOMIC_NUMBER["IN"]=49
ATOMIC_NUMBER["SN"]=50
ATOMIC_NUMBER["SB"]=51
ATOMIC_NUMBER["TE"]=52
ATOMIC_NUMBER["I"]=53
ATOMIC_NUMBER["XE"]=54
ATOMIC_NUMBER["CS"]=55
ATOMIC_NUMBER["BA"]=56
ATOMIC_NUMBER["LA"]=57
ATOMIC_NUMBER["CE"]=58
ATOMIC_NUMBER["PR"]=59
ATOMIC_NUMBER["ND"]=60
ATOMIC_NUMBER["PM"]=61
ATOMIC_NUMBER["SM"]=62
ATOMIC_NUMBER["EU"]=63
ATOMIC_NUMBER["GD"]=64
ATOMIC_NUMBER["TB"]=65
ATOMIC_NUMBER["DY"]=66
ATOMIC_NUMBER["HO"]=67
ATOMIC_NUMBER["ER"]=68
ATOMIC_NUMBER["TM"]=69
ATOMIC_NUMBER["YB"]=70
ATOMIC_NUMBER["LU"]=71
ATOMIC_NUMBER["HF"]=72
ATOMIC_NUMBER["TA"]=73
ATOMIC_NUMBER["W"]=74
ATOMIC_NUMBER["RE"]=75
ATOMIC_NUMBER["OS"]=76
ATOMIC_NUMBER["IR"]=77
ATOMIC_NUMBER["PT"]=78
ATOMIC_NUMBER["AU"]=79
ATOMIC_NUMBER["HG"]=80
ATOMIC_NUMBER["TL"]=81
ATOMIC_NUMBER["PB"]=82
ATOMIC_NUMBER["BI"]=83
ATOMIC_NUMBER["PO"]=84
ATOMIC_NUMBER["AT"]=85
ATOMIC_NUMBER["RN"]=86
ATOMIC_NUMBER["FR"]=87
ATOMIC_NUMBER["RA"]=88
ATOMIC_NUMBER["AC"]=89
ATOMIC_NUMBER["TH"]=90
ATOMIC_NUMBER["PA"]=91
ATOMIC_NUMBER["U"]=92
ATOMIC_NUMBER["NP"]=93
ATOMIC_NUMBER["PU"]=94
ATOMIC_NUMBER["AM"]=95
ATOMIC_NUMBER["CM"]=96
ATOMIC_NUMBER["BK"]=97
ATOMIC_NUMBER["CF"]=98
ATOMIC_NUMBER["ES"]=99
ATOMIC_NUMBER["FM"]=100
ATOMIC_NUMBER["MD"]=101
ATOMIC_NUMBER["NO"]=102
ATOMIC_NUMBER["LR"]=103
ATOMIC_NUMBER["RF"]=104
ATOMIC_NUMBER["DB"]=105
ATOMIC_NUMBER["SG"]=106
ATOMIC_NUMBER["BH"]=107
ATOMIC_NUMBER["HS"]=108
ATOMIC_NUMBER["MT"]=109
ATOMIC_NUMBER["DS"]=110
ATOMIC_NUMBER["RG"]=111
ATOMIC_NUMBER["CN"]=112
ATOMIC_NUMBER["UUT"]=113
ATOMIC_NUMBER["FL"]=114
ATOMIC_NUMBER["UUP"]=115
ATOMIC_NUMBER["LV"]=116
ATOMIC_NUMBER["UUS"]=117
ATOMIC_NUMBER["UUO"]=118

##
# @brief
##
function sortXYZbyAtomicNumber(){
	local xyzFile=$1
	
	declare -A idSort

	pairs=`gawk '(NR>2){print $1}' $xyzFile | nl | \
		while read id sym; do echo "$id ${ATOMIC_NUMBER[$sym]}" ; done | \
		sort -k2 -n | nl | awk '{printf $1":"$2"  "}'`
		
	for pair in $pairs
	do
		key=`echo $pair | sed 's/:.*//g'`
		value=`echo $pair | sed 's/.*://g'`
		
		idSort[$key]=$value
	done
	
	gawk '(NR<=2){print $0}' $xyzFile
	for k in `seq 1 ${#idSort[@]}`
	do
		gawk '(NR>2){print $0}' $xyzFile | nl | gawk '('${idSort[$k]}'==$1){print $2,$3,$4,$5}'
	done
}

##
# @brief
##
function formatLAMMPSoutput()
{
	local oFileRaw=$1
	local iFile=$2
	
	cat $oFileRaw
	
	dumpFile=`gawk '($0~/^[[:blank:]]*dump /){ print $NF; exit }' $iFile`
	
	nAtoms="NONE"
	if [ -f "$dumpFile" ]
	then
		nAtoms=`head -n1 $dumpFile`
		echo "NATOMS = $nAtoms"
		echo ""
		echo "LAST GEOMETRY"
		tail -n$(( $nAtoms )) $dumpFile
	fi
	
	if [ "$nAtoms" = "NONE" ]
	then
		echo "### ERROR ### Parsing final geometry from $dumpFile"
		exit 1
	fi

	echo ""
}

##
# @brief
##
function runLAMMPS()
{
	local iFile=$1
	local nProcShared=$2
	
	mkdir -p $M3C_LAMMPS_SCRATCH/${iFile%.*}
	
	awk '{ if( $0~/<<<<<<<<<<<<</ ) exit; else print $0 }' $iFile > $M3C_LAMMPS_SCRATCH/${iFile%.*}/lammps.in
	awk 'BEGIN{loc=0}{ if(loc==1) print $0; if( $0~/<<<<<<<<<<<<</ ) loc=1 }' $iFile > $M3C_LAMMPS_SCRATCH/${iFile%.*}/geometry.lmp
	
	latteLammpsInterface=`gawk '($0~/^[[:blank:]]*fix .* latte .*$/){print 1}' $iFile`
	if [ "$latteLammpsInterface" -eq 1 ]
	then
		paramPath=`grep "PARAMPATH=" $TEMPLATE_PATH/latte.in | gawk '{print $2}' | sed '{s/[\"]//g}'`
		paramPath="$TEMPLATE_PATH/$paramPath"
		
		cp $TEMPLATE_PATH/latte.in $M3C_LAMMPS_SCRATCH/${iFile%.*}/
		
        paramPathTMP=`echo "$paramPath" | sed 's/\//@@/g'`
        
        sed -i 's/PARAMPATH= .*$/PARAMPATH= "'$paramPathTMP'"/g' $M3C_LAMMPS_SCRATCH/${iFile%.*}/latte.in
        sed -i 's/@@/\//g' $M3C_LAMMPS_SCRATCH/${iFile%.*}/latte.in
	fi
	
	pushd . &> /dev/null
	cd $M3C_LAMMPS_SCRATCH/${iFile%.*}/
	
	export OMP_NUM_THREADS=1
	if [ -n "$nProcShared" ]
	then
		export OMP_NUM_THREADS=$nProcShared
	fi
	
	$M3C_LAMMPS_HOME/lmp_serial < lammps.in &> rawOutput
	formatLAMMPSoutput rawOutput lammps.in
	
	rm -rf $M3C_LAMMPS_SCRATCH/${iFile%.*}/
	popd &> /dev/null
}

##
# @brief
##
function xyz2LAMMPS_full()
{
	local iFile=$1
	
	# Ref: https://lammps.sandia.gov/doc/read_data.html
	# full:	atom-ID molecule-ID atom-type q x y z
	sortXYZbyAtomicNumber $iFile \
	| gawk '
		BEGIN{
			atomType=0
			lastsymbol="X"
			id=1
		}
		(NR>2){
			if( lastsymbol!=$1 ){
				lastsymbol=$1
				atomType+=1
			}
			printf( "%10d%5d%5d%5.1f%15.5f%15.5f%15.5f\n", id,1,atomType,0.0,$2,$3,$4 )
			id+=1
		}' 
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
		header="Geometry from LAMMPS"
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
	
	COORDSFILE=`grep -E "^[[:blank:]]*read_data " $template | awk '{print $2}' | sed 's/\"//g'`
	
	cat $template > template$SID
	echo "<<<<<<<<<<<<<"  >> template$SID
	cat $TEMPLATE_REALPATH/$COORDSFILE >> template$SID
	
	xyz2LAMMPS_full $xyzFile > .geom$SID
	nAtoms=`cat .geom$SID | wc -l`
	
	gawk '{
		if( $0~/@GEOMETRY_LAMMPS_FULL/ ){
			while( ( getline line < "'.geom$SID'" ) > 0 ) {
				print line
			}
		}else{
			print $0
		}
	}
	' template$SID > .tmp$SID
	
	sed -i 's/read_data .*$/read_data geometry.lmp/g' .tmp$SID
	
	nAtomTypes=`gawk '(NR>2){map[$1]=0}END{print length(map)}' $xyzFile`
	
	sed -i 's/@NATOMS/'$nAtoms'/g' .tmp$SID
	sed -i 's/@NATOM_TYPES/'$nAtomTypes'/g' .tmp$SID
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
function optgLAMMPSTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	TEMPLATE_PATH=`realpath $template`
	TEMPLATE_PATH=${TEMPLATE_PATH%/*}
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.lammps
			
		runLAMMPS input$SID.lammps $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.lammps ${xyzFile%.*}.lammps 2> /dev/null
		
		if [ -z "`grep "Total wall time:" input$SID.out`" ]
		then
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		else
			energy=`grep -A1 "Energy initial, next-to-last, final" input$SID.out | tail -n1 | awk '{print $NF*0.0367493088244753}'`
			awk 'BEGIN{loc=0}($0!~/^[[:blank:]]*$/){if(loc==1) print $0; if($0~/LAST GEOMETRY/) loc=1}' input$SID.out > .finalGeom$SID
			geom2xyz .finalGeom$SID $energy
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.lammps input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsLAMMPSTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	echo "### Error ### freqsLAMMPSTemplate is not implemented yet"
	exit
}

##
# @brief
##
function ienerLAMMPSTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	echo "### Error ### ienerLAMMPSTemplate is not implemented yet"
	exit
}

