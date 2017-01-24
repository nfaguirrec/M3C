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
function runGAMESS()
{
	local iFile=$1
	local nProcShared=$2
	
	$M3C_GAMESS_HOME/rungms $iFile 01 $nProcShared
}

##
# @brief
##
function xyz2geom()
{
	local iFile=$1
	
	local SID="-$iFile$RANDOM"
	
	local label=""
	local id=""
	local nAtoms=""
	local i=""
	
	cp $iFile .tmp$SID
	
	for i in "${!ATOMIC_NUMBER[@]}"
	do
		label=$i
		id=${ATOMIC_NUMBER[$i]}
		
		sed -i 's/^[[:blank:]]*'$label' / '$label'   '$id'./g' .tmp$SID
	done
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' .tmp$SID`
	tail -n $nAtoms .tmp$SID
	
	rm .tmp$SID
}

##
# @brief
##
function geom2xyz()
{
	local iFile=$1
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}($0!~/^[[:blank:]]*$/){i++}END{print i}' $iFile`
	
	echo "$nAtoms"
	echo ""
	gawk '{print $1"   "$3"   "$4"   "$5}' $iFile
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
	local oshell=$5
	
	local SID="-$xyzFile$RANDOM"
	
	if [ -z "$charge" ]
	then
		charge="0"
	fi
	
	if [ -z "$mult" ]
	then
		mult="1"
	fi
	
	xyz2geom $xyzFile > .geom$SID
	
	gawk '{
		if( $0~/@GEOMETRY/ ){
			while( ( getline line < "'.geom$SID'" ) > 0 ) {
				print line
			}
		}else{
			print $0
		}
	}
	' $template > .tmp$SID
	
	sed -i 's/@CHARGE/'$charge'/g' .tmp$SID
	
	if [ "$mult" -eq 1 ]
	then
		sed -i 's/@MULT/'$mult' scftyp=rhf/g' .tmp$SID
	else
		if [ -z "$oshell" ]
		then
			sed -i 's/@MULT/'$mult' scftyp=uhf/g' .tmp$SID
		else
			sed -i 's/@MULT/'$mult' scftyp='$oshell'/g' .tmp$SID
		fi
	fi
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID
}

##
# @brief
##
function optgGAMESSTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	local oshell=$6
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		rm -rf input$SID.dat 2> /dev/null
		fillTemplate $template $xyzFile $charge $mult $oshell > input$SID.inp
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${xyzFile%.*}.inp 2> /dev/null
		
		if grep "EXECUTION OF GAMESS TERMINATED NORMALLY" input$SID.out > /dev/null
		then
			grep -A $(( $nAtoms + 3 )) "***** EQUILIBRIUM GEOMETRY LOCATED *****" input$SID.out \
				| tail -n $nAtoms > .finalGeom$SID
			
			geom2xyz .finalGeom$SID
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.inp input$SID.dat input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsGAMESSTemplate()
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
		rm -rf input$SID.dat
		fillTemplate $template $xyzFile $charge $mult > input$SID.inp
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${xyzFile%.*}.inp 2> /dev/null
		
		energy=`grep "FINAL .* ENERGY IS" input$SID.out | head -n1 | gawk '{print $5}'`
		fv=`molecule.fv $xyzFile`
		
		echo $nAtoms
		echo "Energy = $energy"
		cat $xyzFile | gawk '(NR>2){print $0}'
		echo ""
		grep "       FREQUENCY:" input$SID.out \
			| gawk '
				BEGIN{n=0}{for(i=2;i<=NF;i++) if($i>100.0){ val[n]=$i; n++ }}
				END{
					asort(val)
					fv='`echo $fv`'
					print "FREQUENCIES   ", fv
					for (i = n; i > sqrt( (fv-n)**2 ); i--) print val[i]
				}'
	fi
	
	rm -rf input$SID.dat input$SID.inp input$SID.out input$SID.rst
}

##
# @brief
##
function ienerGAMESSTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$rxyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		rm -rf input$SID.dat
		
		rxyz2xyz $rxyzFile > .xyzFile$SID
		fillTemplate $template .xyzFile$SID $charge $mult rohf > input$SID.inp
		rm .xyzFile$SID
		
		runGAMESS input$SID.inp $nProcShared > input$SID.out 2>&1
		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
		cp input$SID.inp ${rxyzFile%.*}.inp 2> /dev/null
		
		if grep "EXECUTION OF GAMESS TERMINATED NORMALLY" input$SID.out > /dev/null
		then
			energy=`grep "COUPLED-CLUSTER ENERGY" input$SID.out | gawk '{print $5}'`  # RHF
			
			if [ -z "$energy" ]
			then
				energy=`grep "CCSD ENERGY" input$SID.out | gawk '{print $3}'`  # ROHF
			fi
			
			sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
		else
			echo "***** FAILURE CONVERGE *****"
		fi
	fi
	
	rm -rf input$SID.dat input$SID.inp input$SID.out
}
