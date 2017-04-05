#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2017-2017)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2017-2017)
#          nfaguirrec@lanl.gov
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

# awk '($2~/Number/){number=$4}($2~/Symbol/){map[number]="ATOMIC_SYMBOL["number"]=\""$4"\""}END{for(k in map) print map[k]}' NIST.dat
ATOMIC_SYMBOL[1]="H"
ATOMIC_SYMBOL[2]="He"
ATOMIC_SYMBOL[3]="Li"
ATOMIC_SYMBOL[4]="Be"
ATOMIC_SYMBOL[5]="B"
ATOMIC_SYMBOL[6]="C"
ATOMIC_SYMBOL[7]="N"
ATOMIC_SYMBOL[8]="O"
ATOMIC_SYMBOL[9]="F"
ATOMIC_SYMBOL[10]="Ne"
ATOMIC_SYMBOL[11]="Na"
ATOMIC_SYMBOL[12]="Mg"
ATOMIC_SYMBOL[13]="Al"
ATOMIC_SYMBOL[14]="Si"
ATOMIC_SYMBOL[15]="P"
ATOMIC_SYMBOL[16]="S"
ATOMIC_SYMBOL[17]="Cl"
ATOMIC_SYMBOL[18]="Ar"
ATOMIC_SYMBOL[19]="K"
ATOMIC_SYMBOL[20]="Ca"
ATOMIC_SYMBOL[21]="Sc"
ATOMIC_SYMBOL[22]="Ti"
ATOMIC_SYMBOL[23]="V"
ATOMIC_SYMBOL[24]="Cr"
ATOMIC_SYMBOL[25]="Mn"
ATOMIC_SYMBOL[26]="Fe"
ATOMIC_SYMBOL[27]="Co"
ATOMIC_SYMBOL[28]="Ni"
ATOMIC_SYMBOL[29]="Cu"
ATOMIC_SYMBOL[30]="Zn"
ATOMIC_SYMBOL[31]="Ga"
ATOMIC_SYMBOL[32]="Ge"
ATOMIC_SYMBOL[33]="As"
ATOMIC_SYMBOL[34]="Se"
ATOMIC_SYMBOL[35]="Br"
ATOMIC_SYMBOL[36]="Kr"
ATOMIC_SYMBOL[37]="Rb"
ATOMIC_SYMBOL[38]="Sr"
ATOMIC_SYMBOL[39]="Y"
ATOMIC_SYMBOL[40]="Zr"
ATOMIC_SYMBOL[41]="Nb"
ATOMIC_SYMBOL[42]="Mo"
ATOMIC_SYMBOL[43]="Tc"
ATOMIC_SYMBOL[44]="Ru"
ATOMIC_SYMBOL[45]="Rh"
ATOMIC_SYMBOL[46]="Pd"
ATOMIC_SYMBOL[47]="Ag"
ATOMIC_SYMBOL[48]="Cd"
ATOMIC_SYMBOL[49]="In"
ATOMIC_SYMBOL[50]="Sn"
ATOMIC_SYMBOL[51]="Sb"
ATOMIC_SYMBOL[52]="Te"
ATOMIC_SYMBOL[53]="I"
ATOMIC_SYMBOL[54]="Xe"
ATOMIC_SYMBOL[55]="Cs"
ATOMIC_SYMBOL[56]="Ba"
ATOMIC_SYMBOL[57]="La"
ATOMIC_SYMBOL[58]="Ce"
ATOMIC_SYMBOL[59]="Pr"
ATOMIC_SYMBOL[60]="Nd"
ATOMIC_SYMBOL[61]="Pm"
ATOMIC_SYMBOL[62]="Sm"
ATOMIC_SYMBOL[63]="Eu"
ATOMIC_SYMBOL[64]="Gd"
ATOMIC_SYMBOL[65]="Tb"
ATOMIC_SYMBOL[66]="Dy"
ATOMIC_SYMBOL[67]="Ho"
ATOMIC_SYMBOL[68]="Er"
ATOMIC_SYMBOL[69]="Tm"
ATOMIC_SYMBOL[70]="Yb"
ATOMIC_SYMBOL[71]="Lu"
ATOMIC_SYMBOL[72]="Hf"
ATOMIC_SYMBOL[73]="Ta"
ATOMIC_SYMBOL[74]="W"
ATOMIC_SYMBOL[75]="Re"
ATOMIC_SYMBOL[76]="Os"
ATOMIC_SYMBOL[77]="Ir"
ATOMIC_SYMBOL[78]="Pt"
ATOMIC_SYMBOL[79]="Au"
ATOMIC_SYMBOL[80]="Hg"
ATOMIC_SYMBOL[81]="Tl"
ATOMIC_SYMBOL[82]="Pb"
ATOMIC_SYMBOL[83]="Bi"
ATOMIC_SYMBOL[84]="Po"
ATOMIC_SYMBOL[85]="At"
ATOMIC_SYMBOL[86]="Rn"
ATOMIC_SYMBOL[87]="Fr"
ATOMIC_SYMBOL[88]="Ra"
ATOMIC_SYMBOL[89]="Ac"
ATOMIC_SYMBOL[90]="Th"
ATOMIC_SYMBOL[91]="Pa"
ATOMIC_SYMBOL[92]="U"
ATOMIC_SYMBOL[93]="Np"
ATOMIC_SYMBOL[94]="Pu"
ATOMIC_SYMBOL[95]="Am"
ATOMIC_SYMBOL[96]="Cm"
ATOMIC_SYMBOL[97]="Bk"
ATOMIC_SYMBOL[98]="Cf"
ATOMIC_SYMBOL[99]="Es"
ATOMIC_SYMBOL[100]="Fm"
ATOMIC_SYMBOL[101]="Md"
ATOMIC_SYMBOL[102]="No"
ATOMIC_SYMBOL[103]="Lr"
ATOMIC_SYMBOL[104]="Rf"
ATOMIC_SYMBOL[105]="Db"
ATOMIC_SYMBOL[106]="Sg"
ATOMIC_SYMBOL[107]="Bh"
ATOMIC_SYMBOL[108]="Hs"
ATOMIC_SYMBOL[109]="Mt"
ATOMIC_SYMBOL[110]="Ds"
ATOMIC_SYMBOL[111]="Rg"
ATOMIC_SYMBOL[112]="Cn"
ATOMIC_SYMBOL[113]="Uut"
ATOMIC_SYMBOL[114]="Fl"
ATOMIC_SYMBOL[115]="Uup"
ATOMIC_SYMBOL[116]="Lv"
ATOMIC_SYMBOL[117]="Uus"
ATOMIC_SYMBOL[118]="Uuo"

##
# @brief
##
function runNWCHEM()
{
	local iFile=$1
	local nProcShared=$2
	
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$M3C_NWCHEM_HOME/lib/LINUX64
	if [ ! -d  "$M3C_NWCHEM_SCRATCH" ]
	then
		mkdir -p $M3C_NWCHEM_SCRATCH
	fi
	
	awk '
	BEGIN{
		line=1
	}
	( $1!~/SCRATCH_DIR/ || $1!~/scratch_dir/ ){
		map[line]=$0
		line++
	}
	END{
		print "SCRATCH_DIR '$SCRATCH_DIR'"
		for(i=1;i<line;i+=1)
			print map[i]
	}
	' $iFile > .$$tmp-$iFile
	mv .$$tmp-$iFile $iFile
	
	echo "@@@ M3C:WARNING @@@ \$nProcShared option is not enabled yet"
	
	$M3C_NWCHEM_HOME/bin/LINUX64/nwchem $iFile
}

##
# @brief
##
function xyz2geom()
{
	local iFile=$1
	
	gawk 'BEGIN{i=0}( NR>2 && $0!~/^[[:blank:]]*$/ ){ print $0 }' $iFile
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
		header="Geometry from NWCHEM"
	else
		header="Energy = $energy"
	fi
	
	nAtoms=`gawk 'BEGIN{i=0}($0!~/^[[:blank:]]*$/){i++}END{print i}' $iFile`

	echo "$nAtoms"
	echo "$header"
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
	sed -i 's/@MULT/'$mult'/g' .tmp$SID
	
	cat .tmp$SID
	
	rm .geom$SID .tmp$SID
}

##
# @brief
##
function optgNWCHEMTemplate()
{
	local template=$1
	local nProcShared=$2
	local xyzFile=$3
	local charge=$4
	local mult=$5
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		fillTemplate $template $xyzFile $charge $mult > input$SID.nw
			
		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.nw ${xyzFile%.*}.nw 2> /dev/null
		
		if grep "Normal termination" input$SID.out > /dev/null
		then
			energy=`grep "SCF Done" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
# 			grep -A$(( $nAtoms+4 )) "Standard orientation:" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			grep -A$(( $nAtoms+4 )) "Input orientation:" input$SID.out | tail -n$nAtoms > .finalGeom$SID
			
			geom2xyz .finalGeom$SID $energy
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	else
		cat $xyzFile
	fi
	
	rm -rf .finalGeom$SID input$SID.nw input$SID.out .finalGeom$SID
}

##
# @brief
##
function freqsNWCHEMTemplate()
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
		fillTemplate $template $xyzFile $charge $mult > input$SID.nw
		
		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
		cp input$SID.out ${xyzFile%.*}.out 2> /dev/null
		cp input$SID.nw ${xyzFile%.*}.nw 2> /dev/null
		
		if grep "Normal termination" input$SID.out > /dev/null
		then
			energy=`grep "SCF Done" input$SID.out | tail -n 1 | cut -d "=" -f 2 | cut -d "A" -f 1`
			
			cat /dev/null > .freqs$SID
			
			grep "Frequencies" input$SID.out | while read a1 a2 freq1 freq2 freq3
			do
				if [ -n "$freq1" ]
				then
					echo $freq1 >> .freqs$SID
				fi
				
				if [ -n "$freq2" ]
				then
					echo $freq2 >> .freqs$SID
				fi
				
				if [ -n "$freq3" ]
				then
					echo $freq3 >> .freqs$SID
				fi
			done
			fv=`cat .freqs$SID | wc -l`
			
			echo $nAtoms
			echo "Energy = $energy"
			cat $xyzFile | gawk '(NR>2){print $0}'
			echo ""
			
			echo "FREQUENCIES $fv"
			cat .freqs$SID
			echo ""
			
			if [ "$nAtoms" -eq 1  ]
			then
				echo "SYMMETRY R3"
				echo "ELECTRONIC_STATE ??"
			else
				group=`grep "Full point group" input$SID.out | gawk '{print $4}'`
				if [ "$group" = "Nop" -o "$group" = "NOp" ]
				then
					echo "SYMMETRY ??"
				else
					echo "SYMMETRY $group"
				fi
				
				state=`grep "The electronic state is" input$SID.out | gawk '{print $5}' | sed 's/\.//'`
				if [ -n "$state" ]
				then
					echo "ELECTRONIC_STATE $state"
				else
					echo "ELECTRONIC_STATE ??"
				fi
			fi
		else
			echo "***** FAILURE TO LOCATE STATIONARY POINT, TOO MANY STEPS TAKEN *****"
		fi
	fi
	
	rm -rf .freqs$SID input$SID.nw input$SID.out
}

##
# @brief
##
function ienerNWCHEMTemplate()
{
	local template=$1
	local nProcShared=$2
	local rxyzFile=$3
	local charge=$4
	local mult=$5
	
        echo "### Error ### ienerNWCHEMTemplate is not implemented yet"
        exit
	
# 	local SID="-$rxyzFile$RANDOM"
# 	
# 	local nAtoms=""
# 	
# 	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $rxyzFile`
# 	
# 	if [ "$nAtoms" -gt 0  ]
# 	then
# 		rxyz2xyz $rxyzFile > .xyzFile$SID
# 		fillTemplate $template .xyzFile$SID $charge $mult > input$SID.nw
# 		rm .xyzFile$SID
# 			
# 		runNWCHEM input$SID.nw $nProcShared &> input$SID.out
# 		cp input$SID.out ${rxyzFile%.*}.out 2> /dev/null
# 		cp input$SID.nw ${rxyzFile%.*}.nw 2> /dev/null
# 		
# 		if grep "Normal termination" input$SID.out > /dev/null
# 		then
# 			energy=`grep -E "^[[:blank:]]+CCSD\(T\)= " input$SID.out | sed 's/D/E/g' | gawk '{printf "%.10f\n", $2}'`
# 			
# 			sed 's/Energy = .*/Energy = '$energy'/g' $rxyzFile
# 		else
# 			echo "***** FAILURE CONVERGE *****"
# 		fi
# 	fi
# 	
# 	rm -rf .finalGeom$SID input$SID.nw input$SID.out .finalGeom$SID
}
