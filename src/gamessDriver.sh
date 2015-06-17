#!/bin/bash

declare -A ATOMIC_NUMBER

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

##
# @brief
##
function runGAMESS()
{
	local iFile=$1
	
	rungms $iFile 01 1
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
	
	gawk '{if($0~/^[[:blank:]]*$/) exit; print $0}' $iFile
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
	local xyzFile=$2
	local charge=$3
	local mult=$4
	local oshell=$5
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 1  ]
	then
		rm -rf input$SID.dat 2> /dev/null
		fillTemplate $template $xyzFile $charge $mult $oshell > input$SID.inp
		
		runGAMESS input$SID.inp > input$SID.out 2>&1
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
	local xyzFile=$2
	local charge=$3
	local mult=$4
	
	local SID="-$xyzFile$RANDOM"
	
	local nAtoms=""
	local energy=""
	local fv="" # Vibrational degrees of freedom
	
	nAtoms=`gawk 'BEGIN{i=0}(NR>2 && $0!~/^[[:blank:]]*$/){i++}END{print i}' $xyzFile`
	
	if [ "$nAtoms" -gt 0  ]
	then
		rm -rf input$SID.dat
		fillTemplate $template $xyzFile $charge $mult > input$SID.inp
		
		runGAMESS input$SID.inp > input$SID.out 2>&1
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
	local rxyzFile=$2
	local charge=$3
	local mult=$4
	
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
		
		runGAMESS input$SID.inp > input$SID.out 2>&1
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
