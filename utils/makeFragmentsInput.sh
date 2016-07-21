#/bin/bash

maxCharge=$1

[ -z "$maxCharge" ] && maxCharge=1

printf "#%19s%8s%5s\n" "XYZfile" "charge" "mult"

for q in `seq 0 $maxCharge`
do
	for f in `ls init/*.xyz`
	do
		name=`echo $f | sed 's/init\\///g'`
		mult=`molecule.minMult $f`
		
		printf "%20s%8d%5d\n" $name $q $mult
	done
	
	echo ""
	
	for f in `ls init/*.xyz`
	do
		name=`echo $f | sed 's/init\\///g'`
		mult=`molecule.minMult $f`
		
		printf "%20s%8d%5d\n" $name $q $(( $mult + 2 ))
	done
	
	echo ""
done