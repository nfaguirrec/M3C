#/bin/bash

iFile=$1
direction=$2

function usage()
{
	echo "Usage: oADF2xyz-noImagFreqs.sh ADF_output [ direction ]"
	echo "                                               +1      "
	echo ""
	exit -1
}

[ -z "$iFile" ] && usage
[ -z "$direction" ] && direction=1
[ "$direction" != "-1" -a "$direction" != "+1" -a "$direction" != "1" ] && usage

optg=`grep "Coordinates in Geometry Cycle" $iFile | wc -l`
awk '
BEGIN{
        loc=0
        n = 0
}
{
        if($0~/^ </) loc=0
        if(loc==1 &&$1~/^[[:digit:]]+./){
                gsub("^[[:digit:]]+.","",$1)
                line[n] = $0
                n++
        }
        if( '$optg' == 0 && $0~/>>>> FRAGM/){
                loc=1
                n = 0
        }else if( '$optg' != 0 && $0~/Coordinates in Geometry Cycle/){
                loc=1
                n = 0
        }
}
END{
        for( k in line )
                print line[k]
}
' $iFile > .finalGeom$$
nAtoms=`cat .finalGeom$$ | wc -l`

awk '
BEGIN{
	loc=0
}
{
	if($0~/List of All Frequencies/)
		loc=0
	if(loc==1 && $1~/^[[:digit:]]+.[[:alpha:]]+$/)
		print $0
	if(loc==1 && $1~/^[-]?[[:digit:]]+.[[:digit:]]+$/)
		for( i=1; i<=NF; i++ ) print $i > "/dev/stderr"
	if($0~/Vibrations and Normal Modes/)
		loc=1
}
' $iFile > .modes$$ 2> .freqs$$

python -c "
import numpy

nAtoms=int('$nAtoms')

coords = []
iFile = open('.finalGeom$$')
for n,line in enumerate(iFile.readlines()):
	arr = line.split()
	coords.append( [ arr[0], numpy.array(map(float, line.split()[1:])) ] )
iFile.close()

freqs = []
iFile = open('.freqs$$')
for n,line in enumerate(iFile.readlines()):
	freqs.append( float(line) )
iFile.close()

iFile = open('.modes$$')
modes = len(freqs)*[ None ]

m=0
for n,line in enumerate(iFile.readlines()):
	arr = map(float, line.split()[1:])
	modesPerLine = len(arr)/3
	
	j = 0
	for i in range(0,len(arr),3):
		if( modes[m+j] is None ):
			modes[m+j] = []
		modes[m+j].append( numpy.array(arr[i:i+3]) )
# 		print m, j, m+j, modes[m+j]
		j += 1
		
	if( (n+1)%nAtoms == 0 ):
		m += modesPerLine
		
iFile.close()

first = True
for i,mode in enumerate(modes):
	if( mode is not None and freqs[i]<0.0 ):
		if( first ): # Largest negative frequency
			for j,atom in enumerate(coords):
# 				print j, i, freqs[i], mode[j][0:3]
				atom[1] = atom[1] + 0.3*float('$direction')*mode[j]
		else:
			for j,atom in enumerate(coords):
# 				print j, i, freqs[i], mode[j][0:3]
				atom[1] = atom[1] + 0.3*mode[j]
# 			print ""

		first = False
			
print len(coords)
print ""
for atom in coords:
	print atom[0], atom[1][0], atom[1][1], atom[1][2]
"

rm .finalGeom$$ .freqs$$ .modes$$
