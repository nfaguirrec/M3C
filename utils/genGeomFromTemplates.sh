#!/bin/bash

#  Este programa utiliza todos los esqueletos moleculares que encuentre en el
#  directorio TEMPLATES_DIR y utiliza esas geometrias paero cambiando la identidad
#  de los atomos utilizando todas las combinaciones que se dan en ATOMS_LIST
# 
# genGeomFromTemplates.sh directory O,O,O,O,Th,Th > salida.xyz
# splitXYZ.sh salida.xyz && rm salida.xyz

TEMPLATES_DIR=$1
ATOMS_LIST=$2

for f in `ls $TEMPLATES_DIR/*.xyz`
do
	nAtomsTemplate=`head -n1 $f`
	nAtomsTarget=`echo $ATOMS_LIST | sed 's/,/ /g' | wc -w`
	
	if [ "$nAtomsTemplate" -eq "$nAtomsTarget" ]
	then
		echo "$f --> $nAtomsTemplate, $nAtomsTarget" > /dev/stderr
		
		gawk '
		function buildPermutations( n,     i, j, k, s, a, id ) {
			id = 1
			
			if ((!k) || (k>n)) k=n
				
			for (j=1; j<=n; j++) {
				a[j]=0
				i[j]=0
			}
			
			a[n+1]=0 # n+1 is always free
			
			# initialize the selected field with 0
			j=1
			while ( j != 0 ) {
				s=i[j]
				
				if (s==0) i[j]=1
					
				while (a[i[j]]) i[j]++
					
				if (i[j]==n+1) {
					i[j]=0
					a[s]=0
					j--
					continue
				}
				
				if ( s!=0 ) a[s]=0
				a[i[j]]=1
				
				if (j==k) {
					value = ""
					for (l=1; l<k; l++) value = value""i[l]","
					value = value""i[k]
					
					permutations[id] = value
					id++
					
					continue
				}
				
				j++
			}
		}
		
		BEGIN{
			i=1
			split( "'$ATOMS_LIST'", atomsList, "," )
			
			buildPermutations( length(atomsList) )
		}
		
		(NR==1){
			nAtoms=$1
		}
		
		(NR>2 && i<=nAtoms){
			coordinates[i]=$2"  "$3"  "$4; i++
		}
		
		END{
			for( p in permutations ){
				split( permutations[p], pos, "," )
				
				print nAtoms
				print "Permutation ", permutations[p]
				for( i=1; i<=nAtoms; i++ )
					print atomsList[pos[i]], coordinates[i]
			}
		}
		' $f
	fi
done

