#!/bin/bash

iFile=$1

gawk '
	BEGIN{
		n=1
	}
	
	( NR==1 ){
		nAtoms=$1
	}
	
	(NR==2){
		title=$0
		print nAtoms
		print title
	}
	
	( NR>2 && n<=nAtoms ){
		print $0
		n++
	}
' $iFile