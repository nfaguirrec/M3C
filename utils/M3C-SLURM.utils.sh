#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2013-2016)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2013-2016)
#          nestor.aguirre@uam.es
#    * Dr. Sergio Díaz-Tendero (2015-2015)
#          sergio.diaztendero@uam.es
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

test -z "$M3CSLURM_GAUSSIAN_HOME"     && M3CSLURM_GAUSSIAN_HOME=$HOME/.gaussian
test -z "$M3CSLURM_GAUSSIAN_SCRATCH"  && M3CSLURM_GAUSSIAN_SCRATCH=/scratch/$USER/gaussian

test -z "$M3CSLURM_GAMESS_HOME"       && M3CSLURM_GAMESS_HOME=$HOME/.gamess
test -z "$M3CSLURM_GAMESS_SCRATCH"    && M3CSLURM_GAMESS_SCRATCH=/scratch/$USER/gamess

##
# BASIC CONFIGURATION
#
function SLURM.buildHead()
{
	local queueParams=$1
	local name=$2
	
	local partition=""
	local nTask=""
	local account=""
	local ttime=""
	local values=""
	
	local partition=`echo $queueParams | awk 'BEGIN{FS=","}{ print $1 }'`
	local nTask=`echo $queueParams | awk 'BEGIN{FS=","}{ print $2 }'`
	local account=`echo $queueParams | awk 'BEGIN{FS=","}{ print $3 }'`
	local ttime=`echo $queueParams | awk 'BEGIN{FS=","}{ print $4 }'`
	local qos=`echo $queueParams | awk 'BEGIN{FS=","}{ print $5 }'`
	
	if [ -z "$partition" ]
	then
		echo "### ERROR ### partition parameter is requiered (queue name)" > /dev/stderr
		echo "              COMMAND <partition> [OPTIONS]" > /dev/stderr
		return 1
	fi
	
	if [ -f "$HOME/.SLURM.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $HOME/.SLURM.properties`
	elif [ -f "$M3C_HOME/utils/SLURM.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $M3C_HOME/utils/SLURM.properties`
	else
		echo "### ERROR ### SLURM default parameters file not found" > /dev/stderr
		return 1
	fi
	
	[ -z "$partition" ] && partition=`echo $values | awk '{print $1}'`
	[ -z "$nTask" ]     && nTask=`echo $values | awk '{print $2}'`
	[ -z "$account" ]   && account=`echo $values | awk '{print $3}'`
	[ -z "$ttime" ]     && ttime=`echo $values | awk '{print $4}'`
	[ -z "$qos" ]       && qos=`echo $values | awk '{print $5}'`

	[ -z "$qos" ] && qos=$partition
	
	local jobdir=`echo $PWD | sed s/.*$USER/~/`
	
	[ ! -d log ] && mkdir log
	
	cat << EOF
#!/bin/bash
#SBATCH --partition=$partition
#SBATCH --ntasks=$nTask
#SBATCH --account=$account
#SBATCH --time=$ttime
#SBATCH --qos=$qos
#SBATCH --job-name=$jobdir/$name
#SBATCH -o log/$name.slurm.log
#SBATCH -e log/$name.slurm.err

EOF
	
	return 0
}

##
# M3C-gamess.geniso CONFIGURATION
#
function SLURM.M3C-gamess.geniso()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export PATH=$M3CSLURM_GAMESS_HOME:\$PATH

M3C-gamess.geniso $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gamess.optg CONFIGURATION
#
function SLURM.M3C-gamess.optg()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export PATH=$M3CSLURM_GAMESS_HOME:\$PATH

M3C-gamess.optg $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gamess.freqs CONFIGURATION
#
function SLURM.M3C-gamess.freqs()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export PATH=$M3CSLURM_GAMESS_HOME:\$PATH

M3C-gamess.freqs $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.geniso CONFIGURATION
#
function SLURM.M3C-gaussian.geniso()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export GAUSS_EXEDIR=$M3CSLURM_GAUSSIAN_HOME
export PATH=\$GAUSS_EXEDIR:\$PATH
export GAUSS_SCRDIR=$M3CSLURM_GAUSSIAN_SCRATCH
if [ ! -d  "\$GAUSS_SCRDIR" ]
then
        mkdir -p \$GAUSS_SCRDIR
fi

M3C-gaussian.geniso $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.optg CONFIGURATION
#
function SLURM.M3C-gaussian.optg()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export GAUSS_EXEDIR=$M3CSLURM_GAUSSIAN_HOME
export PATH=\$GAUSS_EXEDIR:\$PATH
export GAUSS_SCRDIR=$M3CSLURM_GAUSSIAN_SCRATCH
if [ ! -d  "\$GAUSS_SCRDIR" ]
then
        mkdir -p \$GAUSS_SCRDIR
fi

M3C-gaussian.optg $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.freqs CONFIGURATION
#
function SLURM.M3C-gaussian.freqs()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
export GAUSS_EXEDIR=$M3CSLURM_GAUSSIAN_HOME
export PATH=\$GAUSS_EXEDIR:\$PATH
export GAUSS_SCRDIR=$M3CSLURM_GAUSSIAN_SCRATCH
if [ ! -d  "\$GAUSS_SCRDIR" ]
then
        mkdir -p \$GAUSS_SCRDIR
fi

M3C-gaussian.freqs $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C CONFIGURATION
#
function SLURM.M3C()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C.p CONFIGURATION
#
function SLURM.M3C.p()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C.p $* > SLURM.log 2> SLURM.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

