#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2016-2017)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2016-2017)
#          nestor.aguirre@uam.es
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

#############################
# M3C_GAUSSIAN_HOME
# M3C_GAUSSIAN_SCRATCH
# M3C_GAMESS_HOME
# M3C_GAMESS_SCRATCH
#############################

##
# BASIC CONFIGURATION
#
function PBS.buildHead()
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
	
	if [ -f "$HOME/.PBS.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $HOME/.PBS.properties`
	elif [ -f "$M3C_HOME/utils/PBS.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $M3C_HOME/utils/PBS.properties`
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
#PBS -q $partition
#PBS -l nodes=1:ppn=$nTask,walltime=$ttime
###PBS --account=$account # There is not equivalence for QOS in PBS
###PBS --qos=$qos # There is not equivalence for QOS in PBS
#PBS -N $jobdir/$name
#PBS -o log/$name.pbs.log
#PBS -e log/$name.pbs.err

cd \$PBS_O_WORKDIR

EOF

#SLURM --partition=$partition
#SLURM --ntasks=$nTask
#SLURM --account=$account
#SLURM --time=$ttime
#SLURM --qos=$qos
#SLURM --job-name=$jobdir/$name
#SLURM -o log/$name.slurm.log
#SLURM -e log/$name.slurm.err
#SLURM --nodes=1-1
	
	return 0
}

##
# M3C-gamess.geniso CONFIGURATION
#
function PBS.M3C-gamess.geniso()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gamess.geniso $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gamess.optg CONFIGURATION
#
function PBS.M3C-gamess.optg()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gamess.optg $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gamess.freqs CONFIGURATION
#
function PBS.M3C-gamess.freqs()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gamess.freqs $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gamess.iener CONFIGURATION
#
function PBS.M3C-gamess.iener()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gamess.iener $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.geniso CONFIGURATION
#
function PBS.M3C-gaussian.geniso()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.geniso $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.optg CONFIGURATION
#
function PBS.M3C-gaussian.optg()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.optg $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.freqs CONFIGURATION
#
function PBS.M3C-gaussian.freqs()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.freqs $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.symmetrize CONFIGURATION
#
function PBS.M3C-gaussian.symmetrize()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.symmetrize $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.iener CONFIGURATION
#
function PBS.M3C-gaussian.iener()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.iener $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-gaussian.genpot CONFIGURATION
#
function PBS.M3C-gaussian.genpot()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-gaussian.genpot $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

# ##
# # M3C-adf.geniso CONFIGURATION
# #
# function PBS.M3C-adf.geniso()
# {
# 	local queueParams=$1
# 	local name=""
# 	
# 	PBS.buildHead $queueParams $name > run$$.pbs
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.pbs << EOF
# M3C-adf.geniso $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
# EOF
# 	
# 	qsub run$$.pbs
# 	
# 	cp run$$.pbs log/
# 	rm run$$.pbs
# }

##
# M3C-adf.optg CONFIGURATION
#
function PBS.M3C-adf.optg()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-adf.optg $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C-adf.freqs CONFIGURATION
#
function PBS.M3C-adf.freqs()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C-adf.freqs $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

# ##
# # M3C-adf.symmetrize CONFIGURATION
# #
# function PBS.M3C-adf.symmetrize()
# {
# 	local queueParams=$1
# 	local name=""
# 	
# 	PBS.buildHead $queueParams $name > run$$.pbs
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.pbs << EOF
# M3C-adf.symmetrize $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
# EOF
# 	
# 	qsub run$$.pbs
# 	
# 	cp run$$.pbs log/
# 	rm run$$.pbs
# }

# ##
# # M3C-adf.iener CONFIGURATION
# #
# function PBS.M3C-adf.iener()
# {
# 	local queueParams=$1
# 	local name=""
# 	
# 	PBS.buildHead $queueParams $name > run$$.pbs
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.pbs << EOF
# M3C-adf.iener $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
# EOF
# 	
# 	qsub run$$.pbs
# 	
# 	cp run$$.pbs log/
# 	rm run$$.pbs
# }

##
# M3C.check CONFIGURATION
#
function PBS.M3C.check()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C.check $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C CONFIGURATION
#
function PBS.M3C()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

##
# M3C.p CONFIGURATION
#
function PBS.M3C.p()
{
	local queueParams=$1
	local name=""
	
	PBS.buildHead $queueParams $name > run$$.pbs
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.pbs << EOF
M3C.p $* > PBS-\$PBS_JOBID.log 2> PBS-\$PBS_JOBID.err
EOF
	
	qsub run$$.pbs
	
	cp run$$.pbs log/
	rm run$$.pbs
}

