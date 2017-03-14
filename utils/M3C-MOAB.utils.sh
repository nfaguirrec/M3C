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

#############################
# M3C_GAUSSIAN_HOME
# M3C_GAUSSIAN_SCRATCH
# M3C_GAMESS_HOME
# M3C_GAMESS_SCRATCH
#############################

export M3C_SCHEDULER_NAME="MOAB"
export M3C_SCHEDULER_SUBMIT="msub"
# export M3C_SCHEDULER_JOBID="\$PBS_JOBID"
export M3C_SCHEDULER_JOBID="\$SLURM_JOB_ID"  # << Esto puede ser algo especifico de LANL, porque estan migrando todo a SLURM
export M3C_NTHREADS="" # See SCHEDULER.buildHead()

##
# BASIC CONFIGURATION
#
function SCHEDULER.buildHead()
{
	local queueParams=$1
	
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
	
	if [ -f "$HOME/.${M3C_SCHEDULER_NAME}.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $HOME/.${M3C_SCHEDULER_NAME}.properties`
	elif [ -f "$M3C_HOME/utils/${M3C_SCHEDULER_NAME}.properties" ]
	then
		values=`awk '( $1!~/^#/ ){ if( $1 == "'$partition'" ){ print $0; exit } }' $M3C_HOME/utils/${M3C_SCHEDULER_NAME}.properties`
	else
		echo "### ERROR ### ${M3C_SCHEDULER_NAME} default parameters file not found" > /dev/stderr
		return 1
	fi
	
	[ -z "$partition" ] && partition=`echo $values | awk '{print $1}'`
	[ -z "$nTask" ]     && nTask=`echo $values | awk '{print $2}'`
	[ -z "$account" ]   && account=`echo $values | awk '{print $3}'`
	[ -z "$ttime" ]     && ttime=`echo $values | awk '{print $4}'`
	[ -z "$qos" ]       && qos=`echo $values | awk '{print $5}'`

	[ -z "$qos" ] && qos=$partition
	
	local jobdir=`echo $PWD | sed s/.*$USER/~/`
	
	export M3C_NTHREADS=$nTask
	
	[ ! -d log ] && mkdir log
	
	cat << EOF
#!/bin/bash
#MSUB -q $partition
#MSUB -A $account
#MSUB -l walltime=$ttime
#MSUB -N M3C
#MSUB -o ${M3C_SCHEDULER_NAME}.log
#MSUB -e ${M3C_SCHEDULER_NAME}.err
#MSUB -l nodes=1:ppn=$nTask

EOF

###SLURM --qos=$qos   <<< There is not equivalent in MOAB
	
	return 0
}

source $M3C_HOME/utils/M3C-SCHEDULER.utils.sh

alias MOAB.M3C-gamess.geniso="SCHEDULER.M3C-gamess.geniso"
alias MOAB.M3C-gamess.optg="SCHEDULER.M3C-gamess.optg"
alias MOAB.M3C-gamess.freqs="SCHEDULER.M3C-gamess.freqs"
alias MOAB.M3C-gamess.iener="SCHEDULER.M3C-gamess.iener"
alias MOAB.M3C-gaussian.geniso="SCHEDULER.M3C-gaussian.geniso"
alias MOAB.M3C-gaussian.optg="SCHEDULER.M3C-gaussian.optg"
alias MOAB.M3C-gaussian.freqs="SCHEDULER.M3C-gaussian.freqs"
alias MOAB.M3C-gaussian.symmetrize="SCHEDULER.M3C-gaussian.symmetrize"
alias MOAB.M3C-gaussian.iener="SCHEDULER.M3C-gaussian.iener"
alias MOAB.M3C-gaussian.genpot="SCHEDULER.M3C-gaussian.genpot"
# alias MOAB.M3C-adf.geniso="SCHEDULER.M3C-adf.geniso"
alias MOAB.M3C-adf.optg="SCHEDULER.M3C-adf.optg"
alias MOAB.M3C-adf.freqs="SCHEDULER.M3C-adf.freqs"
# alias MOAB.M3C-adf.symmetrize ="SCHEDULER.M3C-adf.symmetrize"
# alias MOAB.M3C-adf.iener="SCHEDULER.M3C-adf.iener"
alias MOAB.M3C.check="SCHEDULER.M3C.check"
alias MOAB.M3C="SCHEDULER.M3C"
alias MOAB.M3C.p="SCHEDULER.M3C.p"
