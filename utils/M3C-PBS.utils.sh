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
# M3C_GAMESS_HOME
# M3C_GAMESS_SCRATCH
# M3C_GAUSSIAN_HOME
# M3C_GAUSSIAN_SCRATCH
# M3C_ADF_HOME
# M3C_ADF_SCRATCH
# M3C_NWCHEM_HOME
# M3C_NWCHEM_SCRATCH
# M3C_LATTE_HOME
# M3C_LATTE_SCRATCH
# M3C_LAMMPS_HOME
# M3C_LAMMPS_SCRATCH
#############################

export M3C_SCHEDULER_NAME="PBS"
export M3C_SCHEDULER_SUBMIT="qsub"
export M3C_SCHEDULER_JOBID="\$PBS_JOBID"
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
	
	cat << EOF
#!/bin/bash
#PBS -q $partition
#PBS -l nodes=1:ppn=$nTask,walltime=$ttime
###PBS --account=$account # There is not equivalence for QOS in PBS
###PBS --qos=$qos # There is not equivalence for QOS in PBS
#PBS -N $jobdir/$name
#PBS -o ${M3C_SCHEDULER_NAME}.log
#PBS -e ${M3C_SCHEDULER_NAME}.err

cd \$PBS_O_WORKDIR

EOF

	return 0
}

source $M3C_HOME/utils/M3C-SCHEDULER.utils.sh

alias PBS.M3C-gamess.geniso="SCHEDULER.M3C-gamess.geniso"
alias PBS.M3C-gamess.optg="SCHEDULER.M3C-gamess.optg"
alias PBS.M3C-gamess.freqs="SCHEDULER.M3C-gamess.freqs"
alias PBS.M3C-gamess.iener="SCHEDULER.M3C-gamess.iener"

alias PBS.M3C-gaussian.geniso="SCHEDULER.M3C-gaussian.geniso"
alias PBS.M3C-gaussian.optg="SCHEDULER.M3C-gaussian.optg"
alias PBS.M3C-gaussian.freqs="SCHEDULER.M3C-gaussian.freqs"
alias PBS.M3C-gaussian.symmetrize="SCHEDULER.M3C-gaussian.symmetrize"
alias PBS.M3C-gaussian.iener="SCHEDULER.M3C-gaussian.iener"
alias PBS.M3C-gaussian.genpot="SCHEDULER.M3C-gaussian.genpot"

# alias PBS.M3C-adf.geniso="SCHEDULER.M3C-adf.geniso"
alias PBS.M3C-adf.optg="SCHEDULER.M3C-adf.optg"
alias PBS.M3C-adf.freqs="SCHEDULER.M3C-adf.freqs"
# alias PBS.M3C-adf.symmetrize ="SCHEDULER.M3C-adf.symmetrize"
# alias PBS.M3C-adf.iener="SCHEDULER.M3C-adf.iener"

# alias PBS.M3C-nwchem.geniso="SCHEDULER.M3C-nwchem.geniso"
alias PBS.M3C-nwchem.optg="SCHEDULER.M3C-nwchem.optg"
alias PBS.M3C-nwchem.freqs="SCHEDULER.M3C-nwchem.freqs"
# alias PBS.M3C-nwchem.symmetrize="SCHEDULER.M3C-nwchem.symmetrize"
# alias PBS.M3C-nwchem.iener="SCHEDULER.M3C-nwchem.iener"

alias PBS.M3C.check="SCHEDULER.M3C.check"
alias PBS.M3C="SCHEDULER.M3C"
alias PBS.M3C.p="SCHEDULER.M3C.p"
