#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
#                                                                                   #
# Copyright (C) by authors (2017-2018)                                              #
#                                                                                   #
#    Authors:                                                                       #
#       * Néstor F. Aguirre (2017-2018)                                             #
#         nfaguirrec@gmail.com                                                      #
#                                                                                   #
# Copyright (c) 2013-2016 Departamento de Química                                   #
#                         Universidad Autónoma de Madrid                            #
#                         All rights reserved.                                      #
#                                                                                   #
#                         * Néstor F. Aguirre (2013-2016)                           #
#                           nestor.aguirre@uam.es                                   #
#                         * Sergio Díaz-Tendero (2013-2016)                         #
#                           sergio.diaztendero@uam.es                               #
#                         * M. Paul-Antoine Hervieux (2013-2015)                    #
#                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  #
#                         * Manuel Alcamí (2013-2016)                               #
#                           manuel.alcami@uam.es                                    #
#                         * Fernando Martín (2013-2016)                             #
#                           fernando.martin@uam.es                                  #
#                                                                                   #
#  Redistribution and use in source and binary forms, with or without               #
#  modification, are permitted provided that the following conditions are met:      #
#                                                                                   #
#  1. Redistributions of source code must retain the above copyright notice, this   #
#     list of conditions and the following disclaimer.                              #
#  2. Redistributions in binary form must reproduce the above copyright notice,     #
#     this list of conditions and the following disclaimer in the documentation     #
#     and/or other materials provided with the distribution.                        #
#  3. Neither the name of the copyright holders nor the names of its contributors   #
#     may be used to endorse or promote products derived from this software         #
#     without specific prior written permission.                                    #
#                                                                                   #
#  The copyright holders provide no reassurances that the source code provided      #
#  does not infringe any patent, copyright, or any other intellectual property      #
#  rights of third parties.  The copyright holders disclaim any liability to any    #
#  recipient for claims brought against recipient by any third party for            #
#  infringement of that parties intellectual property rights.                       #
#                                                                                   #
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  #
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    #
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           #
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR  #
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     #
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      #
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       #
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    #
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     #
#                                                                                   #
#####################################################################################

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

export M3C_SCHEDULER_NAME="SLURM"
export M3C_SCHEDULER_SUBMIT="sbatch"
export M3C_SCHEDULER_JOBID="\$SLURM_JOBID"
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
##SBATCH --partition=$partition
##SBATCH --account=$account
##SBATCH --time=$ttime
##SBATCH --qos=$qos
#SBATCH --job-name=$jobdir/$name
#SBATCH -o ${M3C_SCHEDULER_NAME}.log
#SBATCH -e ${M3C_SCHEDULER_NAME}.err
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=$nTask

export SLURM_NTASKS=\$(( \$SLURM_NNODES*\$SLURM_CPUS_ON_NODE/\$SLURM_CPUS_PER_TASK ))

EOF
	
	return 0
}

source $M3C_HOME/utils/M3C-SCHEDULER.utils.sh

alias SLURM.M3C-gamess.geniso="SCHEDULER.M3C-gamess.geniso"
alias SLURM.M3C-gamess.optg="SCHEDULER.M3C-gamess.optg"
alias SLURM.M3C-gamess.freqs="SCHEDULER.M3C-gamess.freqs"
alias SLURM.M3C-gamess.iener="SCHEDULER.M3C-gamess.iener"

alias SLURM.M3C-gaussian.geniso="SCHEDULER.M3C-gaussian.geniso"
alias SLURM.M3C-gaussian.optg="SCHEDULER.M3C-gaussian.optg"
alias SLURM.M3C-gaussian.freqs="SCHEDULER.M3C-gaussian.freqs"
alias SLURM.M3C-gaussian.symmetrize="SCHEDULER.M3C-gaussian.symmetrize"
alias SLURM.M3C-gaussian.iener="SCHEDULER.M3C-gaussian.iener"
alias SLURM.M3C-gaussian.genpot="SCHEDULER.M3C-gaussian.genpot"

# alias SLURM.M3C-adf.geniso="SCHEDULER.M3C-adf.geniso"
alias SLURM.M3C-adf.optg="SCHEDULER.M3C-adf.optg"
alias SLURM.M3C-adf.freqs="SCHEDULER.M3C-adf.freqs"
# alias SLURM.M3C-adf.symmetrize ="SCHEDULER.M3C-adf.symmetrize"
# alias SLURM.M3C-adf.iener="SCHEDULER.M3C-adf.iener"

# alias SLURM.M3C-nwchem.geniso="SCHEDULER.M3C-nwchem.geniso"
alias SLURM.M3C-nwchem.optg="SCHEDULER.M3C-nwchem.optg"
alias SLURM.M3C-nwchem.freqs="SCHEDULER.M3C-nwchem.freqs"
# alias SLURM.M3C-nwchem.symmetrize="SCHEDULER.M3C-nwchem.symmetrize"
# alias SLURM.M3C-nwchem.iener="SCHEDULER.M3C-nwchem.iener"

# alias SLURM.M3C-latte.geniso="SCHEDULER.M3C-latte.geniso"
alias SLURM.M3C-latte.optg="SCHEDULER.M3C-latte.optg"
alias SLURM.M3C-latte.freqs="SCHEDULER.M3C-latte.freqs"
# alias SLURM.M3C-latte.symmetrize="SCHEDULER.M3C-latte.symmetrize"
# alias SLURM.M3C-latte.iener="SCHEDULER.M3C-latte.iener"

# alias SLURM.M3C-lammps.geniso="SCHEDULER.M3C-lammps.geniso"
alias SLURM.M3C-lammps.optg="SCHEDULER.M3C-lammps.optg"
# alias SLURM.M3C-lammps.freqs="SCHEDULER.M3C-lammps.freqs"
# alias SLURM.M3C-lammps.symmetrize="SCHEDULER.M3C-lammps.symmetrize"
# alias SLURM.M3C-lammps.iener="SCHEDULER.M3C-lammps.iener"

alias SLURM.M3C.check="SCHEDULER.M3C.check"
alias SLURM.M3C="SCHEDULER.M3C"
alias SLURM.M3C.p="SCHEDULER.M3C.p"
