#!/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
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
# M3C_GAUSSIAN_HOME
# M3C_GAUSSIAN_SCRATCH
# M3C_GAMESS_HOME
# M3C_GAMESS_SCRATCH
#############################

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
#SBATCH --nodes=1-1

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
M3C-gamess.geniso $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C-gamess.optg $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C-gamess.freqs $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gamess.iener CONFIGURATION
#
function SLURM.M3C-gamess.iener()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C-gamess.iener $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C-gaussian.geniso $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C-gaussian.optg $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C-gaussian.freqs $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.symmetrize CONFIGURATION
#
function SLURM.M3C-gaussian.symmetrize()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C-gaussian.symmetrize $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.iener CONFIGURATION
#
function SLURM.M3C-gaussian.iener()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C-gaussian.iener $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

##
# M3C-gaussian.genpot CONFIGURATION
#
function SLURM.M3C-gaussian.genpot()
{
	local queueParams=$1
	local name=""
	
	SLURM.buildHead $queueParams $name > run$$.slurm
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.slurm << EOF
M3C-gaussian.genpot $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
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
M3C.p $* > SLURM-\$SLURM_JOB_ID.log 2> SLURM-\$SLURM_JOB_ID.err
EOF
	
	sbatch run$$.slurm
	
	cp run$$.slurm log/
	rm run$$.slurm
}

