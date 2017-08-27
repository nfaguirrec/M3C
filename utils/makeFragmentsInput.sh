#/bin/bash
#####################################################################################
#                                                                                   #
# This file is part of M3C project                                                  #
# Copyright (c) 2012-2016 Departamento de Química                                   #
#                         Universidad Autónoma de Madrid                            #
#                         All rights reserved.                                      #
#                                                                                   #
#                         * Néstor F. Aguirre (2012-2016)                           #
#                           nestor.aguirre@uam.es                                   #
#                         * Sergio Díaz-Tendero (2012-2016)                         #
#                           sergio.diaztendero@uam.es                               #
#                         * M. Paul-Antoine Hervieux (2012-2015)                    #
#                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  #
#                         * Manuel Alcamí (2012-2016)                               #
#                           manuel.alcami@uam.es                                    #
#                         * Fernando Martín (2012-2016)                             #
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

maxCharge=$1

[ -z "$maxCharge" ] && maxCharge=1

printf "#%19s%8s%5s\n" "XYZfile" "charge" "mult"

for q in `seq 0 $maxCharge`
do
	for f in `ls init/*.xyz`
	do
		name=`echo $f | sed 's/init\\///g'`
		mult=`molecule.minMult $f $q`
		
		printf "%20s%8d%5d\n" $name $q $mult
	done
	
	echo ""
	
	for f in `ls init/*.xyz`
	do
		name=`echo $f | sed 's/init\\///g'`
		mult=`molecule.minMult $f $q`
		
		printf "%20s%8d%5d\n" $name $q $(( $mult + 2 ))
	done
	
	echo ""
done