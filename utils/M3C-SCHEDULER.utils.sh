#!/bin/bash
##################################################################
#
#  This file is part of M3C
#  Copyright (C) by authors (2017-2017)
#  
#  Authors:
#    * Dr. Néstor F. Aguirre (2017-2017)
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

##
# M3C-gamess.geniso CONFIGURATION
#
function SCHEDULER.M3C-gamess.geniso()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gamess.geniso $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gamess.optg CONFIGURATION
#
function SCHEDULER.M3C-gamess.optg()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gamess.optg $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gamess.freqs CONFIGURATION
#
function SCHEDULER.M3C-gamess.freqs()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gamess.freqs $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gamess.iener CONFIGURATION
#
function SCHEDULER.M3C-gamess.iener()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gamess.iener $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.geniso CONFIGURATION
#
function SCHEDULER.M3C-gaussian.geniso()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.geniso $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.optg CONFIGURATION
#
function SCHEDULER.M3C-gaussian.optg()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.optg $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.freqs CONFIGURATION
#
function SCHEDULER.M3C-gaussian.freqs()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.freqs $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.symmetrize CONFIGURATION
#
function SCHEDULER.M3C-gaussian.symmetrize()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.symmetrize $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.iener CONFIGURATION
#
function SCHEDULER.M3C-gaussian.iener()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.iener $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-gaussian.genpot CONFIGURATION
#
function SCHEDULER.M3C-gaussian.genpot()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-gaussian.genpot $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

# ##
# # M3C-adf.geniso CONFIGURATION
# #
# function SCHEDULER.M3C-adf.geniso()
# {
# 	local queueParams=$1
# 	
# 	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
# M3C-adf.geniso $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
# EOF
# 	
# 	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
# 	
# 	cp run$$.${M3C_SCHEDULER_NAME} log/
# 	rm run$$.${M3C_SCHEDULER_NAME}
# }

##
# M3C-adf.optg CONFIGURATION
#
function SCHEDULER.M3C-adf.optg()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-adf.optg $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C-adf.freqs CONFIGURATION
#
function SCHEDULER.M3C-adf.freqs()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C-adf.freqs $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

# ##
# # M3C-adf.symmetrize CONFIGURATION
# #
# function SCHEDULER.M3C-adf.symmetrize()
# {
# 	local queueParams=$1
# 	
# 	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
# M3C-adf.symmetrize $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
# EOF
# 	
# 	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
# 	
# 	cp run$$.${M3C_SCHEDULER_NAME} log/
# 	rm run$$.${M3C_SCHEDULER_NAME}
# }

# ##
# # M3C-adf.iener CONFIGURATION
# #
# function SCHEDULER.M3C-adf.iener()
# {
# 	local queueParams=$1
# 	
# 	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
# 	[ "$?" -eq 1 ] && return 0
# 	
# 	shift # $1 will be discarded
# 	
# 	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
# M3C-adf.iener $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
# EOF
# 	
# 	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
# 	
# 	cp run$$.${M3C_SCHEDULER_NAME} log/
# 	rm run$$.${M3C_SCHEDULER_NAME}
# }

##
# M3C.check CONFIGURATION
#
function SCHEDULER.M3C.check()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C.check $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C CONFIGURATION
#
function SCHEDULER.M3C()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

##
# M3C.p CONFIGURATION
#
function SCHEDULER.M3C.p()
{
	local queueParams=$1
	
	SCHEDULER.buildHead $queueParams > run$$.${M3C_SCHEDULER_NAME}
	[ "$?" -eq 1 ] && return 0
	
	shift # $1 will be discarded
	
	cat >> run$$.${M3C_SCHEDULER_NAME} << EOF
M3C.p $* > ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.log 2> ${M3C_SCHEDULER_NAME}-${M3C_SCHEDULER_JOBID}.err
EOF
	
	$M3C_SCHEDULER_SUBMIT run$$.${M3C_SCHEDULER_NAME}
	
	cp run$$.${M3C_SCHEDULER_NAME} log/
	rm run$$.${M3C_SCHEDULER_NAME}
}

