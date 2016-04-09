#!/bin/bash

type=$1  # SOURCE?

export M3C_HOME="`dirname ${BASH_SOURCE[0]}`"

if [ "$type" = "SOURCE" ]
then
	export PATH=$M3C_HOME/src:$PATH
	export PATH=$M3C_HOME/utils:$PATH
	
	# GAMESS configuration
	export M3C_GAMESS_HOME=$HOME/.gamess
	export M3C_GAMESS_SCRATCH=/scratch/$USER/gamess

	# GAUSSIAN configuration
	export M3C_GAUSSIAN_HOME=$HOME/.gaussian09
	export M3C_GAUSSIAN_SCRATCH=/scratch/$USER/gaussian
else
	export PATH=$M3C_HOME/bin:$PATH
	export PATH=$M3C_HOME/utils:$PATH
fi
