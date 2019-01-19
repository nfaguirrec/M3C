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
	
	# ADF configuration
	export M3C_ADF_HOME=$HOME/.adf
	export M3C_ADF_SCRATCH=/scratch/$USER/adf
	
	# NWCHEM configuration
	export M3C_NWCHEM_HOME=$HOME/.nwchem
	export M3C_NWCHEM_SCRATCH=/scratch/$USER/nwchem
	
	# LATTE configuration
	export M3C_LATTE_HOME=$HOME/bin
	export M3C_LATTE_SCRATCH=/scratch/$USER/latte
else
	export PATH=$M3C_HOME/bin:$PATH
	export PATH=$M3C_HOME/utils:$PATH
fi
