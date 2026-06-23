#!/bin/bash

build="$1"
if [ -z "$build" ]
then
    echo "### ERROR ### Build directory required as first argument (e.g., build_gfortran)"
    return 1
fi

export M3C_HOME="$(dirname "$(realpath "${BASH_SOURCE[0]}")")/$build"
export PATH=$PATH:$M3C_HOME/src
export PATH=$PATH:$M3C_HOME/utils

# GAMESS configuration
export M3C_GAMESS_HOME=$HOME/.gamess
export M3C_GAMESS_SCRATCH=/scratch/$USER/gamess

# GAUSSIAN configuration
export M3C_GAUSSIAN_HOME=$HOME/.gaussian09
export M3C_GAUSSIAN_SCRATCH=/scratch/$USER/gaussian

return 0
