#!/bin/bash

type=$1  # SOURCE?

export M3C_HOME="`dirname ${BASH_SOURCE[0]}`"

if [ "$type" = "SOURCE" ]
then
	export PATH=$M3C_HOME/src:$PATH
	export PATH=$M3C_HOME/utils:$PATH
else
	export PATH=$M3C_HOME/bin:$PATH
	export PATH=$M3C_HOME/utils:$PATH
fi
