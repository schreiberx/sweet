#! /bin/bash

if [ -z "$1" ]; then
	IDS=$(squeue --cluster=mpp2 | tail -n +3 | sed "s/ mpp2.*//")
	echo ${IDS}
else
	IDS=$@
fi

scancel --clusters=mpp2 $IDS
