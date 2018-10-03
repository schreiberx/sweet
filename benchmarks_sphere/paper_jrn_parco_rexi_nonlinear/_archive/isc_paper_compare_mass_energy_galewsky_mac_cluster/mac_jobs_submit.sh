#! /bin/bash


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	sbatch "$i/run.sh"
done

