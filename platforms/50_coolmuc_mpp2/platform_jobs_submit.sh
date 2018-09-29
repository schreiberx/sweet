#! /bin/bash


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	sbatch --cluster mpp2 "$i/run.sh"
done

