#! /bin/bash


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	qsub "$i/run.sh"
done

