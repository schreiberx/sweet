#! /bin/bash


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	test -d "$i" || continue
	bsub < "$i/run.sh" || exit 2>&1
done


