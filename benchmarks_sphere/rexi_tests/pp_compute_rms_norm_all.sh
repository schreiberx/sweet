#! /bin/bash

if [ -z "$1" ]; then
	echo "Please provide reference filename"
	exit 1
fi


FILE=`basename $1`

for i in script_*; do
	test -d "$i" || continue

	if [ -e "$i/$FILE" ]; then
		./pp_compute_rms_norm.py "$1" "$i/$FILE"
	else
		echo "$i/$FILE  [not_found]"
	fi
done
