#! /bin/bash


THISDIR=`pwd`

cd "../../../"

source ./local_software/env_vars.sh || exit 1

cd "$THISDIR"


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	test -d "$i" || continue

	cd "$i"
	./run.sh 2>"output.err" | tee "output.out"
	cd ".."

done
