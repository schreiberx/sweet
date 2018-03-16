#! /bin/bash


THISDIR=`pwd`

cd "../../../"

source ./local_software/env_vars.sh || exit 1

cd "$THISDIR"

if [ -z "$1" ]; then
	DIRS=script_*na_sl*etdrk*
else
	DIRS=$@
fi


for i in $DIRS; do
	test -d "$i" || continue

	cd "$i"
	./run.sh | tee "../$i.out"
	#test ${PIPESTATUS[0]} -eq 0 || exit 2>&1
	cd ".."

done
