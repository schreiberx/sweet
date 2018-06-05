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
    j="$i"".out"
    if [ ! -f "$j" ]; then
	cd "$i"
	echo $i
	./run.sh | tee "../$i.out"
	#test ${PIPESTATUS[0]} -eq 0 || exit 2>&1
	cd ".."
    fi

done
