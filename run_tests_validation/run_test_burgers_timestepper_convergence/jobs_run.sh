#! /bin/bash


THISDIR=`pwd`

cd "../../"

#source ./local_software/env_vars.sh || exit 1
source ~/.bashrc
sweetenv

cd "$THISDIR"


if [ -z "$1" ]; then
	DIRS=script_*burgers*
else
	DIRS=$@
fi


for i in $DIRS; do
	test -d "$i" || continue

	cd "$i"
	echo "$i"
	./run.sh 2>&1 > "../$i.out" || exit 1
	cd ".."

done
