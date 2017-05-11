#! /bin/bash


THISDIR=`pwd`

cd "../../"

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
	#./run.sh | tee "../$i.out" || exit 2>&1
	{ { ./run.sh; } > >(tee ../$i.out ); } \
                      2> >(tee ../$i.err )

	cd ".."

done
