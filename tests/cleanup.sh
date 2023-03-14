#! /bin/bash

source ../local_software/env_vars.sh

if [[ -z "$1" ]]; then
	JOBDIRS="$(ls -1 -d ??_* 2>/dev/null) $(ls -1 -d _??_* 2>/dev/null)"
else
	JOBDIRS=$@
fi

echo_info "Job dirs: ${JOBDIRS}"

RETDIR=$(pwd)
for i in $JOBDIRS; do
	echo_info $i
	cd "$RETDIR"
	cd "$i"

	mule.benchmark.cleanup_all || exit 1

	JOBDIRS2=$(ls -1 -d ??_* 2>/dev/null) $(ls -1 -d _??_* 2>/dev/null)

	RETDIR2="$(pwd)"
	for i2 in $JOBDIRS2; do
		echo_info $i2
		cd "$RETDIR2"
		cd "$i2"

		mule.benchmark.cleanup_all || exit 1
	done
done

echo_success_hline
echo_success "Cleanup finished"
echo_success_hline
