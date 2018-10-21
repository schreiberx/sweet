#! /bin/bash

cd "$(basename $0)"

source ../local_software/env_vars.sh

if [[ -z "$1" ]]; then
	JOBDIRS=$(ls -1 -d ??_*)
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
done

echo_success_hline
echo_success "Cleanup finished"
echo_success_hline
