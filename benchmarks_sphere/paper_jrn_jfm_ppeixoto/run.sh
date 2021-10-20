#! /usr/bin/env bash

#
# Execute all jobs by directly running them from the current shell
#


if [[ -z "$1" ]]; then
	DIRS=bench*
else
	DIRS=$@
fi

P="$(pwd)"
for JOBDIR in $DIRS; do
	cd "$JOBDIR" || exit 1
    echo $JOBDIR
    mule.benchmark.jobs_run_directly
	cd "$P"
done