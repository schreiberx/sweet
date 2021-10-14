#! /usr/bin/env bash

#
# Execute all jobs by directly running them from the current shell
#


if [[ -z "$1" ]]; then
	DIRS=job_bench*
else
	DIRS=$@
fi

P="$(pwd)"
for JOBDIR in $DIRS; do
	cd "$JOBDIR" || exit 1
	echo $JOBDIR

	if test -f "output.out";then
	    echo "Already executed"
	else
	    ./run.sh 2>"output.err" > "./output.out"

	    EXIT_CODE=$?

	    if [[ 0 -ne $EXIT_CODE ]]; then
		echo_error_hline
		echo_error_exit "ERROR"
		echo_error_hline
		cat "output.out"
		cat "output.err"
		exit 1
	    fi
	fi
	cd "$P"
done
