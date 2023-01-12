#! /bin/bash

cd "$(dirname $0)"

for i in $(ls -1 -d ??_*/ | sort); do
	i=$(basename "$i")
	if [[ $i == *_no_test ]]; then
		continue
	fi
	echo_info_hline
	echo_info "Running compile tests for $i"
	echo_info_hline

	PWD_BACKUP="$(pwd)"

	cd "$i"
	OUTFILE="output_${i/\//}.out"

	if [[ "$?" != "0" ]]; then
		cat "$OUTFILE"
		exit 1
	fi

	# Cleanup subbenchmark
	mule.benchmark.cleanup_all

	cd "$PWD_BACKUP"

done
