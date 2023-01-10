#! /bin/bash

cd "$(dirname $0)"

for i in $(ls -1 -d ??_*/); do
	i=$(basename "$i")
	echo_info_hline
	echo_info "Running compile tests for $i"
	echo_info_hline
	OUTFILE=output_"${i/\//}"".out"
	./$i/test.* > $OUTFILE 2>&1 || { cat "$OUTFILE"; exit 1; }
done
