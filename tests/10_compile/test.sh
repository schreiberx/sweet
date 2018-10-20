#! /bin/bash

cd "$(dirname $0)"

for i in $(ls -1 -d *); do
	echo_info_hline
	echo_info "Runing compile tests for $i"
	echo_info_hline
	$i/test.* || exit 1
done
