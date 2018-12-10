#! /bin/bash

if [[ -z "$1" ]]; then
	echo ""
	echo "Usage:"
	echo "	$0 [procname1] [procname2] ..."
	echo ""
fi
for PROCNAME in $@; do
	PIDS=$(pgrep $@)
	$(dirname $0)/get_thread_affinities_by_pid.sh $PIDS
done
