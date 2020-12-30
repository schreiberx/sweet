#! /bin/bash

if [[ -z "$1" ]]; then
	echo ""
	echo "Usage:"
	echo "	$0 [procname1] [procname2] ..."
	echo ""
fi
for PROCNAME in $@; do
	PIDS=$(pgrep $@ || exit 1)
	$(dirname $0)/get_thread_affinities_by_pid.sh $PIDS || exit 1
done
