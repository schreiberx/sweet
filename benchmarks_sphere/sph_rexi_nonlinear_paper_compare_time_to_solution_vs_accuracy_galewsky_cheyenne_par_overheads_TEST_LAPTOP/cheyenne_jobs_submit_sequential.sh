#! /bin/bash


if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi


for i in $DIRS; do
	#
	# Wait until no job is in the queue
	#

	while true; do
		QSTATWCL=$((`qstat -n -u $USER | wc -l`))
		echo "qstat wcl: $QSTATWCL"
		test "$((QSTATWCL))" -eq "0" && break
		sleep 1
	done

	qsub "$i/run.sh"
done

