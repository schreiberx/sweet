#! /bin/bash


if [ -z "$1" ]; then
	echo "Usage:"
	echo "	$0 [pid]"
	echo
	exit 1
fi

TIDS=$(ps --no-headers -mo tid -p $1 | tail -n +1)
TNR=0
for tid in $TIDS; do
	# detect '-' symbol
	test $((${tid}1+1)) -eq 0 && continue

	MASK=$(taskset -p $tid | sed "s/.*mask: //")
	echo "Thread ${TNR} using mask ${MASK}"
	TNR=$((TNR+1))
done
