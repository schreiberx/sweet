#! /bin/bash


if [[ -z "$1" ]]; then
	echo "Usage:"
	echo "	$0 [pid]"
	echo
	exit 1
fi

PID=$1

getmask()
(
	OUTPUT=$(taskset -p $1 2>&1)
	test $? -ne 0 && echo_error_exit "${OUTPUT}"

	OUTPUT="${OUTPUT/*mask: /}"
	echo -n "$OUTPUT"
)

MASK=$(eval getmask $PID)
test $? -ne 0 && { echo "$MASK"; exit 1; }

echo "Master process using mask ${MASK}"

TIDS=$(ps --no-headers -mo tid -p $PID | tail -n +1 || exit 1)
TNR=0
for TID in $TIDS; do
	# detect '-' symbol
	test "#$TID" = "#-" && continue

	MASK=$(eval getmask $TID)
	test $? -ne 0 && { echo "$MASK"; exit 1; }

	echo "Thread ${TNR} using mask ${MASK}"
	TNR=$((TNR+1))
done
