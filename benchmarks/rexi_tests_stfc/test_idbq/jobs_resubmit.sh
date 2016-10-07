#! /bin/bash

for f in `./generate_csv.py  | grep "\\.sh"`; do
	JOBID=${f/.sh/}
	JOBID=${JOBID/run_/}
	C=`bjobs -w | grep $JOBID | wc -l`

	if [ $C -eq 0 ]; then
		echo "Resubmitting $f"
		if [ -e "$f" ]; then
			bsub < "$f"
		fi
	fi
done

