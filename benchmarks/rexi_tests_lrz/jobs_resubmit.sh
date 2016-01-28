#! /bin/bash

for f in `./generate_csv.py  | grep "^run.*\\.sh"`; do
	JOBNAME=${f/.sh/}
	JOBNAME=${JOBNAME/run_/}
	C=`squeue -M mpp2 -n $JOBNAME | wc -l`

	if [ $C -eq 2 ]; then
		echo "Resubmitting $f"
		if [ -e "$f" ]; then
			echo sbatch "$f"
			sbatch "$f"
		fi
	fi
done

