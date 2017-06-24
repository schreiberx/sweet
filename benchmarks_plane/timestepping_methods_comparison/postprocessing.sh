#! /bin/bash


DATAFILE="output_prog_h_t00000000000.10000000.csv"

for group in l1 l2 ln1 ln2; do
	REFDIR=`ls -d1 "script_""$group""_ref"*`
	echo
	echo "Using reference: $REFDIR"

	for DIR in "script_""$group"*; do
		test -d "$DIR" || continue
		./pp_compute_max_and_rms_errors.py "$REFDIR/$DATAFILE" "$DIR/$DATAFILE" "$DIR"
	done
done

