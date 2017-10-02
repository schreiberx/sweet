#! /bin/bash


CSVFILE="output_prog_eta_t00000000120.00000000.csv"
OUTFILE="output_prog_eta_t00000000120.00000000.png"


for i in script_*; do
	echo "$i"
	./pp_plot_csv.py $i/$CSVFILE
	mv "$i/$OUTFILE" "$i""_""$OUTFILE"
done
