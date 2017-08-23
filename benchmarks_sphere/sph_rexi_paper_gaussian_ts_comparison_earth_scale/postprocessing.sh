#! /bin/bash

REFFILE=output_prog_h_t00000129600.00000000.csv
REF="reference output_rk4_ts100/$REFFILE"
#REF=""

if true; then
	./pp_plot_csv.py $REF output_*/$REFFILE

	for i in output_*/${REFFILE/.csv/.png}; do
		cp "$i" "result_${i//\//__}"
	done
fi

if true; then
	./pp_plot_csv_pdf.py $REF output_*/$REFFILE

	for i in output_*/${REFFILE/.csv/.pdf}; do
		cp "$i" "result_${i//\//__}"
	done
fi
