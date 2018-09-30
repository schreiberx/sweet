#! /bin/bash

REF="reference output_rk4_ts100_r1_pde1/output_prog_h_t00000129600.00000000.csv"
REF=""

#./pp_plot_csv.py $REF output_*/output_prog_h_t00000129600.00000000.csv
./pp_plot_csv.py $REF output_rexi*/output_prog_h_t00000129600.00000000.csv
./pp_plot_csv.py $REF output_rk*/output_prog_h_t00000129600.00000000.csv
./pp_plot_csv.py $REF output_lf*/output_prog_h_t00000129600.00000000.csv
./pp_plot_csv.py $REF output_cn*/output_prog_h_t00000129600.00000000.csv


for i in output_*/output_prog_h_t00000129600.00000000.png; do
	cp "$i" "result_${i//\//__}"
done
