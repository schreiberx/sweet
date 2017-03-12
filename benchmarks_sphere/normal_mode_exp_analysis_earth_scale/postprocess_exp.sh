#! /bin/bash


BASENAME="`pwd`"
for i in script_*; do 
	test -d "$i" || continue
	cd "$i"

	../normal_modes_compute_exp.py ./output_normal_modes_physical_t???????????.????????.csv
	../normal_modes_plot_and_analyse.py ./output_normal_modes_physical_t???????????.????????.csv_evalues_complex.csv "$i"

	cd "$BASENAME"
	for k in lambda exp_lambda freq_low freq_high all_evalues; do
		cp "$i/output_$k.png" "a_output_""$k""_""$i"".png"
	done
done
