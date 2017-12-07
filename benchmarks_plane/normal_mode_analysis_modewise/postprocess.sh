#! /bin/bash


BASENAME="`pwd`"
for i in script_*; do 
	cd "$i"

	../normal_modes_compute.py ./output_normal_modes_physical_t000????????.????????.csv
	../normal_modes_plot_and_analyse.py ./output_normal_modes_physical_t000????????.????????.csv_evalues_complex.csv "$i"

	cd "$BASENAME"
	cp "$i/output_lambda.png" "a_lambda_$i.png"
	cp "$i/output_exp_lambda.png" "a_exp_lambda_$i.png"
done
