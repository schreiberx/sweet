#! /bin/bash


BASENAME="`pwd`"
for i in output_*; do 
	cd "$i"

	
	../../../scripts/normal_modes_compute.py ./output_normal_modes_physical_t000????????.????????.csv
	../../../scripts/normal_modes_plot_and_analyse.py ./output_normal_modes_physical_t000????????.????????.csv_evalues_complex.csv "$i"

	cd "$BASENAME"
	cp "$i/output_e_lambda.png" "a_output_e_lambda_$i.png"
	cp "$i/output_lambda.png" "a_output_lambda_$i.png"
done
