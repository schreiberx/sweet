#! /bin/bash


BASENAME="`pwd`"
for i in output_*; do 
	cd "$i"

	
	../../../scripts/normal_modes_compute.py ./output_normal_modes_physical_t000????????.????????.csv
	../../../scripts/normal_modes_plot.py ./output_normal_modes_physical_t000????????.????????.csv_evalues_complex.csv "$i"

	cd "$BASENAME"
	cp "$i/output.png" "a_$i.png"
done
