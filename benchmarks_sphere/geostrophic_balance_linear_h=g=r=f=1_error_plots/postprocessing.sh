#! /bin/bash

F0=output_ref_diff_h_t00000000000.00000000.csv
F=output_ref_diff_h_t???????????.????????.csv

for i in script_modes064_bench10_nonlin0*robert1; do
	echo "*************************************"
	echo "$i"
	echo "*************************************"

	rm -f "$i/$F0"
	FILE=$(ls $i/$F)
	FILE_OUT1=${FILE/.csv/.pdf}
	FILE_OUT2="output_${FILE_OUT1/\//_}"

	./pp_plot_csv_pdf.py $FILE

	mv "$FILE_OUT1" "$FILE_OUT2"
done

