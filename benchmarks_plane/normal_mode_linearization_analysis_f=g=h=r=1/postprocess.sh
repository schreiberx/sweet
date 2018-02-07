#! /bin/bash


BASENAME="`pwd`"

if [ -z "$1" ]; then
	DIRS=script_*
else
	DIRS=$@
fi

for i in $DIRS; do 
	./normal_modes_compute.py  $i/output_normal_modes_physical_t00000000000.00001000.csv
	./normal_modes_plot_and_analyse.py $i/output_normal_modes_physical_t00000000000.00001000.csv_evalues_complex.csv

#	test -d "$i" || continue
#	cd "$i"

#	cp "$i/output.png" "a_output_""$I"".png"
#	cp "$i/output.pdf" "a_output_""$I"".pdf"
done

