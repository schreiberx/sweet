#! /bin/bash


BASENAME="`pwd`"
for i in script_*; do 
	test -d "$i" || continue

	echo "*********"
	echo "PROCESSING $i"
	echo "*********"
	cd "$i"

	#../normal_modes_compute_exp.py ./output_normal_modes_physical_t???????????.????????.csv || exit
	../normal_modes_plot_and_analyse.py ./output_normal_modes_physical_t???????????.????????.csv_evalues_complex.csv "$i" "png" || exit
	ls

	cd "$BASENAME"
	for k in freq_high freq_high_diff; do
		cp "$i/output_$k.png" "a_output_""$k""_""$i"".png"
		cp "$i/output_$k.pdf" "a_output_""$k""_""$i"".pdf"
	done
done
