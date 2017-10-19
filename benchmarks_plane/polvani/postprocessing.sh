#! /bin/bash

for i in script_swe_plane_rexi_polvani_*/output_prog_vort_t00000000900.00000000.csv; do

	../../pp_plot_quad_bwr_csv.py "$i"

	j="${i/.csv/.png}"

	mv "$j" "${j/\//_}"

done



