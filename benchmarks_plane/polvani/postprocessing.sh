#! /bin/bash

for i in script_swe_plane_polvani_*/output_diag_vort_t00000001000.00000000.csv; do

	../../scripts/pp_plot_plane_bwr_csv.py "$i"

	j="${i/.csv/.png}"

	mv "$j" "${j/\//_}"

done



