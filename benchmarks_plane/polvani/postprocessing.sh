#! /bin/bash

for i in script_swe_plane_rexi_polvani_*01/output_prog_vort_t00000001000.00000000.csv; do
	../../scripts/pp_plot_plane_csv.py "$i"
	mv "$i" "${i/\//_}"
done



