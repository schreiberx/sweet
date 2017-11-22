#! /bin/bash


#./postprocessing_h_compare_ref_output_text.sh
#./postprocessing_vort_compare_ref_output_text.sh

echo "It's wise to first filter out the just generated .tex files and then to progress with the plotting"

./postprocessing_h_plot_err_vs_dt.py
./postprocessing_h_plot_err_vs_simtime.py
