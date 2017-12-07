#! /bin/bash

./pp_plot_csv.py output_prog_vort_t00000000*csv || exit

./pp_create_mp4.sh output_prog_vort_ output_prog_vort.mp4

