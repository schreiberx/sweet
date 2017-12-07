#! /bin/bash

./postprocessing_h.py > ./postprocessing_output_h.txt
./postprocessing_output_h_err_vs_dt.py
./postprocessing_output_h_err_vs_simtime.py
./postprocessing_output_h.py

./postprocessing_vort.py > ./postprocessing_output_vort.txt
./postprocessing_output_vort_err_vs_dt.py
./postprocessing_output_vort_err_vs_simtime.py
./postprocessing_output_vort.py

