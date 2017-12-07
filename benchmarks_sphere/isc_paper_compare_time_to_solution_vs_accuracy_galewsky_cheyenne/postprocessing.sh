#! /bin/bash

#./postprocessing_h.py > ./postprocessing_output_h.txt
#./postprocessing_potvort.py > ./postprocessing_output_potvort.txt
#./postprocessing_vort.py > ./postprocessing_output_vort.txt

for i in h potvort vort; do
	echo ./postprocessing_plot.py dt "postprocessing_output_""$i"".txt" "postprocessing_output_""$i""_err_vs_dt.pdf"
	./postprocessing_plot.py dt "postprocessing_output_""$i"".txt" "postprocessing_output_""$i""_err_vs_dt.pdf"

	echo ./postprocessing_plot.py wallclocktime "postprocessing_output_""$i"".txt" "postprocessing_output_""$i""_err_vs_wallclocktime.pdf"
	./postprocessing_plot.py wallclocktime "postprocessing_output_""$i"".txt" "postprocessing_output_""$i""_err_vs_wallclocktime.pdf"
done


#./postprocessing_output_h_err_vs_dt.py
#./postprocessing_output_h_err_vs_wallclocktime.py
#./postprocessing_output_h.py

#./postprocessing_output_potvort_err_vs_dt.py
#./postprocessing_output_potvort_err_vs_wallclocktime.py
#./postprocessing_output_potvort.py

#./postprocessing_output_vort_err_vs_dt.py
#./postprocessing_output_vort_err_vs_wallclocktime.py
#./postprocessing_output_vort.py

