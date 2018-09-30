#! /bin/bash


echo "Extracing H"
./postprocessing_h.py > ./output_h.txt

echo "Extracing PotVort"
./postprocessing_potvort.py > ./output_potvort.txt

echo "Extracing Vort"
./postprocessing_vort.py > ./output_vort.txt

for i in h potvort vort; do
	echo ./postprocessing_plot.py dt "output_""$i"".txt" "output_""$i""_err_vs_dt.pdf"
	./postprocessing_plot.py dt "output_""$i"".txt" "output_""$i""_err_vs_dt.pdf"

	echo ./postprocessing_plot.py wallclocktime "output_""$i"".txt" "output_""$i""_err_vs_wallclocktime.pdf"
	./postprocessing_plot.py wallclocktime "output_""$i"".txt" "output_""$i""_err_vs_wallclocktime.pdf"
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

