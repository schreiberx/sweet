#! /bin/bash

echo "******************************************"
echo "* PREPROCESS"
echo "******************************************"

for i in h potvort vort; do

	OUTFILE="./output_benchmark_${i}.txt"
	if [ ! -f $OUTFILE ]; then
		echo "Extracting data and storing it to '${OUTFILE}'"
		EXEC="./postprocessing_script_output_2_output_benchmark.py $i"
		echo "****************************************************"
		echo "* EXEC: $EXEC"
		echo "****************************************************"
		$EXEC > $OUTFILE || exit 1
	else
		echo "Skipping creation of '${OUTFILE}' (remove file to recompute data)"
	fi
done

echo "******************************************"
echo "* PLOTTING"
echo "******************************************"

for i in h potvort vort; do
	EXEC="./postprocessing_output_benchmark_2_plot.py dt output_benchmark_${i}.txt output_benchmark_${i}_err_vs_dt.pdf"
	echo "$EXEC"
	$EXEC

	EXEC="./postprocessing_output_benchmark_2_plot.py wallclocktime output_benchmark_${i}.txt output_benchmark_${i}_err_vs_wallclocktime.pdf"
	echo "****************************************************"
	echo "* EXEC: $EXEC"
	echo "****************************************************"
	$EXEC
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

