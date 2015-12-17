#! /bin/bash

for i in output_*h.csv; do
	./generate_heatmap_plot.py $i blarg_${i/csv/png} -10 10;
done

