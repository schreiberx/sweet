#! /bin/bash

for i in `squeue -M mpp2 -u di69fol -o "%i" | tail -n +2`; do
	echo "CANCEL $i"
	scancel -M mpp2 $i;
done
