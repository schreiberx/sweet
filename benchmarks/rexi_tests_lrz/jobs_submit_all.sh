#! /bin/bash

for i in run_*.sh; do
	sbatch "$i"
done

