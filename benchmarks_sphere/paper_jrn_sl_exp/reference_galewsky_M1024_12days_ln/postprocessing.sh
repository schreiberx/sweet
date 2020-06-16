#! /usr/bin/env bash

JOB_DIR=job_benchref_RT_u0.0_tsm_ln_erk_tso4_tsob4_dt00010.00_W-00001

for i in $JOB_DIR/output_prog_phi_pert_*.sweet; do

	FILE_PHI="$i"
	FILE_VRT="${i/phi_pert/vrt}"
	FILE_DIV="${i/phi_pert/div}"

	./postprocessing_plot_solution.py "$JOB_DIR" "$FILE_PHI" "$FILE_VRT" "$FILE_DIV"
done


for i in kinetic_energy output_prog_h_pert output_prog_vrt output_prog_div; do
	pdftk $JOB_DIR/plot_$i*.pdf output ./output_$i.pdf
done
