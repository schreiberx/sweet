#! /usr/bin/env bash

cd "$(dirname $(realpath $0))"

JOB_DIR=job_benchref_RT_u0.0_tsm_ln_erk_tso4_tsob4_dt00010.00_W-00001

for i in $JOB_DIR/output_prog_phi_pert_*.sweet; do

	FILE_PHI="$i"
	FILE_VRT="${i/phi_pert/vrt}"
	FILE_DIV="${i/phi_pert/div}"

	../postprocessing_generate_spectrum_plots.py "$JOB_DIR" "$FILE_PHI" "$FILE_VRT" "$FILE_DIV" || exit 1
done


for i in spectrum_kinetic_energy spectrum_phi_pert spectrum_vrt spectrum_div output_prog_h_pert output_prog_vrt output_prog_div; do
	pdftk $JOB_DIR/plot_$i*.pdf output ./output_$i.pdf || exit 1
done
