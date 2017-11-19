#! /bin/bash

REFFILE="output_prog_h_t00000129600.00000000.csv"
REFSCRIPT="script_g9.80616_h10000_f7.292e-05_a6371220_u0_U0_fsph0_tsm_l_erk_tso4_tsob1_C0050_REXITER_m00000256_h0.15_nrm1_hlf0_bf0_ext00_M0128"
REF="reference $REFSCRIPT/$REFFILE"

#if false; then
if true; then
	./pp_plot_csv.py $REF script_*/$REFFILE

	for i in script_*/${REFFILE/.csv/.png}; do
		cp "$i" "result_${i//\//__}"
	done
fi

#if false; then
if true; then
	./pp_plot_csv_pdf.py $REF script_*/$REFFILE

	for i in script_*/${REFFILE/.csv/.pdf}; do
		cp "$i" "result_${i//\//__}"
	done
fi
