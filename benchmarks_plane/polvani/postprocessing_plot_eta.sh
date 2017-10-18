#! /bin/bash


#CSVFILE="output_prog_eta_t00000000120.00000000.csv"
#OUTFILE="output_prog_eta_t00000000120.00000000.png"
#FILEA="output_prog_eta_"

CSVFILE="output_prog_vort_t00000001000.00000000.csv"
OUTFILE="output_prog_vort_t00000001000.00000000.png"
FILEA="output_prog_eta_"

FILEEXT="png"


DIRREF="script_ln2_ref_b100_g9.81_h10000_f7.2921e-05_p0_a6371220_u0.0_rob1_fsph0_tsm_ln_erk_tso4_tsob4_REXICI_n00000064_sx50.0_sy50.0_mu0.0_prcircle_nrm0_hlf1_pre1_ext00_C00000040_M0128"

for i in script_*; do
	echo "$i"
	./pp_plot_csv.py "$i/$CSVFILE" "$DIRREF/$CSVFILE"
	mv "$i/$OUTFILE" "$FILEA""_""$i"".""$FILEEXT"
done
