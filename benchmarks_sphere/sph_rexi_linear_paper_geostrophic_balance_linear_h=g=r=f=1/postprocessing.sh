#! /bin/bash

RK2_FILE="script_b10_g1_h1_f1_a1_u0_U0_fsph0_tsm_l_erk_tso2_tsob1_C0.01_REXITER_m00000256_h0.15_nrm0_hlf1_bf0_ext02_M0064/output.err"

for i in 0 1; do
	for j in 0 1; do
		for ext in 0 2; do
			./pp_generate_plots_accumulated.py output "output_norm""$i""_""half""$j""_""ext""$ext".pdf ./*_tsm_l_rexi_*nrm$i*hlf$j*ext0$ext*/output.err "$RK2_FILE"
#			./pp_generate_plots_accumulated_nolog.py output "output_norm""$i""_""half""$j""_nolog.pdf" ./*_tsm100*norm$i*half$j*.err "$RK2_FILE"
		done
	done
done
