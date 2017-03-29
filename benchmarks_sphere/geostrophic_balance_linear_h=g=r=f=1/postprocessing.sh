#! /bin/bash

RK2_FILE="script_modes064_bench10_nonlin0_g1_h1_f1_a1_u0_pdeid1_tsm1_tso2_rexim00256_rexih0.15_rexinorm0_rexihalf1_rexiextmodes02_rexipar1_C00000.01_t000000.1_o00000.01_robert1.err"

for i in 0 1; do
	for j in 0 1; do
		./pp_generate_plots_accumulated.py output "output_norm""$i""_""half""$j".pdf ./*_tsm100*norm$i*half$j*.err "$RK2_FILE"
#		./pp_generate_plots_accumulated_nolog.py output "output_norm""$i""_""half""$j""_nolog.pdf" ./*_tsm100*norm$i*half$j*.err "$RK2_FILE"
	done
done
