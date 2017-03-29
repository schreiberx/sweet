#! /bin/bash

RK2_FILE="script_modes064_bench11_nonlin0_g-1_h-1_f-1_a-1_u0_pdeid1_tsm1_tso2_rexim00000_rexih0.15_rexinorm0_rexihalf0_rexiextmodes00_rexipar1_C00000001_t00000010_o00000001_robert1.err"

for i in 0 1; do
	for j in 0 1; do
		./pp_generate_plots_accumulated.py output "output_norm""$i""_""half""$j".pdf ./*_tsm100*norm$i*half$j*.err "$RK2_FILE"
		./pp_generate_plots_accumulated_nolog.py output "output_norm""$i""_""half""$j""_nolog.pdf" ./*_tsm100*norm$i*half$j*.err "$RK2_FILE"
	done
done
