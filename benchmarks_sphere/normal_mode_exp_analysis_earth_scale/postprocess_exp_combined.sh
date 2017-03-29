#! /bin/bash


BASEDIR=`pwd`
#./normal_modes_compute_eigendecomp.py reference_modes016_bench4_nonlin0_g1_h100000_f0.00014584_a6371220_u0_fsphere0_pdeid1_robert1_tsm001_tso1_tnr001_rexim00128_rexih0.15_rexihalf1_rexiextmodes02_rexipar1_C00000010.csv
#./normal_modes_compute_eigendecomp.py reference_modes016_bench4_nonlin0_g1_h100000_f0.00014584_a6371220_u0_fsphere1_pdeid1_robert1_tsm001_tso1_tnr001_rexim00128_rexih0.15_rexihalf1_rexiextmodes02_rexipar1_C00000010.csv

for i in script_*; do 
	test -d "$i" || continue

	echo "******************************************************"
	echo "* PROCESSING $i"
	echo "******************************************************"
	cd "$i"

	#../normal_modes_compute_exp.py ./output_normal_modes_physical_t???????????.????????.csv || exit
	cd "$BASEDIR"
done


for F in 0 1; do
	./normal_modes_plot_and_analyse_combined.py \
		a_output_combined_fsphere"$F".pdf	\
		reference_modes016_bench4_nonlin0_g1_h100000_f0.00014584_a6371220_u0_fsphere"$F"_pdeid1_robert1_tsm001_tso1_tnr001_rexim00128_rexih0.15_rexihalf1_rexiextmodes02_rexipar1_C00000010.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm001_tso1_rexim00000000_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm001_tso2_rexim00000000_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm001_tso4_rexim00000000_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm003_tso2_rexim00000000_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm100_tso0_rexim00000064_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000512_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000256_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00001024_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \


	echo \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm100_tso0_rexim00000016_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm100_tso0_rexim00000032_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00000400_Tn005_t-0000001_o000.0001_tsm100_tso0_rexim00000128_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		\
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000016_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000032_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000064_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000128_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00000256_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.00014584_a6371220_u0_robert1_pdeid1_fsphere"$F"_C00002000_Tn001_t-0000001_o000.0001_tsm100_tso0_rexim00001024_rexih0.15_rexihalf1_rexiextmodes02/output_normal_modes_physical_t00000002000.00000000.csv_evalues_complex.csv

done
