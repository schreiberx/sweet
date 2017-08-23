#! /bin/bash


BASEDIR=`pwd`

for i in script_*; do 
	test -d "$i" || continue

	echo "******************************************************"
	echo "* PROCESSING $i"
	echo "******************************************************"
	cd "$i"

	#../normal_modes_compute_exp.py ./output_normal_modes_physical_t???????????.????????.csv || exit
	cd "$BASEDIR"
done


for F in 1; do
	./normal_modes_plot_and_analyse_combined.py \
		a_output_combined_fsphere"$F".pdf	\
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso1_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso2_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso4_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_cn_tso2_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00000064_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00000128_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00000512_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00001024_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00002048_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00004096_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00008192_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \
		script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere"$F"_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00016384_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv \

done
