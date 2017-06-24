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


ASDF="output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv"


./normal_modes_relative_dispersion_errors.py \
	plot_reference	\
	a_output_combined_fsphere1_rexihalf1_rexinorm0_rexiextmodes02_STDTS.pdf	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso1_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso2_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso4_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_cn_tso2_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\



./normal_modes_relative_dispersion_errors.py \
	dont_plot_reference	\
	a_output_combined_fsphere1_rexihalf1_rexinorm0_rexiextmodes02_REXI.pdf	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_erk_tso4_rexim00000000_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00000256_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00000512_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00001024_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00002048_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00004096_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00008192_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_u0_robert1_pdeid0_fsphere1_C00000400_Tn001_t00000400_o000.0001_tsm_l_rexi_tso0_rexim00016384_rexih0.15_rexinorm1_rexihalf1_rexiprealloc0_rexiextmodes02/$ASDF
