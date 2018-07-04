#! /bin/bash


BASEDIR=`pwd`

for i in script_*fsph1*; do 
	test -d "$i" || continue

	echo "******************************************************"
	echo "* PROCESSING $i"
	echo "******************************************************"
	cd "$i"

	../normal_modes_compute_exp.py ./output_normal_modes_physical_t???????????.????????.csv || exit
	cd "$BASEDIR"
done


ASDF="output_normal_modes_physical_t00000000400.00000000.csv_evalues_complex.csv"


./normal_modes_relative_dispersion_errors.py \
	plot_reference	\
	a_output_combined_fsphere1_rexihalf0_rexinorm1_STDTS.pdf	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_cn_tso2_tsob1_C000400_T001_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_erk_tso1_tsob1_C000400_T001_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_erk_tso2_tsob1_C000400_T001_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_erk_tso4_tsob1_C000400_T001_M0016/$ASDF	\




./normal_modes_relative_dispersion_errors.py \
	dont_plot_reference	\
	a_output_combined_fsphere1_rexihalf0_rexinorm1_REXI.pdf	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_erk_tso4_tsob1_C000400_T001_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000128_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000256_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000512_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00001024_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00002048_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00004096_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\




./normal_modes_relative_dispersion_errors.py \
	dont_plot_reference	\
	a_output_combined_fsphere1_rexihalf0_rexinorm0_REXI.pdf	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_erk_tso4_tsob1_C000400_T001_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000128_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000256_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000512_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00001024_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00002048_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\
	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00004096_h0.15_nrm0_hlf0_bf0_ext02_M0016/$ASDF	\


#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000016_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000032_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000064_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\

#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000016_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000032_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000064_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00000128_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00008192_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF	\
#	script_g1_h100000_f0.000145842_a6371220_fsph1_u0_U0_tsm_l_rexi_tso0_tsob1_C000400_T001_REXITER_m00016384_h0.15_nrm1_hlf0_bf0_ext02_M0016/$ASDF

