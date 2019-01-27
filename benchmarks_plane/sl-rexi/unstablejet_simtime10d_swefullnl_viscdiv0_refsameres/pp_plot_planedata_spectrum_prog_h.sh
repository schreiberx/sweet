#! /bin/bash

if true; then

	# Which job directories to include in the data?
	DIRS=""
	DIRS+=" job_bench_RT_tsm_l_cn_na_sl_nd_settls_tso2_tsob2_C1.406e+01_S864000"
	DIRS+=" job_bench_RT_tsm_ln_erk_tso2_tsob2_C1.406e+01_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_na_sl_nd_etdrk_tso2_tsob2_C1.406e+01_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C1.406e+01_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_n_erk_tso2_tsob2_C1.406e+01_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_n_etdrk_tso2_tsob2_C1.406e+01_S864000"

	DIRS+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000"

	# Which tag in pickled job data to use to generate spectrum plots?
	SPECTAG="output_prog_h_pert_t00000864000.00000000_spectrum"

	# Output plot file
	OUTPUTFILE="plot_prog_h_pert_spectrum_simtime_10days_C1.406e+01.pdf"

	# Plot!
	../pp_plot_planedata_spectrum.py $SPECTAG $OUTPUTFILE $DIRS

fi


if true; then

	# Which job directories to include in the data?
	DIRS=""
	DIRS+=" job_bench_RT_tsm_l_cn_na_sl_nd_settls_tso2_tsob2_C1.125e+02_S864000"
# unstable
#	DIRS+=" job_bench_RT_tsm_ln_erk_tso2_tsob2_C1.125e+02_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_na_sl_nd_etdrk_tso2_tsob2_C1.125e+02_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C1.125e+02_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_n_erk_tso2_tsob2_C1.125e+02_S864000"
	DIRS+=" job_bench_RT_tsm_l_rexi_n_etdrk_tso2_tsob2_C1.125e+02_S864000"

	DIRS+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000"

	# Which tag in pickled job data to use to generate spectrum plots?
	SPECTAG="output_prog_h_pert_t00000864000.00000000_spectrum"

	# Output plot file
	OUTPUTFILE="plot_prog_h_pert_spectrum_simtime_10days_C1.125e+02.pdf"

	# Plot!
	../pp_plot_planedata_spectrum.py $SPECTAG $OUTPUTFILE $DIRS

fi
