#! /usr/bin/env bash


DATA=""

# Add jobs like this 
DATA+=""

# Reference solution
DATA+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000 "

# SL-SI-SETTLS_dt0112.5
DATA+=" job_bench_RT_tsm_l_cn_na_sl_nd_settls_tso2_tsob2_C1.125e+02_S864000"

# SL-SI-SETTLS_dt0112.5
DATA+=" job_bench_RT_tsm_l_cn_na_sl_nd_settls_tso2_tsob2_C1.125e+02_S864000"

# SL-EXP-SETTLS_dt0112.5
DATA+=" job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C1.125e+02_S864000"

# SL-ETD2RK_dt0112.5
DATA+=" job_bench_RT_tsm_l_rexi_na_sl_nd_etdrk_tso2_tsob2_C1.125e+02_S864000"

# SL-ETD2RK_dt0112.5
DATA+=" job_bench_RT_tsm_l_rexi_na_sl_nd_etdrk_tso2_tsob2_C1.125e+02_S864000"

# ETD2RK_dt0112.5
DATA+=" job_bench_RT_tsm_l_rexi_n_etdrk_tso2_tsob2_C1.125e+02_S864000"



# Plot Energy
../pp_plot_kineticenergy_spectrum_single.py		\
	plane_data_kinetic_energy_t00000864000.00000000		\
	plot_kineticenergy_spectrum_simtime_10days.pdf 		\
	$DATA
	
