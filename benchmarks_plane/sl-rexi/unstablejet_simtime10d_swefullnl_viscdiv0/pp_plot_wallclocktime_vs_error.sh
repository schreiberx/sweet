#! /bin/bash

DATA=""

# SL-SI-SETTLS
DATA+=" *tsm_l_cn_na_sl_nd_settls_tso2*"

# RK-FDC
DATA+=" *tsm_ln_erk_tso2*"

# SL-ETD2RK
DATA+=" *tsm_l_rexi_na_sl_nd_etdrk_tso2*"

# SL-EXP-SETTLS
DATA+=" *tsm_l_rexi_na_sl_nd_settls_tso2*"

# ETD2RK
DATA+=" *tsm_l_rexi_n_etdrk_tso2*"

# ERK
DATA+=" *tsm_ln_erk_tso2*"

# REXI+ERK
DATA+=" *tsm_l_rexi_n_erk_tso2*"


#../pp_plot_errors_single.py			\

../pp_plot_wallclocktime_vs_error.py	\
	plane_data_diff_prog_h_pert.norm_linf	\
	plot_errors_wallclocktime_prog_h_10days.pdf			\
	$DATA

