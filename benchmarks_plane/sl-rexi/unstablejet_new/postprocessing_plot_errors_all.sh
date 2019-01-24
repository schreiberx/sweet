#! /usr/bin/env bash


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


./postprocessing_plot_errors_single.py plot_errors.pdf $DATA

