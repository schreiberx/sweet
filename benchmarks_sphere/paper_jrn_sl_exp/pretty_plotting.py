


latex_pretty_names = {
	'ln_erk' : "$L(U)+N(U)$: RK",

	'l_irk_n_erk_ver0'	: "$L(U)$: CN, $N(U)$: RK, SSv0",
	'l_irk_n_erk_ver1'	: "$L(U)$: CN, $N(U)$: RK, SSv1",

	'lg_irk_lc_n_erk_ver0'	: "$L_g(U)$: CN, $L_c(U)+N(U)$: RK, SSv0",
	'lg_irk_lc_n_erk_ver1'	: "$L_g(U)$: CN, $L_c(U)+N(U)$: RK, SSv1",

	'l_rexi_n_erk_ver0'	: "$L(U)$: REXI, $N(U)$: RK, SSv0",
	'l_rexi_n_erk_ver1'	: "$L(U)$: REXI, $N(U)$: RK, SSv1",
	'lg_rexi_lc_n_erk_ver0'	: "$L_g(U)$: REXI, $L_c(U)+N(U)$: RK, SSv0",
	'lg_rexi_lc_n_erk_ver1'	: "$L_g(U)$: REXI, $L_c(U)+N(U)$: RK, SSv1",
	'l_rexi_n_etdrk'	: "$L(U)$: REXI, $N(U)$: ETD2RK",
	'lg_rexi_lc_n_etdrk'	: "$L_g(U)$: REXI, $L_c(U)+N(U)$: ETD2RK",

	'l_irk_na_sl_nd_settls'	: "$L(U)$: CN, $N_a$: SL, $N_d$: SETTLS",
	'lg_irk_na_sl_lc_nd_settls'	: "$L_g(U)$: CN, $N_a$: SL, $L_C+N_d$: SETTLS",

	'l_exp_na_sl_nd_settls'	: "$L(U)$: EXP, $N_a$: SL, $N_d$: SETTLS",
	'lg_exp_na_sl_lc_nd_settls'	: "$L_g(U)$: EXP, $N_a$: SL, $L_C+N_d$: SETTLS",

	'sphere_data_diff_prog_h.res_norm_l1'	: '$L_1$ norm surface height (res.\,normalized)',
	'sphere_data_diff_prog_h.res_norm_l2'	: '$L_2$ norm surface height (res.\,normalized)',
	'sphere_data_diff_prog_h.res_norm_linf'	: '$L_{\infty}$ norm surface height',

	'sphere_data_diff_prog_phi_pert.res_norm_l1'	: '$L_1$ norm geopotential (res.\,normalized)',
	'sphere_data_diff_prog_phi_pert.res_norm_l2'	: '$L_2$ norm geopotential (res.\,normalized)',
	'sphere_data_diff_prog_phi_pert.res_norm_linf'	: '$L_{\infty}$ norm surface height',

	'sphere_data_diff_prog_vrt.res_norm_l1'	: '$L_1$ norm vorticity (res.\,normalized)',
	'sphere_data_diff_prog_vrt.res_norm_l2'	: '$L_2$ norm vorticity (res.\,normalized)',
	'sphere_data_diff_prog_vrt.res_norm_linf'	: '$L_{\infty}$ norm vorticity',

	'sphere_data_diff_prog_div.res_norm_l1'	: '$L_1$ norm divergence (res.\,normalized)',
	'sphere_data_diff_prog_div.res_norm_l2'	: '$L_2$ norm divergence (res.\,normalized)',
	'sphere_data_diff_prog_div.res_norm_linf'	: '$L_{\infty}$ norm divergence',
}


def get_pretty_name(name):
	if name in latex_pretty_names:
		return latex_pretty_names[name]
	else:
		return name.replace('_', '\\_')


# order of time stepping methods
tms_order = [
	'ln_erk',

	'l_irk_n_erk_ver0',
	'l_irk_n_erk_ver1',

	'lg_irk_lc_n_erk_ver0',
	'lg_irk_lc_n_erk_ver1',

	'l_rexi_n_erk_ver0',
	'l_rexi_n_erk_ver1',
	'lg_rexi_lc_n_erk_ver0',
	'lg_rexi_lc_n_erk_ver1',
	'l_rexi_n_etdrk',
	'lg_rexi_lc_n_etdrk',

	'l_irk_na_sl_nd_settls',
	'lg_irk_na_sl_lc_nd_settls',

	'l_exp_na_sl_nd_settls',
	'lg_exp_na_sl_lc_nd_settls',
]


def get_pretty_name_order(name):
	try:
		pos = tms_order.index(name)
	except:
		pos = 666

	return str(pos).zfill(3)
