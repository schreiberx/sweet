


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

	'sphere_data_diff_prog_h.res_norm_l1'	: '$L_1$ norm surface height (res.\,normalized)',
	'sphere_data_diff_prog_h.res_norm_l2'	: '$L_2$ norm surface height (res.\,normalized)',
	'sphere_data_diff_prog_h.res_norm_linf'	: '$L_{\infty}$ norm surface height',
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
]


def get_pretty_name_order(name):
	pos = tms_order.index(name)
	if pos >= 0:
		return str(pos).zfill(3)
	return 666
