


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

        'parareal_errors.prog_phi_pert.t1.0.norm_l1' : '$L_1$ normalized error geopotential pert.',
        'parareal_errors.prog_vrt.t1.0.norm_l1' : '$L_1$ normalized error vorticity',
        'parareal_errors.prog_div.t1.0.norm_l1' : '$L_1$ normalized error divergence',

        'parareal_errors.prog_phi_pert.t1.0.norm_l2' : '$L_2$ normalized error geopotential pert.',
        'parareal_errors.prog_vrt.t1.0.norm_l2' : '$L_2$ normalized error vorticity',
        'parareal_errors.prog_div.t1.0.norm_l2' : '$L_2$ normalized error divergence',

        'parareal_errors.prog_phi_pert.t1.0.norm_linf' : '$L_{\infty}$ normalized error geopotential pert.',
        'parareal_errors.prog_vrt.t1.0.norm_linf' : '$L_{\infty}$ normalized error vorticity',
        'parareal_errors.prog_div.t1.0.norm_linf' : '$L_{\infty}$ normalized error divergence',

        'parareal_errors.spec.prog_phi_pert.t1.0.norm_linf_rnorm_16': '$L_{\infty, R_{norm} = 16}$ normalized error geopotential pert. (spectral)',
        'parareal_errors.spec.prog_phi_pert.t1.0.norm_linf_rnorm_32': '$L_{\infty, R_{norm} = 32}$ normalized error geopotential pert. (spectral)',
        'parareal_errors.spec.vrt.t1.0.norm_linf_rnorm_16': '$L_{\infty, R_{norm} = 16}$ normalized error vorticity (spectral)',
        'parareal_errors.spec.vrt.t1.0.norm_linf_rnorm_32': '$L_{\infty, R_{norm} = 32}$ normalized error vorticity (spectral)',
        'parareal_errors.spec.div.t1.0.norm_linf_rnorm_16': '$L_{\infty, R_{norm} = 16}$ normalized error divergence (spectral)',
        'parareal_errors.spec.div.t1.0.norm_linf_rnorm_32': '$L_{\infty, R_{norm} = 32}$ normalized error divergence (spectral)'

}


abbrev_param_names = {

        'runtime.parareal_coarse_timestepping_method': 'pCtsm',
        'runtime.parareal_coarse_timestepping_order': 'pCtso',
        'runtime.parareal_coarse_timestep_size': 'pCtss',
        'runtime.parareal_coarse_slices': 'pCslc',
        'runtime.iteration': 'pIter'
}


def get_pretty_name(name):
	if name in latex_pretty_names:
		return latex_pretty_names[name]
	else:
		return name.replace('_', '\\_')

def get_abbrev_param_name(name):
	if name in abbrev_param_names:
		return abbrev_param_names[name]
	else:
		return name

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
