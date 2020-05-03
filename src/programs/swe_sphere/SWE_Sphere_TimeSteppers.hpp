/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "../swe_sphere/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "SWE_Sphere_TS_interface.hpp"

/*
 * Linear only
 */
#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_cn.hpp"
#include "SWE_Sphere_TS_lg_cn.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_erk_ver01.hpp"

/*
 * Full nonlinear
 */
#include "SWE_Sphere_TS_ln_erk.hpp"
#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"
#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_l_exp_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_exp_n_etdrk.hpp"
#include "SWE_Sphere_TS_l_irk_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"

/*
 * Almost full nonlinear
 */

#include "SWE_Sphere_TS_l_irk_na_erk_vd_ver01.hpp"


/*
 * Exponential
 */
#include "SWE_Sphere_TS_l_irk_na_sl_nr_settls_vd_only.hpp"
#include "SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only.hpp"
#include "SWE_Sphere_TS_lg_exp_lc_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_lg_exp_lc_n_etdrk.hpp"
#include "SWE_Sphere_TS_ln_settls_vd.hpp"
#include "SWE_Sphere_TS_ln_sl_exp_settls_vd.hpp"


/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
	/*
	 * Linear
	 */
	SWE_Sphere_TS_l_erk *l_erk = nullptr;
	SWE_Sphere_TS_lg_erk *lg_erk = nullptr;
	SWE_Sphere_TS_lg_irk *lg_irk = nullptr;
	SWE_Sphere_TS_lg_cn *lg_cn = nullptr;

	SWE_Sphere_TS_l_irk *l_irk = nullptr;
	SWE_Sphere_TS_l_exp *l_rexi = nullptr;
	SWE_Sphere_TS_l_exp *lg_rexi = nullptr;

	SWE_Sphere_TS_l_cn *l_cn = nullptr;

	SWE_Sphere_TS_lg_erk_lc_erk *lg_erk_lc_erk = nullptr;
	SWE_Sphere_TS_lg_irk_lc_erk *lg_irk_lc_erk = nullptr;

	/*
	 * Full nonlinear
	 */
	SWE_Sphere_TS_ln_erk *ln_erk = nullptr;

	SWE_Sphere_TS_ln_erk_split_uv *ln_erk_split_uv = nullptr;
	SWE_Sphere_TS_ln_erk_split_vd *ln_erk_split_vd = nullptr;

	SWE_Sphere_TS_lg_erk_lc_n_erk *lg_erk_lc_n_erk = nullptr;
	SWE_Sphere_TS_lg_irk_lc_n_erk *lg_irk_lc_n_erk = nullptr;

	SWE_Sphere_TS_l_erk_n_erk *l_erk_n_erk = nullptr;
	SWE_Sphere_TS_l_irk_n_erk *l_irk_n_erk = nullptr;

	/*
	 * A little bit nonlinear
	 */

	SWE_Sphere_TS_l_irk_na_erk_vd *l_irk_na_erk_vd = nullptr;

	/*
	 * EXP
	 */
	SWE_Sphere_TS_l_exp_n_erk *l_rexi_n_erk = nullptr;
	SWE_Sphere_TS_l_exp_n_etdrk *l_rexi_n_etdrk = nullptr;
	SWE_Sphere_TS_lg_exp_lc_n_erk *lg_rexi_lc_n_erk = nullptr;
	SWE_Sphere_TS_lg_exp_lc_n_etdrk *lg_rexi_lc_n_etdrk = nullptr;


	/*
	 * Semi Lagrangian
	 */
	SWE_Sphere_TS_l_irk_na_sl_nr_settls_vd_only *l_irk_na_sl_nr_settls_vd_only = nullptr;
	SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only *l_irk_na_sl_nr_settls_uv_only = nullptr;
	SWE_Sphere_TS_ln_settls_vd *ln_sl_settls = nullptr;


	/*
	 * EXP SL
	 */
	SWE_Sphere_TS_ln_sl_exp_settls_vd *ln_sl_exp_settls = nullptr;

	SWE_Sphere_TS_interface *master = nullptr;


	SWE_Sphere_TimeSteppers();

	void reset();

	void setup(
			const std::string &i_timestepping_method,
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	);


	~SWE_Sphere_TimeSteppers();
};




#endif
