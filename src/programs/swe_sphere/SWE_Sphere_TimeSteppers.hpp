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

#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"
#include "SWE_Sphere_TS_l_lf.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_cn.hpp"
#include "SWE_Sphere_TS_lg_cn.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"

#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_irk_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_ln_erk.hpp"
#include "SWE_Sphere_TS_l_na_erk.hpp"
#include "SWE_Sphere_TS_l_rexi_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_rexi_n_etdrk.hpp"
#include "SWE_Sphere_TS_lg_rexi_lc_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_lg_rexi_lc_n_etdrk.hpp"


#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_erk_ver01.hpp"

#include "SWE_Sphere_TS_ln_settls.hpp"
#include "SWE_Sphere_TS_ln_sl_exp_settls.hpp"


/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
	SWE_Sphere_TS_l_erk *l_erk = nullptr;

	SWE_Sphere_TS_l_erk_n_erk *l_erk_n_erk = nullptr;
	SWE_Sphere_TS_l_irk_n_erk *l_irk_n_erk = nullptr;

	SWE_Sphere_TS_lg_erk_lc_erk *lg_erk_lc_erk = nullptr;
	SWE_Sphere_TS_lg_irk_lc_erk *lg_irk_lc_erk = nullptr;

	SWE_Sphere_TS_lg_erk_lc_n_erk *lg_erk_lc_n_erk = nullptr;
	SWE_Sphere_TS_lg_irk_lc_n_erk *lg_irk_lc_n_erk = nullptr;

	SWE_Sphere_TS_l_irk *l_irk = nullptr;
	SWE_Sphere_TS_l_lf *l_leapfrog = nullptr;
	SWE_Sphere_TS_l_rexi *l_rexi = nullptr;
	SWE_Sphere_TS_l_rexi *lg_rexi = nullptr;

	SWE_Sphere_TS_l_cn *l_cn = nullptr;
	SWE_Sphere_TS_ln_erk *ln_erk = nullptr;

	SWE_Sphere_TS_l_na_erk *l_na_erk = nullptr;
	SWE_Sphere_TS_l_rexi_n_erk *l_rexi_n_erk = nullptr;
	SWE_Sphere_TS_l_rexi_n_etdrk *l_rexi_n_etdrk = nullptr;
	SWE_Sphere_TS_lg_rexi_lc_n_erk *lg_rexi_lc_n_erk = nullptr;
	SWE_Sphere_TS_lg_rexi_lc_n_etdrk *lg_rexi_lc_n_etdrk = nullptr;

	SWE_Sphere_TS_lg_erk *lg_erk = nullptr;
	SWE_Sphere_TS_lg_irk *lg_irk = nullptr;
	SWE_Sphere_TS_lg_cn *lg_cn = nullptr;

	SWE_Sphere_TS_ln_settls *ln_sl_settls = nullptr;
	SWE_Sphere_TS_ln_sl_exp_settls *ln_sl_exp_settls = nullptr;

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
