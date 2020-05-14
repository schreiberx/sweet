/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "SWE_Sphere_TS_interface.hpp"
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
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
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
#include "SWE_Sphere_TS_ln_settls_uv.hpp"
#include "SWE_Sphere_TS_ln_sl_exp_settls_vd.hpp"
#include "SWE_Sphere_TS_ln_sl_exp_settls_uv.hpp"


/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
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
