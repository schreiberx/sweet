/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "SWE_Sphere_TimeSteppers.hpp"


SWE_Sphere_TimeSteppers::SWE_Sphere_TimeSteppers()
{
}

void SWE_Sphere_TimeSteppers::reset()
{
	delete master;
	master = nullptr;
}



void SWE_Sphere_TimeSteppers::setup(
		const std::string &i_timestepping_method,
		SphereOperators_SphereData &i_op,
		SimulationVariables &i_simVars
)
{
	std::vector <SWE_Sphere_TS_interface*> registered_integrators;

	/*
	 * Register time integrators
	 */
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_erk_n_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_erk_lc_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_irk_lc_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_irk_n_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_irk_na_erk_vd(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_exp_n_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_irk_lc_n_erk(i_simVars, i_op)));

	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_exp_lc_n_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_erk_lc_n_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_erk_split_uv(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_erk_split_vd(i_simVars, i_op)));

	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_exp_n_etdrk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_exp_lc_n_etdrk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_irk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_irk(i_simVars, i_op)));

	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_cn(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_lg_cn(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_exp(i_simVars, i_op)));

	/*
	 * EXP SETTLS VERSION
	 */
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_sl_exp_settls_vd(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_sl_exp_settls_uv(i_simVars, i_op)));

	/*
	 * ONLY SETTLS VERSION without any special variants
	 */
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_irk_na_sl_nr_settls_vd_only(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only(i_simVars, i_op)));

	/*
	 * IRK SETTLS VERSION
	 */
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_settls_vd(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SWE_Sphere_TS_interface*>(new SWE_Sphere_TS_ln_settls_uv(i_simVars, i_op)));



	/*
	 * Find right one
	 */
	master = nullptr;

	for (std::size_t i = 0; i < registered_integrators.size(); i++)
	{
		SWE_Sphere_TS_interface *ts = registered_integrators[i];

		if (ts->implements_timestepping_method(i_timestepping_method))
		{
			if (master != nullptr)
				SWEETError(std::string("Duplicate implementation for method ")+i_timestepping_method);

			ts->setup_auto();
			master = ts;
		}
	}


	// Found integrator, freeing others
	if (master != nullptr)
	{
		for (std::size_t i = 0; i < registered_integrators.size(); i++)
		{
			SWE_Sphere_TS_interface *ts = registered_integrators[i];

			if (ts == master)
				continue;

			delete ts;
		}
	}

	if (master != nullptr)
		return;


	SWEETError(std::string("No valid --timestepping-method '")+i_timestepping_method+std::string("' provided"));
}


SWE_Sphere_TimeSteppers::~SWE_Sphere_TimeSteppers()
{
	reset();
}



#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
