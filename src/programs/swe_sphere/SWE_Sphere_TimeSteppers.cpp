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
	delete l_erk;
	l_erk = nullptr;

	delete l_erk_n_erk;
	l_erk_n_erk = nullptr;

	delete l_irk_n_erk;
	l_irk_n_erk = nullptr;

	delete lg_erk_lc_erk;
	lg_erk_lc_erk = nullptr;

	delete lg_irk_lc_erk;
	lg_irk_lc_erk = nullptr;

	delete lg_erk_lc_n_erk;
	lg_erk_lc_n_erk = nullptr;

	delete lg_irk_lc_n_erk;
	lg_irk_lc_n_erk = nullptr;

	delete l_irk;
	l_irk = nullptr;

	delete l_rexi;
	l_rexi = nullptr;

	delete lg_rexi;
	lg_rexi = nullptr;

	delete l_cn;
	l_cn = nullptr;

	delete ln_erk;
	ln_erk = nullptr;

	delete ln_erk_split_uv;
	ln_erk_split_uv = nullptr;

	delete ln_erk_split_vd;
	ln_erk_split_vd = nullptr;

	delete l_rexi_n_erk;
	l_rexi_n_erk = nullptr;

	delete l_rexi_n_etdrk;
	l_rexi_n_etdrk = nullptr;

	delete lg_rexi_lc_n_erk;
	lg_rexi_lc_n_erk = nullptr;

	delete lg_rexi_lc_n_etdrk;
	lg_rexi_lc_n_etdrk = nullptr;

	delete lg_erk;
	lg_erk = nullptr;

	delete lg_irk;
	lg_irk = nullptr;

	delete lg_cn;
	lg_cn = nullptr;

	delete ln_sl_settls;
	ln_sl_settls = nullptr;

	delete l_irk_na_sl_nr_settls_vd_only;
	l_irk_na_sl_nr_settls_vd_only = nullptr;

	delete l_irk_na_sl_nr_settls_uv_only;
	l_irk_na_sl_nr_settls_uv_only = nullptr;

	delete ln_sl_exp_settls;
	ln_sl_exp_settls = nullptr;
}



void SWE_Sphere_TimeSteppers::setup(
		const std::string &i_timestepping_method,
		SphereOperators_SphereData &i_op,
		SimulationVariables &i_simVars
)
{
	if (SWE_Sphere_TS_l_erk::implements_timestepping_method(i_timestepping_method))
	{
		l_erk = new SWE_Sphere_TS_l_erk(i_simVars, i_op);
		l_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_erk;
		return;
	}

	if (SWE_Sphere_TS_l_erk_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		l_erk_n_erk = new SWE_Sphere_TS_l_erk_n_erk(i_simVars, i_op);
		l_erk_n_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_erk_n_erk;
		return;
	}

	if (SWE_Sphere_TS_lg_erk_lc_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_erk_lc_erk = new SWE_Sphere_TS_lg_erk_lc_erk(i_simVars, i_op);
		lg_erk_lc_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_erk_lc_erk;
		return;
	}

	if (SWE_Sphere_TS_lg_irk_lc_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_irk_lc_erk = new SWE_Sphere_TS_lg_irk_lc_erk(i_simVars, i_op);
		lg_irk_lc_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_erk;
		return;
	}

	if (SWE_Sphere_TS_l_irk_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		l_irk_n_erk = new SWE_Sphere_TS_l_irk_n_erk(i_simVars, i_op);
		l_irk_n_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_irk_n_erk;
		return;
	}

	if (SWE_Sphere_TS_l_exp_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi_n_erk = new SWE_Sphere_TS_l_exp_n_erk(i_simVars, i_op);
		l_rexi_n_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_erk;
		return;
	}


	if (SWE_Sphere_TS_lg_irk_lc_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_irk_lc_n_erk = new SWE_Sphere_TS_lg_irk_lc_n_erk(i_simVars, i_op);
		lg_irk_lc_n_erk->setup_auto();

		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_n_erk;
		return;
	}


	if (SWE_Sphere_TS_lg_exp_lc_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_rexi_lc_n_erk = new SWE_Sphere_TS_lg_exp_lc_n_erk(i_simVars, i_op);
		lg_rexi_lc_n_erk->setup_auto();

		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_erk;
		return;
	}



	if (SWE_Sphere_TS_lg_erk_lc_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_erk_lc_n_erk = new SWE_Sphere_TS_lg_erk_lc_n_erk(i_simVars, i_op);
		lg_erk_lc_n_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_erk_lc_n_erk;
		return;
	}

	if (SWE_Sphere_TS_lg_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_erk = new SWE_Sphere_TS_lg_erk(i_simVars, i_op);
		lg_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_erk;
		return;
	}


	if (SWE_Sphere_TS_ln_erk::implements_timestepping_method(i_timestepping_method))
	{
		ln_erk = new SWE_Sphere_TS_ln_erk(i_simVars, i_op);
		ln_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*ln_erk;
		return;
	}



	if (SWE_Sphere_TS_ln_erk_split_uv::implements_timestepping_method(i_timestepping_method))
	{
		ln_erk_split_uv = new SWE_Sphere_TS_ln_erk_split_uv(i_simVars, i_op);
		ln_erk_split_uv->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*ln_erk_split_uv;
		return;
	}



	if (SWE_Sphere_TS_ln_erk_split_vd::implements_timestepping_method(i_timestepping_method))
	{
		ln_erk_split_vd = new SWE_Sphere_TS_ln_erk_split_vd(i_simVars, i_op);
		ln_erk_split_vd->setup_auto();

		master = &(SWE_Sphere_TS_interface&)*ln_erk_split_vd;
		return;
	}

	if (SWE_Sphere_TS_l_exp_n_etdrk::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi_n_etdrk = new SWE_Sphere_TS_l_exp_n_etdrk(i_simVars, i_op);
		l_rexi_n_etdrk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_etdrk;
		return;
	}


	if (SWE_Sphere_TS_lg_exp_lc_n_etdrk::implements_timestepping_method(i_timestepping_method))
	{
		lg_rexi_lc_n_etdrk = new SWE_Sphere_TS_lg_exp_lc_n_etdrk(i_simVars, i_op);
		lg_rexi_lc_n_etdrk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_etdrk;
		return;
	}


	if (SWE_Sphere_TS_l_irk::implements_timestepping_method(i_timestepping_method))
	{
		l_irk = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
		l_irk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_irk;
		return;
	}


	if (SWE_Sphere_TS_lg_irk::implements_timestepping_method(i_timestepping_method))
	{
		lg_irk = new SWE_Sphere_TS_lg_irk(i_simVars, i_op);
		lg_irk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_irk;
		return;
	}

	if (SWE_Sphere_TS_l_cn::implements_timestepping_method(i_timestepping_method))
	{
		l_cn = new SWE_Sphere_TS_l_cn(i_simVars, i_op);
		l_cn->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_cn;
		return;
	}


	if (SWE_Sphere_TS_lg_cn::implements_timestepping_method(i_timestepping_method))
	{
		lg_cn = new SWE_Sphere_TS_lg_cn(i_simVars, i_op);
		lg_cn->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_cn;
		return;
	}


	if (SWE_Sphere_TS_l_exp::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi = new SWE_Sphere_TS_l_exp(i_simVars, i_op);
		l_rexi->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_rexi;
		return;
	}


	/*
	 * EXP SETTLS VERSION
	 */
	if (SWE_Sphere_TS_ln_sl_exp_settls_vd::implements_timestepping_method(i_timestepping_method))
	{
		ln_sl_exp_settls = new SWE_Sphere_TS_ln_sl_exp_settls_vd(i_simVars, i_op, true);
		master = &(SWE_Sphere_TS_interface&)*ln_sl_exp_settls;
		return;
	}


	/*
	 * ONLY SETTLS VERSION without any special variants
	 */
	if (SWE_Sphere_TS_l_irk_na_sl_nr_settls_vd_only::implements_timestepping_method(i_timestepping_method))
	{
		l_irk_na_sl_nr_settls_vd_only = new SWE_Sphere_TS_l_irk_na_sl_nr_settls_vd_only(i_simVars, i_op, true);
		master = &(SWE_Sphere_TS_interface&)*l_irk_na_sl_nr_settls_vd_only;
		return;
	}

	if (SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only::implements_timestepping_method(i_timestepping_method))
	{
		l_irk_na_sl_nr_settls_uv_only = new SWE_Sphere_TS_l_irk_na_sl_nr_settls_uv_only(i_simVars, i_op, true);
		master = &(SWE_Sphere_TS_interface&)*l_irk_na_sl_nr_settls_uv_only;
		return;
	}



	/*
	 * IRK SETTLS VERSION
	 */
	if (SWE_Sphere_TS_ln_settls_vd::implements_timestepping_method(i_timestepping_method))
	{
		ln_sl_settls = new SWE_Sphere_TS_ln_settls_vd(i_simVars, i_op, true);
		master = &(SWE_Sphere_TS_interface&)*ln_sl_settls;
		return;
	}

	std::cout << i_timestepping_method << std::endl;
	FatalError("No valid --timestepping-method provided");
}


SWE_Sphere_TimeSteppers::~SWE_Sphere_TimeSteppers()
{
	reset();
}



#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
