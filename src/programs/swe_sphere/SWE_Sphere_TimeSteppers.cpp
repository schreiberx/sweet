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
	if (l_erk != nullptr)
	{
		delete l_erk;
		l_erk = nullptr;
	}
	if (l_erk_n_erk != nullptr)
	{
		delete l_erk_n_erk;
		l_erk_n_erk = nullptr;
	}
	if (l_erk_na_sl != nullptr)
	{
		delete l_erk_na_sl;
		l_erk_na_sl = nullptr;
	}
	if (l_irk_n_erk != nullptr)
	{
		delete l_irk_n_erk;
		l_irk_n_erk = nullptr;
	}


	if (lg_erk_lc_erk != nullptr)
	{
		delete lg_erk_lc_erk;
		lg_erk_lc_erk = nullptr;
	}
	if (lg_irk_lc_erk != nullptr)
	{
		delete lg_irk_lc_erk;
		lg_irk_lc_erk = nullptr;
	}


	if (lg_erk_lc_n_erk != nullptr)
	{
		delete lg_erk_lc_n_erk;
		lg_erk_lc_n_erk = nullptr;
	}
	if (lg_irk_lc_n_erk != nullptr)
	{
		delete lg_irk_lc_n_erk;
		lg_irk_lc_n_erk = nullptr;
	}


	if (l_irk != nullptr)
	{
		delete l_irk;
		l_irk = nullptr;
	}
	if (l_rexi != nullptr)
	{
		delete l_rexi;
		l_rexi = nullptr;
	}
	if (lg_rexi != nullptr)
	{
		delete lg_rexi;
		lg_rexi = nullptr;
	}

	if (l_cn != nullptr)
	{
		delete l_cn;
		l_cn = nullptr;
	}
	if (ln_erk != nullptr)
	{
		delete ln_erk;
		ln_erk = nullptr;
	}
	if (ln_erk_split_uv != nullptr)
	{
		delete ln_erk_split_uv;
		ln_erk_split_uv = nullptr;
	}
	if (ln_erk_split_vd != nullptr)
	{
		delete ln_erk_split_vd;
		ln_erk_split_vd = nullptr;
	}


	if (na_erk != nullptr)
	{
		delete na_erk;
		na_erk = nullptr;
	}
	if (na_sl_vd != nullptr)
	{
		delete na_sl_vd;
		na_sl_vd = nullptr;
	}
	if (na_sl_uv != nullptr)
	{
		delete na_sl_uv;
		na_sl_uv = nullptr;
	}

	if (l_rexi_n_erk != nullptr)
	{
		delete l_rexi_n_erk;
		l_rexi_n_erk = nullptr;
	}
	if (l_rexi_n_etdrk != nullptr)
	{
		delete l_rexi_n_etdrk;
		l_rexi_n_etdrk = nullptr;
	}
	if (lg_rexi_lc_n_erk != nullptr)
	{
		delete lg_rexi_lc_n_erk;
		lg_rexi_lc_n_erk = nullptr;
	}
	if (lg_rexi_lc_n_etdrk != nullptr)
	{
		delete lg_rexi_lc_n_etdrk;
		lg_rexi_lc_n_etdrk = nullptr;
	}


	if (lg_erk != nullptr)
	{
		delete lg_erk;
		lg_erk = nullptr;
	}
	if (lg_irk != nullptr)
	{
		delete lg_irk;
		lg_irk = nullptr;
	}
	if (lg_cn != nullptr)
	{
		delete lg_cn;
		lg_cn = nullptr;
	}

	if (ln_sl_settls != nullptr)
	{
		delete ln_sl_settls;
		ln_sl_settls = nullptr;
	}

	if (l_irk_na_sl_settls != nullptr)
	{
		delete l_irk_na_sl_settls;
		l_irk_na_sl_settls = nullptr;
	}

	if (ln_sl_exp_settls != nullptr)
	{
		delete ln_sl_exp_settls;
		ln_sl_exp_settls = nullptr;
	}
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

	if (SWE_Sphere_TS_l_erk_na_sl::implements_timestepping_method(i_timestepping_method))
	{
		l_erk_na_sl = new SWE_Sphere_TS_l_erk_na_sl(i_simVars, i_op);
		l_erk_na_sl->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_erk_na_sl;
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

	if (SWE_Sphere_TS_l_rexi_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi_n_erk = new SWE_Sphere_TS_l_rexi_n_erk(i_simVars, i_op);
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


	if (SWE_Sphere_TS_lg_rexi_lc_n_erk::implements_timestepping_method(i_timestepping_method))
	{
		lg_rexi_lc_n_erk = new SWE_Sphere_TS_lg_rexi_lc_n_erk(i_simVars, i_op);
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


	if (SWE_Sphere_TS_na_erk::implements_timestepping_method(i_timestepping_method))
	{
		na_erk = new SWE_Sphere_TS_na_erk(i_simVars, i_op);
		na_erk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*na_erk;
		return;
	}



	if (SWE_Sphere_TS_na_sl_uv::implements_timestepping_method(i_timestepping_method))
	{
		na_sl_uv = new SWE_Sphere_TS_na_sl_uv(i_simVars, i_op);
		na_sl_uv->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*na_sl_uv;
		return;
	}


	if (SWE_Sphere_TS_na_sl_vd::implements_timestepping_method(i_timestepping_method))
	{
		na_sl_vd = new SWE_Sphere_TS_na_sl_vd(i_simVars, i_op);
		na_sl_vd->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*na_sl_vd;
		return;
	}


	if (SWE_Sphere_TS_l_rexi_n_etdrk::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi_n_etdrk = new SWE_Sphere_TS_l_rexi_n_etdrk(i_simVars, i_op);
		l_rexi_n_etdrk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_etdrk;
		return;
	}


	if (SWE_Sphere_TS_lg_rexi_lc_n_etdrk::implements_timestepping_method(i_timestepping_method))
	if (i_timestepping_method == "lg_rexi_lc_n_etdrk")
	{
		lg_rexi_lc_n_etdrk = new SWE_Sphere_TS_lg_rexi_lc_n_etdrk(i_simVars, i_op);
		l_rexi_n_etdrk->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_etdrk;
		return;
	}


	if (SWE_Sphere_TS_l_irk::implements_timestepping_method(i_timestepping_method))
	if (i_timestepping_method == "l_irk")
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


	if (SWE_Sphere_TS_l_rexi::implements_timestepping_method(i_timestepping_method))
	{
		l_rexi = new SWE_Sphere_TS_l_rexi(i_simVars, i_op);
		l_rexi->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_rexi;
		return;
	}


	if (SWE_Sphere_TS_l_irk_na_sl_settls_only::implements_timestepping_method(i_timestepping_method))
	{
		l_irk_na_sl_settls = new SWE_Sphere_TS_l_irk_na_sl_settls_only(i_simVars, i_op);
		l_irk_na_sl_settls->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*l_irk_na_sl_settls;
		return;
	}


	/*
	 * EXP SETTLS VERSION
	 */
	if (SWE_Sphere_TS_ln_sl_exp_settls::implements_timestepping_method(i_timestepping_method))
	{
		ln_sl_exp_settls = new SWE_Sphere_TS_ln_sl_exp_settls(i_simVars, i_op);
		ln_sl_exp_settls->setup_auto();
		master = &(SWE_Sphere_TS_interface&)*ln_sl_exp_settls;
		return;
	}


	/*
	 * IRK SETTLS VERSION
	 */
	if (SWE_Sphere_TS_ln_settls::implements_timestepping_method(i_timestepping_method))
	{
		ln_sl_settls = new SWE_Sphere_TS_ln_settls(i_simVars, i_op);
		ln_sl_settls->setup_auto();
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
