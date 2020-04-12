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
	if (l_leapfrog != nullptr)
	{
		delete l_erk;
		l_erk = nullptr;
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

	if (l_na_erk != nullptr)
	{
		delete l_na_erk;
		l_na_erk = nullptr;
	}

	if (na_erk != nullptr)
	{
		delete na_erk;
		na_erk = nullptr;
	}
	if (na_sl_pvd != nullptr)
	{
		delete na_sl_pvd;
		na_sl_pvd = nullptr;
	}
	if (na_sl != nullptr)
	{
		delete na_sl;
		na_sl = nullptr;
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
	if (i_timestepping_method == "l_erk")
	{
		l_erk = new SWE_Sphere_TS_l_erk(i_simVars, i_op);
		l_erk->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*l_erk;
	}
	else if (i_timestepping_method == "l_erk_n_erk")
	{
		l_erk_n_erk = new SWE_Sphere_TS_l_erk_n_erk(i_simVars, i_op);
		l_erk_n_erk->setup(i_simVars.disc.timestepping_order, i_simVars.disc.timestepping_order2);

		master = &(SWE_Sphere_TS_interface&)*l_erk_n_erk;
	}
	else if (i_timestepping_method == "lg_erk_lc_erk")
	{
		lg_erk_lc_erk = new SWE_Sphere_TS_lg_erk_lc_erk(i_simVars, i_op);
		lg_erk_lc_erk->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*lg_erk_lc_erk;
	}
	else if (i_timestepping_method == "lg_irk_lc_erk" || i_timestepping_method == "lg_irk_lc_erk_ver0")
	{
		lg_irk_lc_erk = new SWE_Sphere_TS_lg_irk_lc_erk(i_simVars, i_op);
		lg_irk_lc_erk->setup(i_simVars.disc.timestepping_order, 0);

		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_erk;
	}
	else if (i_timestepping_method == "lg_irk_lc_erk_ver1")
	{
		lg_irk_lc_erk = new SWE_Sphere_TS_lg_irk_lc_erk(i_simVars, i_op);
		lg_irk_lc_erk->setup(i_simVars.disc.timestepping_order, 1);

		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_erk;
	}
	else if (
		i_timestepping_method == "l_irk_n_erk" || i_timestepping_method == "l_irk_n_erk_ver0" ||
		i_timestepping_method == "l_cn_n_erk" || i_timestepping_method == "l_cn_n_erk_ver0"
	)
	{
		l_irk_n_erk = new SWE_Sphere_TS_l_irk_n_erk(i_simVars, i_op);
		l_irk_n_erk->setup(i_simVars.disc.timestepping_order, i_simVars.disc.timestepping_order2, 0);

		master = &(SWE_Sphere_TS_interface&)*l_irk_n_erk;
	}
	else if (i_timestepping_method == "l_irk_n_erk_ver1")
	{
		l_irk_n_erk = new SWE_Sphere_TS_l_irk_n_erk(i_simVars, i_op);
		l_irk_n_erk->setup(i_simVars.disc.timestepping_order, i_simVars.disc.timestepping_order2, 1);

		master = &(SWE_Sphere_TS_interface&)*l_irk_n_erk;
	}
	else if (i_timestepping_method == "l_rexi_n_erk" || i_timestepping_method == "l_rexi_n_erk_ver0")
	{
		l_rexi_n_erk = new SWE_Sphere_TS_l_rexi_n_erk(i_simVars, i_op);

		l_rexi_n_erk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size,
				i_simVars.sim.sphere_use_fsphere,
				0	// VERSION 0
			);

		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_erk;
	}
	else if (i_timestepping_method == "l_rexi_n_erk_ver1")
	{
		l_rexi_n_erk = new SWE_Sphere_TS_l_rexi_n_erk(i_simVars, i_op);

		l_rexi_n_erk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size,
				i_simVars.sim.sphere_use_fsphere,
				1	// VERSION 1
			);

		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_erk;
	}
	else if (i_timestepping_method == "lg_irk_lc_n_erk" || i_timestepping_method == "lg_irk_lc_n_erk_ver0")
	{
		lg_irk_lc_n_erk = new SWE_Sphere_TS_lg_irk_lc_n_erk(i_simVars, i_op);
		lg_irk_lc_n_erk->setup(
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				0
			);

		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_irk_lc_n_erk_ver1")
	{
		lg_irk_lc_n_erk = new SWE_Sphere_TS_lg_irk_lc_n_erk(i_simVars, i_op);
		lg_irk_lc_n_erk->setup(
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				1
			);

		master = &(SWE_Sphere_TS_interface&)*lg_irk_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_rexi_lc_n_erk" || i_timestepping_method == "lg_rexi_lc_n_erk_ver0")
	{
		lg_rexi_lc_n_erk = new SWE_Sphere_TS_lg_rexi_lc_n_erk(i_simVars, i_op);
		lg_rexi_lc_n_erk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size,
				0
			);

		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_rexi_lc_n_erk_ver1")
	{
		lg_rexi_lc_n_erk = new SWE_Sphere_TS_lg_rexi_lc_n_erk(i_simVars, i_op);
		lg_rexi_lc_n_erk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size,
				1
			);

		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_erk_lc_n_erk" || i_timestepping_method == "lg_erk_lc_n_erk_ver0")
	{
		lg_erk_lc_n_erk = new SWE_Sphere_TS_lg_erk_lc_n_erk(i_simVars, i_op);
		lg_erk_lc_n_erk->setup(
				i_simVars.disc.timestepping_order,
				0
			);

		master = &(SWE_Sphere_TS_interface&)*lg_erk_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_erk_lc_n_erk_ver1")
	{
		lg_erk_lc_n_erk = new SWE_Sphere_TS_lg_erk_lc_n_erk(i_simVars, i_op);
		lg_erk_lc_n_erk->setup(
				i_simVars.disc.timestepping_order,
				1
			);

		master = &(SWE_Sphere_TS_interface&)*lg_erk_lc_n_erk;
	}
	else if (i_timestepping_method == "lg_erk")
	{
		lg_erk = new SWE_Sphere_TS_lg_erk(i_simVars, i_op);
		lg_erk->setup(
				i_simVars.disc.timestepping_order
			);

		master = &(SWE_Sphere_TS_interface&)*lg_erk;
	}
	else if (i_timestepping_method == "ln_erk")
	{
		ln_erk = new SWE_Sphere_TS_ln_erk(i_simVars, i_op);
		ln_erk->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*ln_erk;
	}
	else if (i_timestepping_method == "l_na_erk")
	{
		l_na_erk = new SWE_Sphere_TS_l_na_erk(i_simVars, i_op);
		l_na_erk->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*l_na_erk;
	}
	else if (i_timestepping_method == "na_sl")
	{
		na_sl = new SWE_Sphere_TS_na_sl(i_simVars, i_op);
		na_sl->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*na_sl;
	}
	else if (i_timestepping_method == "na_erk")
	{
		na_erk = new SWE_Sphere_TS_na_erk(i_simVars, i_op);
		na_erk->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*na_erk;
	}
	else if (i_timestepping_method == "na_sl" || i_timestepping_method == "na_sl_puv")
	{
		na_sl = new SWE_Sphere_TS_na_sl(i_simVars, i_op);
		na_sl->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*na_sl;
	}
	else if (i_timestepping_method == "na_sl_pvd")
	{
		na_sl_pvd = new SWE_Sphere_TS_na_sl_pvd(i_simVars, i_op);
		na_sl_pvd->setup(i_simVars.disc.timestepping_order);

		master = &(SWE_Sphere_TS_interface&)*na_sl_pvd;
	}
	else if (i_timestepping_method == "l_rexi_n_etdrk")
	{
		l_rexi_n_etdrk = new SWE_Sphere_TS_l_rexi_n_etdrk(i_simVars, i_op);
		l_rexi_n_etdrk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size
			);

		if (i_simVars.sim.sphere_use_fsphere)
			FatalError("TODO: Not yet supported");

		master = &(SWE_Sphere_TS_interface&)*l_rexi_n_etdrk;
	}
	else if (i_timestepping_method == "lg_rexi_lc_n_etdrk")
	{
		lg_rexi_lc_n_etdrk = new SWE_Sphere_TS_lg_rexi_lc_n_etdrk(i_simVars, i_op);
		lg_rexi_lc_n_etdrk->setup(
				i_simVars.rexi,
				i_simVars.disc.timestepping_order,
				i_simVars.disc.timestepping_order2,
				i_simVars.timecontrol.current_timestep_size
			);

		if (i_simVars.sim.sphere_use_fsphere)
			FatalError("TODO: Not yet supported");

		master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_n_etdrk;
	}
	else if (i_timestepping_method == "l_irk")
	{
		l_irk = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
		l_irk->setup(i_simVars.disc.timestepping_order, i_simVars.timecontrol.current_timestep_size, i_simVars.rexi.use_sphere_extended_modes);

		master = &(SWE_Sphere_TS_interface&)*l_irk;
	}
	else if (i_timestepping_method == "lg_irk")
	{
		lg_irk = new SWE_Sphere_TS_lg_irk(i_simVars, i_op);
		lg_irk->setup(i_simVars.disc.timestepping_order, i_simVars.timecontrol.current_timestep_size);

		master = &(SWE_Sphere_TS_interface&)*lg_irk;
	}
	else if (i_timestepping_method == "l_lf")
	{
		l_leapfrog = new SWE_Sphere_TS_l_lf(i_simVars, i_op);
		l_leapfrog->setup(i_simVars.disc.timestepping_order, i_simVars.disc.timestepping_leapfrog_robert_asselin_filter);

		master = &(SWE_Sphere_TS_interface&)*l_leapfrog;
	}
	else if (
		i_timestepping_method == "l_cn"	||
		i_timestepping_method == "l_irk"
	)
	{
		l_cn = new SWE_Sphere_TS_l_cn(i_simVars, i_op);
		l_cn->setup(i_simVars.disc.timestepping_crank_nicolson_filter, i_simVars.timecontrol.current_timestep_size, i_simVars.rexi.use_sphere_extended_modes);

		master = &(SWE_Sphere_TS_interface&)*l_cn;
	}
	else if (
		i_timestepping_method == "lg_cn"	||
		i_timestepping_method == "lg_irk"
	)
	{
		lg_cn = new SWE_Sphere_TS_lg_cn(i_simVars, i_op);
		lg_cn->setup(i_simVars.disc.timestepping_crank_nicolson_filter, i_simVars.timecontrol.current_timestep_size);

		master = &(SWE_Sphere_TS_interface&)*lg_cn;
	}
	else if (i_timestepping_method == "l_rexi")
	{
		l_rexi = new SWE_Sphere_TS_l_rexi(i_simVars, i_op);
		l_rexi->setup(
				i_simVars.rexi,
				"phi0",
				i_simVars.timecontrol.current_timestep_size,
				i_simVars.sim.sphere_use_fsphere,
				false
			);
#if 0
		if (i_simVars.misc.verbosity > 2)
		{
			std::cout << "ALPHA:" << std::endl;
			for (std::size_t n = 0; n < l_rexi->rexi_alphas.size(); n++)
				std::cout << l_rexi->rexi_alphas[n] << std::endl;

			std::cout << "BETA:" << std::endl;
			for (std::size_t n = 0; n < l_rexi->rexi_betas.size(); n++)
				std::cout << l_rexi->rexi_betas[n] << std::endl;
		}
#endif
		master = &(SWE_Sphere_TS_interface&)*l_rexi;
	}
	else if (i_timestepping_method == "lg_rexi")
	{
		lg_rexi = new SWE_Sphere_TS_l_rexi(i_simVars, i_op);
		lg_rexi->setup(
				i_simVars.rexi,
				"phi0",
				i_simVars.timecontrol.current_timestep_size,
				i_simVars.sim.sphere_use_fsphere,
				true
			);

		if (i_simVars.misc.verbosity > 2)
		{
			if (i_simVars.rexi.rexi_method != "direct")
			{
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < lg_rexi->rexi_alphas.size(); n++)
					std::cout << lg_rexi->rexi_alphas[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < lg_rexi->rexi_betas.size(); n++)
					std::cout << lg_rexi->rexi_betas[n] << std::endl;
			}
		}

		master = &(SWE_Sphere_TS_interface&)*lg_rexi;
	}
	else if (
			(i_timestepping_method.find("_settls") != std::string::npos)
			&&
			(i_timestepping_method.find("_exp") != std::string::npos)
	)
	{

		SWE_Sphere_TS_ln_sl_exp_settls::LinearGravityTreatment_enum linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::LinearCoriolisTreatment_enum linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::NLAdvectionTreatment_enum nonlinear_advection_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::NLDivergenceTreatment_enum nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_IGNORE;
		bool original_linear_operator_sl_treatment = true;

		if (i_timestepping_method == "ln_settls")
		{
			linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IMPLICIT;
			linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_LINEAR;
			nonlinear_advection_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN;
			nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_NONLINEAR;
			original_linear_operator_sl_treatment = true;
		}
		else
		{
			// Search for implicit or exp treatment of linear parts
			if (i_timestepping_method.find("_irk_") != std::string::npos)
				linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IMPLICIT;
			else if (i_timestepping_method.find("_exp_") != std::string::npos)
				linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL;

			// Search for Coriolis
			if (i_timestepping_method.find("l_irk") != std::string::npos || i_timestepping_method.find("l_exp") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_LINEAR;
			else if (i_timestepping_method.find("lc_na_sl") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN;
			else if (i_timestepping_method.find("lc_") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_NONLINEAR;

			if (i_timestepping_method.find("_na_sl") != std::string::npos)
				nonlinear_advection_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN;

			// Search for Nonlinear divergence
			if (i_timestepping_method.find("_nd_") != std::string::npos)
				nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_NONLINEAR;

			if (i_timestepping_method.find("_ver2") != std::string::npos)
				original_linear_operator_sl_treatment = false;

#if 1
			std::string id_string = "";

			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_LINEAR)
				id_string += "l";
			else
				id_string += "lg";

#if 0
			if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IMPLICIT)
				id_string += "_irk";
			else if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL)
				id_string += "_exp";
#else
			if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL)
				id_string += "_exp";
#endif
			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN)
				id_string += "_lc";

			if (nonlinear_advection_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN)
				id_string += "_na";

			if (
				linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN ||
				nonlinear_advection_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN
			)
				id_string += "_sl";

			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_NONLINEAR)
				id_string += "_lc";

			if (nonlinear_divergence_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_NONLINEAR)
				id_string += "_nd";

			id_string += "_settls";

			if (!original_linear_operator_sl_treatment)
				id_string += "_ver2";


			if (i_timestepping_method != id_string)
			{
				if (!original_linear_operator_sl_treatment)
				{
					std::cerr << "Detected time stepping method: "+id_string << std::endl;
					std::cerr << "Provided time stepping method: "+i_timestepping_method << std::endl;
					FatalError("Autodetection of parts of time stepping methods failed!");
				}

				std::string id_string2 = id_string+"_ver1";
				if (i_timestepping_method != id_string2)
				{
					std::cerr << "Detected time stepping method: "+id_string << std::endl;
					std::cerr << "Provided time stepping method: "+i_timestepping_method << std::endl;
					std::cerr << "Detected alternative time stepping method: "+id_string << std::endl;
					FatalError("Autodetection of parts of time stepping methods failed!");
				}
			}
		}
#endif

		ln_sl_exp_settls = new SWE_Sphere_TS_ln_sl_exp_settls(i_simVars, i_op);
		ln_sl_exp_settls->setup(
				i_simVars.disc.timestepping_order,
				linear_gravity_treatment,
				linear_coriolis_treatment,					// Coriolis treatment
				nonlinear_advection_treatment,		// Nonlinear advection treatment
				nonlinear_divergence_treatment,		// Nonlinear divergence treatment
				original_linear_operator_sl_treatment			// original SL linear operator treatment
			);

		master = &(SWE_Sphere_TS_interface&)*ln_sl_exp_settls;
	}

	/**
	 * SL methods with Coriolis treated as part of linear terms
	 */
	else if (
			(i_timestepping_method.find("_settls") != std::string::npos)
			&&
			!(i_timestepping_method.find("_exp") != std::string::npos)
	)
	{
		SWE_Sphere_TS_ln_settls::LinearTreatment_enum linear_treatment = SWE_Sphere_TS_ln_settls::LINEAR_IGNORE;
		SWE_Sphere_TS_ln_settls::LinearCoriolisTreatment_enum linear_coriolis_treatment = SWE_Sphere_TS_ln_settls::CORIOLIS_IGNORE;
		SWE_Sphere_TS_ln_settls::NLAdvectionTreatment_enum nonlinear_advection_treatment = SWE_Sphere_TS_ln_settls::NL_ADV_IGNORE;
		SWE_Sphere_TS_ln_settls::NLDivergenceTreatment_enum nonlinear_divergence_treatment = SWE_Sphere_TS_ln_settls::NL_DIV_IGNORE;
		bool original_linear_operator_sl_treatment = true;

		if (i_timestepping_method == "ln_settls")
		{
			linear_treatment = SWE_Sphere_TS_ln_settls::LINEAR_IMPLICIT;
			linear_coriolis_treatment = SWE_Sphere_TS_ln_settls::CORIOLIS_LINEAR;
			nonlinear_divergence_treatment = SWE_Sphere_TS_ln_settls::NL_DIV_NONLINEAR;
			original_linear_operator_sl_treatment = true;
		}
		else
		{
			// Search for implicit or exp treatment of linear parts
			if (i_timestepping_method.find("_irk_") != std::string::npos)
				linear_treatment = SWE_Sphere_TS_ln_settls::LINEAR_IMPLICIT;
			else if (i_timestepping_method.find("_exp_") != std::string::npos)
				linear_treatment = SWE_Sphere_TS_ln_settls::LINEAR_EXPONENTIAL;

			// Search for Coriolis
			if (i_timestepping_method.find("l_irk") != std::string::npos || i_timestepping_method.find("l_exp") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_settls::CORIOLIS_LINEAR;
			else if (i_timestepping_method.find("lc_na_sl") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_settls::CORIOLIS_SEMILAGRANGIAN;
			else if (i_timestepping_method.find("lc_") != std::string::npos)
				linear_coriolis_treatment = SWE_Sphere_TS_ln_settls::CORIOLIS_NONLINEAR;

			if (i_timestepping_method.find("_na_sl") != std::string::npos)
				nonlinear_advection_treatment = SWE_Sphere_TS_ln_settls::NL_ADV_SEMILAGRANGIAN;

			// Search for Nonlinear divergence
			if (i_timestepping_method.find("_nd_") != std::string::npos)
				nonlinear_divergence_treatment = SWE_Sphere_TS_ln_settls::NL_DIV_NONLINEAR;

			if (i_timestepping_method.find("_ver2") != std::string::npos)
				original_linear_operator_sl_treatment = false;

#if 1
			std::string id_string = "";

			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_settls::CORIOLIS_LINEAR)
				id_string += "l";
			else
				id_string += "lg";

			if (linear_treatment == SWE_Sphere_TS_ln_settls::LINEAR_IMPLICIT)
				id_string += "_irk";
			else if (linear_treatment == SWE_Sphere_TS_ln_settls::LINEAR_EXPONENTIAL)
				id_string += "_exp";

			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_settls::CORIOLIS_SEMILAGRANGIAN)
				id_string += "_lc";

			if (nonlinear_advection_treatment == SWE_Sphere_TS_ln_settls::NL_ADV_SEMILAGRANGIAN)
				id_string += "_na";

			if (
				linear_coriolis_treatment == SWE_Sphere_TS_ln_settls::CORIOLIS_SEMILAGRANGIAN ||
				nonlinear_advection_treatment == SWE_Sphere_TS_ln_settls::NL_ADV_SEMILAGRANGIAN
			)
			id_string += "_sl";

			if (linear_coriolis_treatment == SWE_Sphere_TS_ln_settls::CORIOLIS_NONLINEAR)
				id_string += "_lc";

			if (nonlinear_divergence_treatment == SWE_Sphere_TS_ln_settls::NL_DIV_NONLINEAR)
				id_string += "_nd";

			id_string += "_settls";

			if (!original_linear_operator_sl_treatment)
				id_string += "_ver2";


			if (i_timestepping_method != id_string)
			{
				if (!original_linear_operator_sl_treatment)
				{
					std::cerr << "Detected time stepping method: "+id_string << std::endl;
					std::cerr << "Provided time stepping method: "+i_timestepping_method << std::endl;
					FatalError("Autodetection of parts of time stepping methods failed!");
				}

				std::string id_string2 = id_string+"_ver1";
				if (i_timestepping_method != id_string2)
				{
					std::cerr << "Detected time stepping method: "+id_string << std::endl;
					std::cerr << "Provided time stepping method: "+i_timestepping_method << std::endl;
					std::cerr << "Detected alternative time stepping method: "+id_string << std::endl;
					FatalError("Autodetection of parts of time stepping methods failed!");
				}
			}
		}
#endif

		ln_sl_settls = new SWE_Sphere_TS_ln_settls(i_simVars, i_op);
		ln_sl_settls->setup(
				i_simVars.disc.timestepping_order,
				linear_treatment,
				linear_coriolis_treatment,				// Coriolis treatment
				nonlinear_advection_treatment,			// Nonlinear advection treatment
				nonlinear_divergence_treatment,			// Nonlinear divergence treatment
				original_linear_operator_sl_treatment	// original SL linear operator treatment
			);

		master = &(SWE_Sphere_TS_interface&)*ln_sl_settls;
	}
	else
	{
		std::cout << i_timestepping_method << std::endl;
		FatalError("No valid --timestepping-method provided");
	}
}


SWE_Sphere_TimeSteppers::~SWE_Sphere_TimeSteppers()
{
	reset();
}



#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
