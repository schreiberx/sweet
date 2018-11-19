/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "../swe_sphere/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "SWE_Sphere_TS_interface.hpp"

#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_l_erk_pvd.hpp"
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


/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
	SWE_Sphere_TS_l_erk *l_erk = nullptr;

	SWE_Sphere_TS_l_erk_pvd *l_erk_pvd = nullptr;
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

	SWE_Sphere_TS_interface *master = nullptr;



	SWE_Sphere_TimeSteppers()
	{
	}

	void reset()
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
	}



	void setup(
			const std::string &i_timestepping_method,
			SphereOperators &i_op,
			SimulationVariables &i_simVars
	)
	{
		if (i_timestepping_method == "l_erk")
		{
			l_erk = new SWE_Sphere_TS_l_erk(i_simVars, i_op);
			l_erk->setup(i_simVars.disc.timestepping_order);

			master = &(SWE_Sphere_TS_interface&)*l_erk;
		}
		else if (i_timestepping_method == "l_erk_pvd")
		{
			l_erk_pvd = new SWE_Sphere_TS_l_erk_pvd(i_simVars, i_op);
			l_erk_pvd->setup(i_simVars.disc.timestepping_order);

			master = &(SWE_Sphere_TS_interface&)*l_erk_pvd;
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
		else if (i_timestepping_method == "l_irk_n_erk" || i_timestepping_method == "l_irk_n_erk_ver0")
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
			l_leapfrog->setup(i_simVars.disc.timestepping_order, i_simVars.disc.leapfrog_robert_asselin_filter);

			master = &(SWE_Sphere_TS_interface&)*l_leapfrog;
		}
		else if (i_timestepping_method == "l_cn")
		{
			l_cn = new SWE_Sphere_TS_l_cn(i_simVars, i_op);
			l_cn->setup(i_simVars.disc.crank_nicolson_filter, i_simVars.timecontrol.current_timestep_size, i_simVars.rexi.use_sphere_extended_modes);

			master = &(SWE_Sphere_TS_interface&)*l_cn;
		}
		else if (i_timestepping_method == "lg_cn")
		{
			lg_cn = new SWE_Sphere_TS_lg_cn(i_simVars, i_op);
			lg_cn->setup(i_simVars.disc.crank_nicolson_filter, i_simVars.timecontrol.current_timestep_size);

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
				for (std::size_t n = 0; n < l_rexi->rexi_alpha.size(); n++)
					std::cout << l_rexi->rexi_alpha[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi_beta.size(); n++)
					std::cout << l_rexi->rexi_beta[n] << std::endl;
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
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi_alpha.size(); n++)
					std::cout << l_rexi->rexi_alpha[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi_beta.size(); n++)
					std::cout << l_rexi->rexi_beta[n] << std::endl;
			}

			master = &(SWE_Sphere_TS_interface&)*lg_rexi;
		}
		else
		{
			std::cout << i_timestepping_method << std::endl;
			FatalError("No valid --timestepping-method provided");
		}
	}


	~SWE_Sphere_TimeSteppers()
	{
		reset();
	}
};




#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
