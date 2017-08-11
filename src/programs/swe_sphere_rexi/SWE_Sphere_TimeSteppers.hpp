/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "SWE_Sphere_TS_interface.hpp"

#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_l_lf.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_cn.hpp"
#include "SWE_Sphere_TS_lg_cn.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"
#include "SWE_Sphere_TS_ln_erk.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"



/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
	SWE_Sphere_TS_l_erk *l_erk = nullptr;
	SWE_Sphere_TS_lg_erk *lg_erk = nullptr;
	SWE_Sphere_TS_ln_erk *ln_erk = nullptr;
	SWE_Sphere_TS_l_irk *l_irk = nullptr;
	SWE_Sphere_TS_lg_irk *lg_irk = nullptr;
	SWE_Sphere_TS_l_cn *l_cn = nullptr;
	SWE_Sphere_TS_lg_cn *lg_cn = nullptr;
	SWE_Sphere_TS_l_lf *l_lf = nullptr;
	SWE_Sphere_TS_l_rexi *l_rexi = nullptr;

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
		if (l_irk != nullptr)
		{
			delete l_irk;
			l_irk = nullptr;
		}
		if (lg_irk != nullptr)
		{
			delete lg_irk;
			lg_irk = nullptr;
		}
		if (l_cn != nullptr)
		{
			delete l_cn;
			l_cn = nullptr;
		}
		if (lg_cn != nullptr)
		{
			delete lg_cn;
			lg_cn = nullptr;
		}
		if (l_lf != nullptr)
		{
			delete l_erk;
			l_erk = nullptr;
		}
		if (lg_erk != nullptr)
		{
			delete lg_erk;
			lg_erk = nullptr;
		}
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		if (l_rexi != nullptr)
		{
			delete l_rexi;
			l_rexi = nullptr;
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
		else if (i_timestepping_method == "lg_erk")
		{
			lg_erk = new SWE_Sphere_TS_lg_erk(i_simVars, i_op);
			lg_erk->setup(i_simVars.disc.timestepping_order);

			master = &(SWE_Sphere_TS_interface&)*lg_erk;
		}
		else if (i_timestepping_method == "ln_erk")
		{
			ln_erk = new SWE_Sphere_TS_ln_erk(i_simVars, i_op);
			ln_erk->setup(i_simVars.disc.timestepping_order);

			master = &(SWE_Sphere_TS_interface&)*ln_erk;
		}
		else if (i_timestepping_method == "l_irk")
		{
			l_irk = new SWE_Sphere_TS_l_irk(i_simVars, i_op);
			l_irk->setup(i_simVars.disc.timestepping_order, -i_simVars.sim.CFL, i_simVars.rexi.use_sphere_extended_modes);

			master = &(SWE_Sphere_TS_interface&)*l_irk;
		}
		else if (i_timestepping_method == "lg_irk")
		{
			lg_irk = new SWE_Sphere_TS_lg_irk(i_simVars, i_op);
			lg_irk->setup(i_simVars.disc.timestepping_order, -i_simVars.sim.CFL);

			master = &(SWE_Sphere_TS_interface&)*lg_irk;
		}
		else if (i_timestepping_method == "l_lf")
		{
			l_lf = new SWE_Sphere_TS_l_lf(i_simVars, i_op);
			l_lf->setup(i_simVars.disc.timestepping_order, i_simVars.disc.leapfrog_robert_asselin_filter);

			master = &(SWE_Sphere_TS_interface&)*l_lf;
		}
		else if (i_timestepping_method == "l_cn")
		{
			l_cn = new SWE_Sphere_TS_l_cn(i_simVars, i_op);
			l_cn->setup(i_simVars.disc.crank_nicolson_filter, -i_simVars.sim.CFL, i_simVars.rexi.use_sphere_extended_modes);

			master = &(SWE_Sphere_TS_interface&)*l_cn;
		}
		else if (i_timestepping_method == "lg_cn")
		{
			lg_cn = new SWE_Sphere_TS_lg_cn(i_simVars, i_op);
			lg_cn->setup(i_simVars.disc.crank_nicolson_filter, -i_simVars.sim.CFL);

			master = &(SWE_Sphere_TS_interface&)*lg_cn;
		}
		else if (i_timestepping_method == "l_rexi")
		{
			l_rexi = new SWE_Sphere_TS_l_rexi(i_simVars, i_op);
			l_rexi->setup(
					i_simVars.rexi.h,
					i_simVars.rexi.M,
					i_simVars.rexi.L,

					-i_simVars.sim.CFL,
					i_simVars.rexi.use_half_poles,
					i_simVars.rexi.use_sphere_extended_modes,
					i_simVars.rexi.normalization,
					i_simVars.sim.f_sphere,
					i_simVars.rexi.sphere_solver_preallocation
				);

			if (i_simVars.misc.verbosity > 2)
			{
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi.alpha.size(); n++)
					std::cout << l_rexi->rexi.alpha[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi.beta_re.size(); n++)
					std::cout << l_rexi->rexi.beta_re[n] << std::endl;
			}

			master = &(SWE_Sphere_TS_interface&)*l_rexi;
		}
/*
		else if (i_timestepping_method == "l_rexi_ns_sl_nd_erk")
		{
			l_rexi_ns_sl_nd_erk = new SWE_Sphere_TS_l_rexi_ns_sl_nd_erk(i_simVars, i_op);

			l_rexi_ns_sl_nd_erk->setup(
					i_simVars.rexi.rexi_h,
					i_simVars.rexi.rexi_M,
					i_simVars.rexi.rexi_L,
					i_simVars.rexi.rexi_use_half_poles,
					i_simVars.rexi.rexi_normalization,
					i_simVars.pde.use_nonlinear_equations
				);

			master = &(SWE_Sphere_TS_interface&)*l_rexi_ns_sl_nd_erk;
		}
		else if (i_timestepping_method == "lg_rexi_lc_erk_nt_sl_nd_erk")
		{
			lg_rexi_lc_erk_nt_sl_nd_erk = new SWE_Sphere_TS_lg_rexi_lc_erk_nt_sl_nd_erk(i_simVars, i_op);

			lg_rexi_lc_erk_nt_sl_nd_erk->setup(
					i_simVars.rexi.rexi_h,
					i_simVars.rexi.rexi_M,
					i_simVars.rexi.rexi_L,
					i_simVars.rexi.rexi_use_half_poles,
					i_simVars.rexi.rexi_normalization,
					i_simVars.pde.use_nonlinear_equations
				);

			master = &(SWE_Sphere_TS_interface&)*lg_rexi_lc_erk_nt_sl_nd_erk;
		}
		else if (i_timestepping_method == "l_direct")
		{
			master = &(SWE_Sphere_TS_interface&)*l_direct;
		}
*/
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
