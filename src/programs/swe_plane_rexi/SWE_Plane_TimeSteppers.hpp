/*
 * SWE_Plane_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_

#include "SWE_Plane_TS_interface.hpp"

#include "SWE_Plane_TS_l_direct.hpp"
#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_l_cn.hpp"
#include "SWE_Plane_TS_l_erk_n_erk.hpp"
#include "SWE_Plane_TS_l_cn_n_erk.hpp"
#include "SWE_Plane_TS_l_rexi_n_erk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_l_irk_n_erk.hpp"
#include "SWE_Plane_TS_l_rexi.hpp"
#include "SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp"
#include "SWE_Plane_TS_l_cn_na_sl_nd_settls.hpp"
#include "SWE_Plane_TS_ln_erk.hpp"



/**
 * SWE Plane time steppers
 */
class SWE_Plane_TimeSteppers
{
public:
	SWE_Plane_TS_ln_erk *ln_erk = nullptr;
	SWE_Plane_TS_l_erk *l_erk = nullptr;
	SWE_Plane_TS_l_cn *l_cn = nullptr;
	SWE_Plane_TS_l_erk_n_erk *l_erk_n_erk = nullptr;
	SWE_Plane_TS_l_cn_n_erk *l_cn_n_erk = nullptr;
	SWE_Plane_TS_l_rexi_n_erk *l_rexi_n_erk = nullptr;
	SWE_Plane_TS_l_irk *l_irk = nullptr;
	SWE_Plane_TS_l_rexi *l_rexi = nullptr;
	SWE_Plane_TS_l_direct *l_direct = nullptr;
	SWE_Plane_TS_l_rexi_na_sl_nd_settls *l_rexi_na_sl_nd_settls = nullptr;
	SWE_Plane_TS_l_cn_na_sl_nd_settls *l_cn_na_sl_nd_settls = nullptr;

	SWE_Plane_TS_l_irk_n_erk *l_irk_n_erk = nullptr;

	SWE_Plane_TS_interface *master = nullptr;

	bool linear_only = false;

	SWE_Plane_TimeSteppers()
	{
	}

	void reset()
	{
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		if (l_erk != nullptr)
		{
			delete l_erk;
			l_erk = nullptr;
		}

		if (l_cn != nullptr)
		{
			delete l_cn;
			l_cn = nullptr;
		}

		if (l_erk_n_erk != nullptr)
		{
			delete l_erk_n_erk;
			l_erk_n_erk = nullptr;
		}

		if (l_cn_n_erk != nullptr)
		{
			delete l_cn_n_erk;
			l_cn_n_erk = nullptr;
		}

		if (l_rexi_n_erk != nullptr)
		{
			delete l_rexi_n_erk;
			l_rexi_n_erk = nullptr;
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

		if (l_direct != nullptr)
		{
			delete l_direct;
			l_direct = nullptr;
		}

		if (l_rexi_na_sl_nd_settls != nullptr)
		{
			delete l_rexi_na_sl_nd_settls;
			l_rexi_na_sl_nd_settls = nullptr;
		}

		if (l_cn_na_sl_nd_settls != nullptr)
		{
			delete l_cn_na_sl_nd_settls;
			l_cn_na_sl_nd_settls = nullptr;
		}

		if (l_irk_n_erk != nullptr)
		{
			delete l_irk_n_erk;
			l_irk_n_erk = nullptr;
		}
	}



	void setup(
			const std::string &i_timestepping_method,
			int &i_timestepping_order,
			int &i_timestepping_order2,

			PlaneOperators &i_op,
			SimulationVariables &i_simVars
	)
	{
		reset();

		if (i_simVars.sim.CFL >= 0)
			FatalError("Only constant time step size supported with REXI, use negative CFL to set constant time step size");

		/// Always allocate analytical solution
		l_direct = new SWE_Plane_TS_l_direct(i_simVars, i_op);

		if (i_timestepping_method == "ln_erk")
		{
//			if (i_simVars.disc.use_staggering)
//				FatalError("Staggering not supported for ln_erk");

			ln_erk = new SWE_Plane_TS_ln_erk(i_simVars, i_op);
			ln_erk->setup(i_timestepping_order);

			master = &(SWE_Plane_TS_interface&)*ln_erk;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_cn")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_cn");

			l_cn= new SWE_Plane_TS_l_cn(i_simVars, i_op);
			l_cn->setup(i_timestepping_order, i_simVars.disc.crank_nicolson_filter);

			master = &(SWE_Plane_TS_interface&)*l_cn;

			linear_only = true;
		}
		else if (i_timestepping_method == "l_erk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_erk");

			l_erk = new SWE_Plane_TS_l_erk(i_simVars, i_op);
			l_erk->setup(i_timestepping_order);

			master = &(SWE_Plane_TS_interface&)*l_erk;

			linear_only = true;
		}
		else if (i_timestepping_method == "l_erk_n_erk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_erk_n_erk");

			l_erk_n_erk = new SWE_Plane_TS_l_erk_n_erk(i_simVars, i_op);
			l_erk_n_erk->setup(i_timestepping_order, i_timestepping_order2);

			master = &(SWE_Plane_TS_interface&)*l_erk_n_erk;
		}
		else if (i_timestepping_method == "l_cn_n_erk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_cn_n_erk");

			l_cn_n_erk = new SWE_Plane_TS_l_cn_n_erk(i_simVars, i_op);
			l_cn_n_erk->setup(i_timestepping_order, i_timestepping_order2, i_simVars.disc.crank_nicolson_filter);

			master = &(SWE_Plane_TS_interface&)*l_cn_n_erk;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_rexi_n_erk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_rexi_n_erk");

			l_rexi_n_erk = new SWE_Plane_TS_l_rexi_n_erk(i_simVars, i_op);
			l_rexi_n_erk->setup(
					i_simVars.rexi,
					i_timestepping_order2
				);

			master = &(SWE_Plane_TS_interface&)*l_rexi_n_erk;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_irk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_irk");

			l_irk = new SWE_Plane_TS_l_irk(i_simVars, i_op);
			l_irk->setup(i_timestepping_order);

			master = &(SWE_Plane_TS_interface&)*l_irk;

			linear_only = true;
		}
		else if (i_timestepping_method == "l_rexi")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_rexi");

			l_rexi = new SWE_Plane_TS_l_rexi(i_simVars, i_op);
			l_rexi->setup(
					i_simVars.rexi
				);
/*

			if (i_simVars.misc.verbosity > 2)
			{
				std::cout << "ALPHA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi.alpha.size(); n++)
					std::cout << l_rexi->rexi.alpha[n] << std::endl;

				std::cout << "BETA:" << std::endl;
				for (std::size_t n = 0; n < l_rexi->rexi.beta_re.size(); n++)
					std::cout << l_rexi->rexi.beta_re[n] << std::endl;
			}
*/

			master = &(SWE_Plane_TS_interface&)*l_rexi;

			linear_only = true;
		}
		else if (i_timestepping_method == "l_rexi_na_sl_nd_settls")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_rexi_na_sl_nd_settls");

			l_rexi_na_sl_nd_settls = new SWE_Plane_TS_l_rexi_na_sl_nd_settls(i_simVars, i_op);

			l_rexi_na_sl_nd_settls->setup(
					i_simVars.rexi,
					i_simVars.pde.use_nonlinear_equations
				);

			master = &(SWE_Plane_TS_interface&)*l_rexi_na_sl_nd_settls;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_cn_na_sl_nd_settls")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_cn_na_sl_nd_settls");

			l_cn_na_sl_nd_settls = new SWE_Plane_TS_l_cn_na_sl_nd_settls(i_simVars, i_op);

			l_cn_na_sl_nd_settls->setup(
					i_simVars.pde.use_nonlinear_equations
				);

			master = &(SWE_Plane_TS_interface&)*l_cn_na_sl_nd_settls;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_irk_n_erk")
		{
			if (i_simVars.disc.use_staggering)
				FatalError("Staggering not supported for l_irk_n_erk");

			l_irk_n_erk = new SWE_Plane_TS_l_irk_n_erk(i_simVars, i_op);

			l_irk_n_erk->setup(
					i_timestepping_order,
					i_timestepping_order2
				);

			master = &(SWE_Plane_TS_interface&)*l_irk_n_erk;

			linear_only = false;
		}
		else if (i_timestepping_method == "l_direct")
		{
			master = &(SWE_Plane_TS_interface&)*l_direct;

			linear_only = true;
		}
		else
		{
			std::cout << "Unknown method: " << i_timestepping_method << std::endl;
			FatalError("No valid --timestepping-method provided");
		}
	}



	~SWE_Plane_TimeSteppers()
	{
		reset();
	}
};




#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
