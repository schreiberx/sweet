/*
 * SWE_Plane_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_

#include <sweet/core/ErrorBase.hpp>

#include "time/PDESWEPlaneTS_BaseInterface.hpp"
#include "time/SWE_Plane_TS_l_cn.hpp"
#include "time/SWE_Plane_TS_l_cn_n_erk.hpp"
#include "time/SWE_Plane_TS_l_direct.hpp"
#include "time/SWE_Plane_TS_l_erk.hpp"
#include "time/SWE_Plane_TS_l_erk_n_erk.hpp"
#include "time/SWE_Plane_TS_l_irk.hpp"
#include "time/SWE_Plane_TS_l_irk_n_erk.hpp"
#include "time/SWE_Plane_TS_l_rexi.hpp"
#include "time/SWE_Plane_TS_l_rexi_n_erk.hpp"
#include "time/SWE_Plane_TS_l_rexi_n_etdrk.hpp"
#include "time/SWE_Plane_TS_ln_erk.hpp"

#include "time/SWE_Plane_TS_l_cn_na_sl_nr_settls.hpp"
#include "time/SWE_Plane_TS_l_rexi_na_sl_nr_etdrk.hpp"
#include "time/SWE_Plane_TS_l_rexi_na_sl_nr_settls.hpp"

#include "time/ShackPDESWEPlaneTimeDiscretization.hpp"

/**
 * SWE Plane time steppers
 */
class PDESWEPlaneTimeSteppers
{
public:
	sweet::ErrorBase error;

	SWE_Plane_TS_l_direct *l_direct = nullptr;
#if 0
	SWE_Plane_TS_ln_erk *ln_erk = nullptr;
	SWE_Plane_TS_l_rexi_n_etdrk *l_rexi_n_etdrk = nullptr;
	SWE_Plane_TS_l_erk *l_erk = nullptr;
	SWE_Plane_TS_l_cn *l_cn = nullptr;
	SWE_Plane_TS_l_erk_n_erk *l_erk_n_erk = nullptr;
	SWE_Plane_TS_l_cn_n_erk *l_cn_n_erk = nullptr;
	SWE_Plane_TS_l_rexi_n_erk *l_rexi_n_erk = nullptr;
	SWE_Plane_TS_l_irk *l_irk = nullptr;
	SWE_Plane_TS_l_rexi *l_rexi = nullptr;
	SWE_Plane_TS_l_rexi_na_sl_nr_settls *l_rexi_na_sl_nd_settls = nullptr;
	SWE_Plane_TS_l_cn_na_sl_nr_settls *l_cn_na_sl_nd_settls = nullptr;
	SWE_Plane_TS_l_rexi_na_sl_nr_etdrk *l_rexi_na_sl_nd_etdrk = nullptr;
	SWE_Plane_TS_l_irk_n_erk *l_irk_n_erk = nullptr;
#endif

	PDESWEPlaneTS_BaseInterface *timestepper = nullptr;

	ShackPDESWEPlaneTimeDiscretization *shackPDESWEPlaneTimeDiscretization;

	bool linear_only = false;

	PDESWEPlaneTimeSteppers()
	{
	}

	void clear()
	{
		if (timestepper != nullptr)
		{
			delete timestepper;
			timestepper = nullptr;
		}

		if (l_direct != nullptr)
		{
			delete l_direct;
			l_direct = nullptr;
		}
#if 0
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		if (l_rexi_n_etdrk != nullptr)
		{
			delete l_rexi_n_etdrk;
			l_rexi_n_etdrk = nullptr;
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


		if (l_rexi_na_sl_nd_settls != nullptr)
		{
			delete l_rexi_na_sl_nd_settls;
			l_rexi_na_sl_nd_settls = nullptr;
		}
		if (l_rexi_na_sl_nd_etdrk != nullptr)
		{
			delete l_rexi_na_sl_nd_etdrk;
			l_rexi_na_sl_nd_etdrk = nullptr;
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
#endif
	}

	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackPDESWEPlaneTimeDiscretization = io_shackDict.getAutoRegistration<ShackPDESWEPlaneTimeDiscretization>();

		SWE_Plane_TS_l_direct dummy;
		dummy.shackRegistration(&io_shackDict);

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(io_shackDict);
		return true;
	}

	bool setup(
			sweet::PlaneOperators *io_ops,
			sweet::ShackDictionary *io_shackDict
	)
	{
		assert(io_ops != nullptr);
		assert(io_shackDict != nullptr);

		const std::string &timestepping_method = shackPDESWEPlaneTimeDiscretization->timestepping_method;

		linear_only = false;

		if (timestepping_method == "ln_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_ln_erk);

			linear_only = false;
		}
		else if (timestepping_method == "l_rexi_n_etdrk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi_n_etdrk);

			linear_only = false;
		}
		else if (timestepping_method == "l_cn")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_cn);

			linear_only = true;
		}
		else if (timestepping_method == "l_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_erk);

			linear_only = true;
		}
		else if (timestepping_method == "l_erk_na_nd2_erk")
		{
			/*
			 * Special case which treats div(u*h) as u.div(h)
			 */
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_erk_n_erk);
		}
		else if (timestepping_method == "l_erk_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_erk_n_erk);
		}
		else if (timestepping_method == "l_cn_na_nd2_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_cn_n_erk);
		}
		else if (timestepping_method == "l_cn_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_cn_n_erk);
		}
		else if (timestepping_method == "l_rexi_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi_n_erk);
		}
		else if (timestepping_method == "l_irk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_irk);

			linear_only = true;
		}
		else if (timestepping_method == "l_rexi")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi);

			linear_only = true;
		}
		else if (timestepping_method == "l_rexi_na_sl_nd_settls")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi_na_sl_nr_settls);
		}
		else if (timestepping_method == "l_rexi_na_sl_nd_etdrk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi_na_sl_nr_etdrk);
		}
		else if (timestepping_method == "l_cn_na_sl_ld2_settls")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_cn_na_sl_nr_settls);
		}
		else if (timestepping_method == "l_cn_na_sl_nd_settls")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_cn_na_sl_nr_settls);
		}
		else if (timestepping_method == "l_irk_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_irk_n_erk);
		}
		else if (timestepping_method == "l_direct")
		{
			timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_direct);
			linear_only = true;
		}
		else //Help menu with list of schemes
		{
			std::cout << "Unknown method: " << timestepping_method << std::endl;
			std::cout << "Available --timestepping-method :"  << std::endl;
			std::cout << "      l_direct       : Linear:     Analytical solution to linear SW operator"  << std::endl;
			std::cout << "      l_erk          : Linear:     Explicit RK scheme (supports FD-C staggering)"  << std::endl;
			std::cout << "      l_erk_n_erk    : Non-linear: Linear RK, Non-linear RK, Strang-split"  << std::endl;
			std::cout << "      l_cn           : Linear:     Crank-Nicolson (CN) scheme"  << std::endl;
			std::cout << "      l_cn_n_erk     : Non-linear: Linear CN, Non-linear RK, Strang-split" << std::endl;
			std::cout << "      l_rexi         : Linear:     Pure REXI, our little dog." << std::endl;
			std::cout << "      l_rexi_n_erk   : Non-linear: Linear REXI, Non-linear RK, Strang-split" << std::endl;
			std::cout << "      l_irk          : Linear:     Implicit Euler"  << std::endl;
			std::cout << "      l_irk_n_erk    : Non-linear: Linear Implicit Euler, Non-linear RK, Strang-split"  << std::endl;
			std::cout << "      ln_erk         : Non-linear: Linear and nonlinear solved jointly with erk (supports FD-C staggering)"  << std::endl;
			std::cout << "      l_rexi_na_sl_nd_settls   : Non-linear: Linear Rexi, Advection: Semi-Lag, Nonlinear-diverg: SETTLS"  << std::endl;
			std::cout << "      l_cn_na_sl_nd_settls     : Non-linear: Linear CN, Advection: Semi-Lag, Nonlinear-diverg: SETTLS"  << std::endl;
			std::cout << "      l_rexi_n_etdrk           : Non-linear: Linear REXI, Non-linear: ETDRK"  << std::endl;
			std::cout << "      l_rexi_na_sl_nd_etdrk    : Non-linear: Linear REXI, Advection: Semi-Lag, Nonlinear-diverg: ETDRK"  << std::endl;

			SWEETError("No valid --timestepping-method provided");
		}

		timestepper->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timestepper);

		timestepper->setup(io_ops);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timestepper);


		/*
		 * Make analytical solution available
		 */
		if (linear_only)
		{
#if 0
			l_direct = new SWE_Plane_TS_l_direct;

			l_direct->shackRegistration(io_shackDict);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*l_direct);

			l_direct->setup(io_ops);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*l_direct);
#endif
		}

		return true;
	}



	~PDESWEPlaneTimeSteppers()
	{
		clear();
	}
};




#endif
