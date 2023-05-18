/*
 * PDESWEPlaneMoriZwanzig_TimeSteppers.hpp
 *
 *  Created on: 12 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TIMESTEPPERS_HPP_

#include <sweet/core/ErrorBase.hpp>

#include "time/PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "time/SWE_Plane_Mori_Zwanzig_TS_l_direct.hpp"
#include "time/SWE_Plane_Mori_Zwanzig_TS_l_exp.hpp"
#include "time/SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk.hpp"
#include "time/SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk.hpp"

#include "time/ShackPDESWEPlaneMoriZwanzigTimeDiscretization.hpp"

/**
 * SWE Plane time steppers
 */
class PDESWEPlaneMoriZwanzigTimeSteppers
{
public:
	sweet::ErrorBase error;

	SWE_Plane_Mori_Zwanzig_TS_l_direct *l_direct = nullptr;

	PDESWEPlaneMoriZwanzigTS_BaseInterface *timestepper = nullptr;

	ShackPDESWEPlaneMoriZwanzigTimeDiscretization *shackPDESWEPlaneMoriZwanzigTimeDiscretization;

	bool linear_only = false;

	std::string equation;

	PDESWEPlaneMoriZwanzigTimeSteppers()
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
	}

	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackPDESWEPlaneMoriZwanzigTimeDiscretization = io_shackDict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzigTimeDiscretization>();

		SWE_Plane_Mori_Zwanzig_TS_l_direct dummy;
		dummy.shackRegistration(&io_shackDict);

		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(io_shackDict);
		return true;
	}

	bool setup(
			sweet::PlaneOperators *io_ops,
			sweet::ShackDictionary *io_shackDict,
			std::string i_equation
	)
	{
		assert(io_ops != nullptr);
		assert(io_shackDict != nullptr);

		equation = i_equation;
		if (equation != "P" && equation != "Q" && equation != "SF")
			SWEETError("Invalid equation: " + equation);

		// TODO: different tsm for each equation
		///const std::string &timestepping_method = shackPDESWEPlaneMoriZwanzigTimeDiscretization->timestepping_method;
		std::string timestepping_method;
		if (equation == "P")
			timestepping_method = shackPDESWEPlaneMoriZwanzigTimeDiscretization->timestepping_method_P;
		else
			timestepping_method = shackPDESWEPlaneMoriZwanzigTimeDiscretization->timestepping_method_Q;

		linear_only = false;

		if (timestepping_method == "l_exp_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneMoriZwanzigTS_BaseInterface*>(new SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk);
		}
		if (timestepping_method == "l_irk_n_erk")
		{
			timestepper = static_cast<PDESWEPlaneMoriZwanzigTS_BaseInterface*>(new SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk);
		}
		////else if (timestepping_method == "l_rexi")
		////{
		////	timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_rexi);

		////	linear_only = true;
		////}
		////else if (timestepping_method == "l_direct")
		////{
		////	timestepper = static_cast<PDESWEPlaneTS_BaseInterface*>(new SWE_Plane_TS_l_direct);
		////	linear_only = true;
		////}
		else //Help menu with list of schemes
		{
			std::cout << "Unknown method: " << timestepping_method << std::endl;
			std::cout << "Available --timestepping-method :"  << std::endl;
			////std::cout << "      l_direct       : Linear:     Analytical solution to linear SW operator"  << std::endl;
			////std::cout << "      l_rexi         : Linear:     Pure REXI, our little dog." << std::endl;
			std::cout << "      l_exp_n_erk   : Non-linear: Linear REXI, Non-linear RK, Strang-split" << std::endl;
			std::cout << "      l_irk_n_erk   : Non-linear: IRK, Non-linear ERK" << std::endl;

			SWEETError("No valid --timestepping-method provided");
		}

		timestepper->shackRegistration(io_shackDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timestepper);

		timestepper->setup(io_ops, equation);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*timestepper);

		return true;
	}

	~PDESWEPlaneMoriZwanzigTimeSteppers()
	{
		clear();
	}
};




#endif
