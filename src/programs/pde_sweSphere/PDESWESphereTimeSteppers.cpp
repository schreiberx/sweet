/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTimeSteppers.hpp"

#include "time/PDESWESphereTS_l_erk.hpp"
#include "time/PDESWESphereTS_l_erk_n_erk.hpp"
#include "time/PDESWESphereTS_l_erk_na_erk_uv.hpp"
#include "time/PDESWESphereTS_l_erk_na_erk_vd.hpp"
#include "time/PDESWESphereTS_l_exp.hpp"
#include "time/PDESWESphereTS_l_exp_direct_special.hpp"
#include "time/PDESWESphereTS_l_exp_n_erk.hpp"
#include "time/PDESWESphereTS_l_exp_n_etdrk.hpp"
#include "time/PDESWESphereTS_l_irk.hpp"
#include "time/PDESWESphereTS_l_irk_n_erk.hpp"
#include "time/PDESWESphereTS_l_irk_na_erk_uv.hpp"
#include "time/PDESWESphereTS_l_irk_na_erk_vd.hpp"
#include "time/PDESWESphereTS_l_irk_na_sl_nr_settls_uv_only.hpp"
#include "time/PDESWESphereTS_l_irk_na_sl_nr_settls_vd_only.hpp"
#include "time/PDESWESphereTS_l_irk_na_sl_settls_uv_only.hpp"
#include "time/PDESWESphereTS_l_irk_na_sl_settls_vd_only.hpp"
#include "time/PDESWESphereTS_lg_erk.hpp"
#include "time/PDESWESphereTS_lg_erk_lc_erk.hpp"
#include "time/PDESWESphereTS_lg_erk_lc_n_erk.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_erk.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_n_erk.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_n_etd_uv.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_n_etd_vd.hpp"
#include "time/PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv.hpp"
#include "time/PDESWESphereTS_lg_exp_na_sl_lc_nr_etdrk_uv.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_n_etdrk.hpp"
#include "time/PDESWESphereTS_lg_irk.hpp"
#include "time/PDESWESphereTS_lg_irk_lc_erk.hpp"
#include "time/PDESWESphereTS_lg_irk_lc_n_erk_ver01.hpp"
#include "time/PDESWESphereTS_lg_irk_lc_na_erk_vd.hpp"
#include "time/PDESWESphereTS_ln_erk.hpp"
#include "time/PDESWESphereTS_ln_erk_split_uv.hpp"
#include "time/PDESWESphereTS_ln_erk_split_vd.hpp"
#include "time/PDESWESphereTS_ln_settls_uv.hpp"
#include "time/PDESWESphereTS_ln_settls_vd.hpp"
#include "time/PDESWESphereTS_ln_sl_exp_settls_uv.hpp"
#include "time/PDESWESphereTS_ln_sl_exp_settls_vd.hpp"
#include "time/PDESWESphereTS_lg_0_lc_n_erk_bv.hpp"
#include "time/PDESWESphereTS_lg_exp_direct.hpp"
#include "time/PDESWESphereTS_lg_exp_lc_taylor.hpp"
#include "time/PDESWESphereTS_ln_imex_sdc.hpp"





PDESWESphereTimeSteppers::PDESWESphereTimeSteppers()
{
}




void PDESWESphereTimeSteppers::integrators_register_all(SphereOperators &i_op, sweet::ShackDictionary &i_shackDict)
{

	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk_lc_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_taylor(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_n_erk(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_n_erk(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_na_erk_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_na_erk_uv(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_erk_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_erk_uv(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_n_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_direct(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_na_erk_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_n_erk(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_erk(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk_lc_n_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk_split_uv(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk_split_vd(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_n_etdrk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etdrk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etd_uv(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etd_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_na_sl_lc_nr_etdrk_uv(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk(i_shackDict, i_op)));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_direct_special(i_shackDict, i_op)));

	/*
	 * IMEX SDC
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_imex_sdc(i_shackDict, i_op)));

	/*
	 * EXP SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_sl_exp_settls_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_sl_exp_settls_uv(i_shackDict, i_op)));

	/*
	 * ONLY SETTLS VERSION without any special variants
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_nr_settls_vd_only(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_nr_settls_uv_only(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_settls_vd_only(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_settls_uv_only(i_shackDict, i_op)));

	/*
	 * IRK SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_settls_vd(i_shackDict, i_op)));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_settls_uv(i_shackDict, i_op)));

	/*
	 * BAROTROPIC VORTICITY EQ
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_0_lc_n_erk_bv(i_shackDict, i_op)));
}



bool PDEAdvectionSphereTimeSteppers::setup_2_shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->shackRegistration(io_shackDict);
	}
	return true;
}



void PDEAdvectionSphereTimeSteppers::printImplementedTimesteppingMethods(
	std::ostream &o_ostream,
	const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	std::string prefix = i_prefix+"  ";
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->printImplementedTimesteppingMethods(o_ostream, prefix);
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << std::endl;
}


bool PDEAdvectionSphereTimeSteppers::setup_3_timestepper(
		const std::string &i_timestepping_method,
		sweet::ShackDictionary *io_shackDict,
		sweet::SphereOperators *io_ops
)
{
	if (i_timestepping_method == "")
	{
		printImplementedTimesteppingMethods();
		return error.set("Please set time stepping method using --timestepping-method=...");
	}
	/*
	 * Find right one
	 */
	timestepper = nullptr;

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->testImplementsTimesteppingMethod(i_timestepping_method))
		{
			if (timestepper != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				return error.set("Duplicate implementation for method "+i_timestepping_method);
			}

			//std::cout << "Found matching time stepping method at " << i+1 << "th element" << std::endl;
			ts->setup(io_ops);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(*ts);
			timestepper = ts;
		}
	}

	if (timestepper == nullptr)
		return error.set("No valid --timestepping-method '"+i_timestepping_method+"' provided");

	// Found integrator, freeing others
	_timesteppersFreeAll(timestepper);

	return true;
}


void PDEAdvectionSphereTimeSteppers::_timesteppersFreeAll(
		PDEAdvectionSphereTS_BaseInterface *i_skip_this_timestepper
)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == i_skip_this_timestepper)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}


void PDEAdvectionSphereTimeSteppers::clear()
{
	delete timestepper;
	timestepper = nullptr;

	_timesteppersFreeAll();
}


PDEAdvectionSphereTimeSteppers::~PDEAdvectionSphereTimeSteppers()
{
	clear();
}
