/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphere_TimeSteppers.hpp"

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





PDESWESphere_TimeSteppers::PDESWESphere_TimeSteppers()
{
}




void PDESWESphere_TimeSteppers::setup_1_registerAllTimesteppers()
{
	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_taylor));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_erk_na_erk_uv));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_erk_uv));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_n_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_direct));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_na_erk_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk_lc_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_erk));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk_lc_n_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk_split_uv));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_erk_split_vd));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_n_etdrk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etdrk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etd_uv));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_lc_n_etd_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_na_sl_lc_nr_etd_uv));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_exp_na_sl_lc_nr_etdrk_uv));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_irk));

	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_exp_direct_special));

	/*
	 * IMEX SDC
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_imex_sdc));

	/*
	 * EXP SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_sl_exp_settls_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_sl_exp_settls_uv));

	/*
	 * ONLY SETTLS VERSION without any special variants
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_nr_settls_vd_only));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_nr_settls_uv_only));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_settls_vd_only));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_l_irk_na_sl_settls_uv_only));

	/*
	 * IRK SETTLS VERSION
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_settls_vd));
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_ln_settls_uv));

	/*
	 * BAROTROPIC VORTICITY EQ
	 */
	_registered_integrators.push_back(static_cast<PDESWESphereTS_BaseInterface*>(new PDESWESphereTS_lg_0_lc_n_erk_bv));
}



bool PDESWESphere_TimeSteppers::setup_2_shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->shackRegistration(io_shackDict);
	}
	return true;
}



void PDESWESphere_TimeSteppers::printImplementedTimesteppingMethods(
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


bool PDESWESphere_TimeSteppers::setup_3_timestepper(
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
		PDESWESphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->implementsTimesteppingMethod(i_timestepping_method))
		{
			if (timestepper != nullptr)
			{
				return error.set("Duplicate implementation for method "+i_timestepping_method);
			}

			ts->setup_auto(i_timestepping_method, io_ops);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*ts);
			timestepper = ts;
		}
	}

	if (timestepper == nullptr)
		return error.set("No valid --timestepping-method '"+i_timestepping_method+"' provided");

	// Found integrator, freeing others
	_timesteppersFreeAll(timestepper);

	return true;
}


void PDESWESphere_TimeSteppers::_timesteppersFreeAll(
		PDESWESphereTS_BaseInterface *i_skip_this_timestepper
)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDESWESphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == i_skip_this_timestepper)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}


void PDESWESphere_TimeSteppers::clear()
{
	delete timestepper;
	timestepper = nullptr;

	_timesteppersFreeAll();
}


PDESWESphere_TimeSteppers::~PDESWESphere_TimeSteppers()
{
	clear();
}
