/*
 * PDESWESphereTS_split_lg_lc_na_nr_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_SPLIT_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_ERK_SPLIT_VD_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"


class PDESWESphereTS_ln_erk_split_vd	: public PDESWESphereTS_BaseInterface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		//timestepping_order2 = shackDict.disc.timestepping_order2;
		if (
				i_timestepping_method == "l_na_erk_split_vd"	||
				i_timestepping_method == "l_na_erk_split_aa_vd"	||
				i_timestepping_method == "l_erk_split_vd"	||
				i_timestepping_method == "l_erk_split_aa_vd"	||
				i_timestepping_method == "ln_erk_split_vd"		||
				i_timestepping_method == "ln_erk_split_aa_vd"
		)
			return true;

		return false;
	}

	std::string string_id()
	{
		return "ln_erk_split_vd";
	}

	void setup_auto()
	{
		/*
		 * l_na
		 */
		if (timestepping_method == "l_na_erk_split_vd")
		{
			setup(timestepping_order, true, true, true, false, false);
			return;
		}

		if (timestepping_method == "l_na_erk_split_aa_vd")
		{
			setup(timestepping_order, true, true, true, false, true);
			return;
		}

		/*
		 * l
		 */
		if (timestepping_method == "l_erk_split_vd")
		{
			setup(timestepping_order, true, true, false, false, false);
			return;
		}

		if (timestepping_method == "l_erk_split_aa_vd")
		{
			setup(timestepping_order, true, true, false, false, true);
			return;
		}

		/*
		 * ln
		 */
		if (timestepping_method == "ln_erk_split_vd")
		{
			setup(timestepping_order, true, true, true, true, false);
			return;
		}

		if (timestepping_method == "ln_erk_split_aa_vd")
		{
			setup(timestepping_order, true, true, true, true, true);
			return;
		}

		SWEETError("Should never happen");
	}

private:
	sweet::ShackDictionary &shackDict;
	
	sweet::SphereOperators &op;

	int timestepping_order;

	bool use_lg = false;
	bool use_lc = false;
	bool use_na = false;
	bool use_nr = false;

	bool anti_aliasing_for_each_term = false;


	// Sampler
	SphereTimestepping_ExplicitRK timestepping_rk;


public:
	void euler_timestep_update_lg(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_pert_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_lc(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_lc_spectral_only(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_na(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_update_nr(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);



public:
	void euler_timestep_set_tendencies(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_set_tendencies_na_only(
			const sweet::SphereData_Spectral &i_U_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_vrt,	///< prognostic variables
			const sweet::SphereData_Spectral &i_U_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_U_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_U_vrt_t,	///< time updates
			sweet::SphereData_Spectral &o_U_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphereTS_ln_erk_split_vd(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			bool i_lg,
			bool i_lc,
			bool i_na,
			bool i_nr,
			bool i_antialiasing_for_each_term
	);

	void run_timestep(
			sweet::SphereData_Spectral &io_U_phi,
			sweet::SphereData_Spectral &io_U_vrt,
			sweet::SphereData_Spectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_na(
			sweet::SphereData_Spectral &io_U_phi,
			sweet::SphereData_Spectral &io_U_vrt,
			sweet::SphereData_Spectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);




	virtual ~PDESWESphereTS_ln_erk_split_vd();
};

#endif
