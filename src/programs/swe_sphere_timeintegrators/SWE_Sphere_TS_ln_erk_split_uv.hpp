/*
 * SWE_Sphere_TS_split_lg_lc_na_nr_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LN_ERK_SPLIT_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LN_ERK_SPLIT_UV_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"


class SWE_Sphere_TS_ln_erk_split_uv	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = simVars.disc.timestepping_order;
		//timestepping_order2 = simVars.disc.timestepping_order2;
		if (
				i_timestepping_method == "l_na_erk_split_uv"	||
				i_timestepping_method == "l_na_erk_split_aa_uv"	||
				i_timestepping_method == "l_erk_split_uv"	||
				i_timestepping_method == "l_erk_split_aa_uv"	||
				i_timestepping_method == "ln_erk_split_uv"		||
				i_timestepping_method == "ln_erk_split_aa_uv"
		)
			return true;

		return false;
	}

	std::string string_id()
	{
		return "ln_erk_split_uv";
	}

	void setup_auto()
	{
		/*
		 * l_na
		 */
		if (timestepping_method == "l_na_erk_split_uv")
		{
			setup(timestepping_order, true, true, true, false, false);
			return;
		}

		if (timestepping_method == "l_na_erk_split_aa_uv")
		{
			setup(timestepping_order, true, true, true, false, true);
			return;
		}

		/*
		 * l
		 */
		if (timestepping_method == "l_erk_split_uv")
		{
			setup(timestepping_order, true, true, false, false, false);
			return;
		}

		if (timestepping_method == "l_erk_split_aa_uv")
		{
			setup(timestepping_order, true, true, false, false, true);
			return;
		}

		/*
		 * ln
		 */
		if (timestepping_method == "ln_erk_split_uv")
		{
			setup(timestepping_order, true, true, true, true, false);
			return;
		}

		if (timestepping_method == "ln_erk_split_aa_uv")
		{
			setup(timestepping_order, true, true, true, true, true);
			return;
		}

		SWEETError("Should never happen");
	}

private:
	SimulationVariables &simVars;

	SphereOperators_SphereData &op;

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
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_lc(
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_na(
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


public:
	void euler_timestep_update_nr(
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


private:
	void euler_timestep_set_tendencies(
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


private:
	void euler_timestep_set_tendencies_na_only(
			const SphereData_Spectral &i_U_phi,
			const SphereData_Spectral &i_U_vrt,
			const SphereData_Spectral &i_U_div,

			SphereData_Spectral &o_U_phi_t,
			SphereData_Spectral &o_U_vrt_t,
			SphereData_Spectral &o_U_div_t,

			double i_simulation_timestamp = -1
	);


public:
	SWE_Sphere_TS_ln_erk_split_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
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
			SphereData_Spectral &io_U_phi,
			SphereData_Spectral &io_U_vrt,
			SphereData_Spectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_na(
			SphereData_Spectral &io_U_phi,
			SphereData_Spectral &io_U_vrt,
			SphereData_Spectral &io_U_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);




	virtual ~SWE_Sphere_TS_ln_erk_split_uv();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
