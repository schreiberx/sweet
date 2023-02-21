/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 *         
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_0_LF_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_0_LF_N_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"



class SWE_Sphere_TS_lg_0_lc_n_erk_bv	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = simVars.disc.timestepping_order;
		timestepping_order2 = simVars.disc.timestepping_order2;
		if (
			i_timestepping_method == "lg_0_lc_n_erk_bv" 
		)
			return true;

		return false;
	}

	std::string string_id()
	{
		std::string s = "lg_0_lc_n_erk_bv";
		return s;
	}

	void setup_auto();
	void print_help();


private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;

	double timestep_size;

	/*
	 * Non-linear time steppers
	 */
	SphereTimestepping_ExplicitRK timestepping_rk;

public:
	SWE_Sphere_TS_lg_0_lc_n_erk_bv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void euler_timestep_update(
		const SphereData_Spectral &i_phi, //prog
		const SphereData_Spectral &i_vrt, //prog
		const SphereData_Spectral &i_div, //prog

		SphereData_Spectral &o_phi_t, //updated with euler
		SphereData_Spectral &o_vrt_t, //updated with euler
		SphereData_Spectral &o_div_t, //updated with euler

		double i_simulation_timestamp
	);

	virtual ~SWE_Sphere_TS_lg_0_lc_n_erk_bv();
};

#endif 
