/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *         
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_0_LF_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_0_LF_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_lg_erk_lc_n_erk.hpp"
#include "PDESWESphereTS_lg_irk.hpp"



class PDESWESphereTS_lg_0_lc_n_erk_bv	: public PDESWESphereTS_BaseInterface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		timestepping_order2 = shackDict.disc.timestepping_order2;
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
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int timestepping_order;

	double timestep_size;

	/*
	 * Non-linear time steppers
	 */
	SphereTimestepping_ExplicitRK timestepping_rk;

public:
	PDESWESphereTS_lg_0_lc_n_erk_bv(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void euler_timestep_update(
		const sweet::SphereData_Spectral &i_phi, //prog
		const sweet::SphereData_Spectral &i_vrt, //prog
		const sweet::SphereData_Spectral &i_div, //prog

		sweet::SphereData_Spectral &o_phi_t, //updated with euler
		sweet::SphereData_Spectral &o_vrt_t, //updated with euler
		sweet::SphereData_Spectral &o_div_t, //updated with euler

		double i_simulation_timestamp
	);

	virtual ~PDESWESphereTS_lg_0_lc_n_erk_bv();
};

#endif 
