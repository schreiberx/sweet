/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"



class PDESWESphereTS_l_erk_n_erk	: public PDESWESphereTS_BaseInterface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);
	std::string string_id();

	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int timestepping_order;
	int timestepping_order2;

	sweet::TimesteppingExplicitRKSphereData timestepping_rk_linear;
	sweet::TimesteppingExplicitRKSphereData timestepping_rk_nonlinear;


public:
	void euler_timestep_update_linear(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			const sweet::SphereData_Spectral &i_h,	///< prognostic variables
			const sweet::SphereData_Spectral &i_u,	///< prognostic variables
			const sweet::SphereData_Spectral &i_v,	///< prognostic variables

			sweet::SphereData_Spectral &o_h_t,	///< time updates
			sweet::SphereData_Spectral &o_u_t,	///< time updates
			sweet::SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

public:
	void euler_timestep_update(
			const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
			const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
			const sweet::SphereData_Spectral &i_div,	///< prognostic variables

			sweet::SphereData_Spectral &o_phi_t,	///< time updates
			sweet::SphereData_Spectral &o_vort_t,	///< time updates
			sweet::SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	PDESWESphereTS_l_erk_n_erk(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2
	);


	void setup_auto();

	void run_timestep(
			sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);




	virtual ~PDESWESphereTS_l_erk_n_erk();
};

#endif
