#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_NA_SL_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_NA_SL_HPP_


#include <limits>

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_l_erk_na_sl	: public SWE_Sphere_TS_interface
{
public:
	static bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
		return i_timestepping_method == "l_erk_na_sl";
	}

public:
	std::string string_id()
	{
		return "l_erk_na_sl";
	}

	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;
	int timestepping_order2;

	SphereTimestepping_ExplicitRK timestepping_rk_linear;
	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	SphereData_Spectral U_phi_pert_prev, U_vrt_prev, U_div_prev;

public:
	void euler_timestep_update_linear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);



public:
	void euler_timestep_update_n(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


	void euler_timestep_update_n_erk_stage(
			const SphereData_Spectral &i_phi,	///< prognostic variables
			const SphereData_Spectral &i_vort,	///< prognostic variables
			const SphereData_Spectral &i_div,	///< prognostic variables

			SphereData_Spectral &o_phi_dt,	///< time updates
			SphereData_Spectral &o_vort_dt,	///< time updates
			SphereData_Spectral &o_div_dt,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_n(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


public:
	void euler_timestep_update(
			const SphereData_Spectral &i_phi,	///< prognostic variables
			const SphereData_Spectral &i_vort,	///< prognostic variables
			const SphereData_Spectral &i_div,	///< prognostic variables

			SphereData_Spectral &o_phi_t,	///< time updates
			SphereData_Spectral &o_vort_t,	///< time updates
			SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);


public:
	SWE_Sphere_TS_l_erk_na_sl(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);


	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2
	);


	void setup_auto();

	void run_timestep_pert(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


	void run_timestep_nonpert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_erk_na_sl();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
