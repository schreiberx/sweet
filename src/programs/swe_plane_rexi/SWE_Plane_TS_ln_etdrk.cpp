/*
 * SWE_Plane_TS_l_phi0_n_edt.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_ln_etdrk.hpp"



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_ln_etdrk::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double &o_dt,
		double i_fixed_dt,
		double i_max_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

	o_dt = i_fixed_dt;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_ln_etdrk::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t	///< time updates
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

}


void SWE_Plane_TS_ln_etdrk::run_timestep(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("SWE_Plane_TS_l_phi0_n_edt: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;


	const PlaneDataConfig *planeDataConfig = io_h_pert.planeDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		PlaneData phi0_Un_h_pert(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				io_h_pert, io_u, io_v,
				phi0_Un_h_pert, phi0_Un_u, phi0_Un_v,
				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);

		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);
		euler_timestep_update_nonlinear(
				io_h_pert, io_u, io_v,
				FUn_h, FUn_u, FUn_v
		);

		PlaneData phi1_FUn_h_pert(planeDataConfig);
		PlaneData phi1_FUn_u(planeDataConfig);
		PlaneData phi1_FUn_v(planeDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h_pert, phi1_FUn_u, phi1_FUn_v,
				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);

		io_h_pert = phi0_Un_h_pert + i_fixed_dt*phi1_FUn_h_pert;
		io_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		io_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;
	}
	else if (timestepping_order == 2)
	{

		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		PlaneData phi0_Un_h_pert(planeDataConfig);
		PlaneData phi0_Un_u(planeDataConfig);
		PlaneData phi0_Un_v(planeDataConfig);
		ts_phi0_rexi.run_timestep(
				io_h_pert, io_u, io_v,
				phi0_Un_h_pert, phi0_Un_u, phi0_Un_v,
				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);

		PlaneData FUn_h(planeDataConfig);
		PlaneData FUn_u(planeDataConfig);
		PlaneData FUn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				io_h_pert, io_u, io_v,
				FUn_h, FUn_u, FUn_v
		);

		PlaneData phi1_FUn_h_pert(planeDataConfig);
		PlaneData phi1_FUn_u(planeDataConfig);
		PlaneData phi1_FUn_v(planeDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h_pert, phi1_FUn_u, phi1_FUn_v,
				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);

		PlaneData A_h_pert = phi0_Un_h_pert + i_fixed_dt*phi1_FUn_h_pert;
		PlaneData A_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		PlaneData A_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		PlaneData FAn_h(planeDataConfig);
		PlaneData FAn_u(planeDataConfig);
		PlaneData FAn_v(planeDataConfig);

		euler_timestep_update_nonlinear(
				A_h_pert, A_u, A_v,
				FAn_h, FAn_u, FAn_v
		);


		PlaneData phi2_X_h_pert(planeDataConfig);
		PlaneData phi2_X_u(planeDataConfig);
		PlaneData phi2_X_v(planeDataConfig);

		ts_phi2_rexi.run_timestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h_pert,
				phi2_X_u,
				phi2_X_v,

				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);

		io_h_pert.physical_set_all(1);
		io_u.physical_set_all(1);
		io_v.physical_set_all(1);

		io_h_pert = A_h_pert + i_fixed_dt*phi2_X_h_pert;
		io_u = A_u + i_fixed_dt*phi2_X_u;
		io_v = A_v + i_fixed_dt*phi2_X_v;

	}
	else
	{
		FatalError("TODO: This order is not implemented, yet!");
	}

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_ln_etdrk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order
)
{
	ts_phi0_rexi.setup(i_rexiSimVars, "phi0");
	ts_phi1_rexi.setup(i_rexiSimVars, "phi1");
	ts_phi2_rexi.setup(i_rexiSimVars, "phi2");

	timestepping_order = i_timestepping_order;
}


SWE_Plane_TS_ln_etdrk::SWE_Plane_TS_ln_etdrk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_phi0_rexi(simVars, op),
		ts_phi1_rexi(simVars, op),
		ts_phi2_rexi(simVars, op)
{
}



SWE_Plane_TS_ln_etdrk::~SWE_Plane_TS_ln_etdrk()
{
}

