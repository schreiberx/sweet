/*
 * SWE_Sphere_TS_ln_edtrk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_exp_n_etdrk.hpp"



bool SWE_Sphere_TS_l_exp_n_etdrk::implements_timestepping_method(const std::string &i_timestepping_method)
{
	if (i_timestepping_method == "l_exp_n_etdrk")
		return true;

	return false;
}

std::string SWE_Sphere_TS_l_exp_n_etdrk::string_id()
{
	return "l_exp_n_etdrk";
}

void SWE_Sphere_TS_l_exp_n_etdrk::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	setup(
			simVars.rexi,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			simVars.timecontrol.current_timestep_size
		);
}



void SWE_Sphere_TS_l_exp_n_etdrk::run_timestep_pert(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_phi.sphereDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		SphereData_Spectral phi0_Un_phi(sphereDataConfig), phi0_Un_vrt(sphereDataConfig), phi0_Un_div(sphereDataConfig);
		ts_phi0_rexi.run_timestep_pert(
				io_phi, io_vrt, io_div,
				phi0_Un_phi, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_phi(sphereDataConfig), FUn_vrt(sphereDataConfig), FUn_div(sphereDataConfig);
		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_vrt, io_div,
				FUn_phi, FUn_vrt, FUn_div,
				i_simulation_timestamp
			);

		SphereData_Spectral phi1_FUn_phi(sphereDataConfig), phi1_FUn_vrt(sphereDataConfig), phi1_FUn_div(sphereDataConfig);
		ts_phi1_rexi.run_timestep_pert(
				FUn_phi, FUn_vrt, FUn_div,
				phi1_FUn_phi, phi1_FUn_vrt, phi1_FUn_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi = phi0_Un_phi + i_fixed_dt*phi1_FUn_phi;
		io_vrt = phi0_Un_vrt + i_fixed_dt*phi1_FUn_vrt;
		io_div = phi0_Un_div + i_fixed_dt*phi1_FUn_div;
		return;
	}

	if (timestepping_order == 2)
	{
		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		SphereData_Spectral phi0_Un_h(sphereDataConfig), phi0_Un_u(sphereDataConfig), phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep_pert(
				io_phi, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_h(sphereDataConfig), FUn_u(sphereDataConfig), FUn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		SphereData_Spectral phi1_FUn_h(sphereDataConfig), phi1_FUn_u(sphereDataConfig), phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.run_timestep_pert(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral A_h = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		SphereData_Spectral A_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		SphereData_Spectral A_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		SphereData_Spectral FAn_h(sphereDataConfig), FAn_u(sphereDataConfig), FAn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		SphereData_Spectral phi2_X_h(sphereDataConfig), phi2_X_u(sphereDataConfig), phi2_X_v(sphereDataConfig);

		ts_phi2_rexi.run_timestep_pert(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi = A_h + i_fixed_dt*phi2_X_h;
		io_vrt = A_u + i_fixed_dt*phi2_X_u;
		io_div = A_v + i_fixed_dt*phi2_X_v;
		return;
	}

	if (timestepping_order == 4)
	{
		double dt = i_fixed_dt;
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */
		SphereData_Spectral phi0_Un_h(sphereDataConfig), phi0_Un_u(sphereDataConfig), phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep_pert(
				io_phi, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_h(sphereDataConfig), FUn_u(sphereDataConfig), FUn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);



		/*
		 * Some commonly shared buffers
		 */

		SphereData_Spectral phi1_h(sphereDataConfig), phi1_u(sphereDataConfig), phi1_v(sphereDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.run_timestep_pert(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral A_h = phi0_Un_h + dt_half*phi1_h;
		SphereData_Spectral A_u = phi0_Un_u + dt_half*phi1_u;
		SphereData_Spectral A_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */

		SphereData_Spectral FAn_h(sphereDataConfig), FAn_u(sphereDataConfig), FAn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep_pert(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral B_h = phi0_Un_h + dt_half*phi1_h;
		SphereData_Spectral B_u = phi0_Un_u + dt_half*phi1_u;
		SphereData_Spectral B_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */

		SphereData_Spectral phi0_An_h(sphereDataConfig), phi0_An_u(sphereDataConfig), phi0_An_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep_pert(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);


		SphereData_Spectral FBn_h(sphereDataConfig), FBn_u(sphereDataConfig), FBn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep_pert(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral C_h = phi0_An_h + dt_half*phi1_h;
		SphereData_Spectral C_u = phi0_An_u + dt_half*phi1_u;
		SphereData_Spectral C_v = phi0_An_v + dt_half*phi1_v;



		/*
		 * R0 - R3
		 */
		SphereData_Spectral FCn_h(sphereDataConfig), FCn_u(sphereDataConfig), FCn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);

		SphereData_Spectral R0_h = io_phi;
		SphereData_Spectral R0_u = io_vrt;
		SphereData_Spectral R0_v = io_div;

		SphereData_Spectral &R1_h = FUn_h;
		SphereData_Spectral &R1_u = FUn_u;
		SphereData_Spectral &R1_v = FUn_v;

		SphereData_Spectral R2_h = FAn_h + FBn_h;
		SphereData_Spectral R2_u = FAn_u + FBn_u;
		SphereData_Spectral R2_v = FAn_v + FBn_v;

		SphereData_Spectral &R3_h = FCn_h;
		SphereData_Spectral &R3_u = FCn_u;
		SphereData_Spectral &R3_v = FCn_v;


		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL)R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */
		ts_ups0_rexi.run_timestep_pert(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

		ts_ups1_rexi.run_timestep_pert(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);

		ts_ups2_rexi.run_timestep_pert(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);

		ts_ups3_rexi.run_timestep_pert(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);

		io_phi = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_vrt = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_div = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
		return;
	}

	SWEETError("TODO: This order is not implemented, yet!");
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_exp_n_etdrk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size
)
{
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	if (timestepping_order != timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	ts_l_erk_n_erk.setup(timestepping_order, timestepping_order2);

	if (timestepping_order == 1)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, false);
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", i_timestep_size, false, false);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, false);
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", i_timestep_size, false, false);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", i_timestep_size, false, false);
	}
	else if  (timestepping_order == 4)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size*0.5, false, false);
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", i_timestep_size*0.5, false, false);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", i_timestep_size*0.5, false, false);

		// phi0, but with a full time step size
		ts_ups0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, false);
		ts_ups1_rexi.setup(i_rexiSimVars, "ups1", i_timestep_size, false, false);
		ts_ups2_rexi.setup(i_rexiSimVars, "ups2", i_timestep_size, false, false);
		ts_ups3_rexi.setup(i_rexiSimVars, "ups3", i_timestep_size, false, false);
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}


SWE_Sphere_TS_l_exp_n_etdrk::SWE_Sphere_TS_l_exp_n_etdrk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		ts_l_erk_n_erk(simVars, op),

		ts_phi0_rexi(simVars, op),
		ts_phi1_rexi(simVars, op),
		ts_phi2_rexi(simVars, op),

		ts_ups0_rexi(simVars, op),
		ts_ups1_rexi(simVars, op),
		ts_ups2_rexi(simVars, op),
		ts_ups3_rexi(simVars, op)
{
}



SWE_Sphere_TS_l_exp_n_etdrk::~SWE_Sphere_TS_l_exp_n_etdrk()
{
}

