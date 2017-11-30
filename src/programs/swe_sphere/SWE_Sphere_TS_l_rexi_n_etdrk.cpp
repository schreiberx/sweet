/*
 * SWE_Sphere_TS_ln_edtrk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Sphere_TS_l_rexi_n_etdrk.hpp"




void SWE_Sphere_TS_l_rexi_n_etdrk::run_timestep(
		SphereData &io_phi,	///< prognostic variables
		SphereData &io_u,	///< prognostic variables
		SphereData &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_phi0_n_edt: Only constant time step size allowed");

	const SphereDataConfig *sphereDataConfig = io_phi.sphereDataConfig;

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		SphereData phi0_Un_h(sphereDataConfig);
		SphereData phi0_Un_u(sphereDataConfig);
		SphereData phi0_Un_v(sphereDataConfig);
		ts_phi0_rexi.run_timestep(
				io_phi, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		SphereData FUn_h(sphereDataConfig);
		SphereData FUn_u(sphereDataConfig);
		SphereData FUn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		SphereData phi1_FUn_h(sphereDataConfig);
		SphereData phi1_FUn_u(sphereDataConfig);
		SphereData phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		io_phi = phi0_Un_h + i_dt*phi1_FUn_h;
		io_u = phi0_Un_u + i_dt*phi1_FUn_u;
		io_v = phi0_Un_v + i_dt*phi1_FUn_v;
	}
	else if (timestepping_order == 2)
	{

		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		SphereData phi0_Un_h(sphereDataConfig);
		SphereData phi0_Un_u(sphereDataConfig);
		SphereData phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep(
				io_phi, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_dt,
				i_simulation_timestamp
			);

		SphereData FUn_h(sphereDataConfig);
		SphereData FUn_u(sphereDataConfig);
		SphereData FUn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		SphereData phi1_FUn_h(sphereDataConfig);
		SphereData phi1_FUn_u(sphereDataConfig);
		SphereData phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_dt,
				i_simulation_timestamp
			);

		SphereData A_h = phi0_Un_h + i_dt*phi1_FUn_h;
		SphereData A_u = phi0_Un_u + i_dt*phi1_FUn_u;
		SphereData A_v = phi0_Un_v + i_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		SphereData FAn_h(sphereDataConfig);
		SphereData FAn_u(sphereDataConfig);
		SphereData FAn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		SphereData phi2_X_h(sphereDataConfig);
		SphereData phi2_X_u(sphereDataConfig);
		SphereData phi2_X_v(sphereDataConfig);

		ts_phi2_rexi.run_timestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_dt,
				i_simulation_timestamp
			);

		io_phi = A_h + i_dt*phi2_X_h;
		io_u = A_u + i_dt*phi2_X_u;
		io_v = A_v + i_dt*phi2_X_v;
	}
	else if (timestepping_order == 4)
	{
		double dt = i_dt;
		double dt_half = dt*0.5;



		/*
		 * Precompute commonly used terms
		 */
		SphereData phi0_Un_h(sphereDataConfig);
		SphereData phi0_Un_u(sphereDataConfig);
		SphereData phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep(
				io_phi, io_u, io_v,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData FUn_h(sphereDataConfig);
		SphereData FUn_u(sphereDataConfig);
		SphereData FUn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				io_phi, io_u, io_v,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);



		/*
		 * Some commonly shared buffers
		 */

		SphereData phi1_h(sphereDataConfig);
		SphereData phi1_u(sphereDataConfig);
		SphereData phi1_v(sphereDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.run_timestep(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData A_h = phi0_Un_h + dt_half*phi1_h;
		SphereData A_u = phi0_Un_u + dt_half*phi1_u;
		SphereData A_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */

		SphereData FAn_h(sphereDataConfig);
		SphereData FAn_u(sphereDataConfig);
		SphereData FAn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData B_h = phi0_Un_h + dt_half*phi1_h;
		SphereData B_u = phi0_Un_u + dt_half*phi1_u;
		SphereData B_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */

		SphereData phi0_An_h(sphereDataConfig);
		SphereData phi0_An_u(sphereDataConfig);
		SphereData phi0_An_v(sphereDataConfig);

		ts_phi0_rexi.run_timestep(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);


		SphereData FBn_h(sphereDataConfig);
		SphereData FBn_u(sphereDataConfig);
		SphereData FBn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.run_timestep(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		SphereData C_h = phi0_An_h + dt_half*phi1_h;
		SphereData C_u = phi0_An_u + dt_half*phi1_u;
		SphereData C_v = phi0_An_v + dt_half*phi1_v;



		/*
		 * R0 - R3
		 */
		SphereData FCn_h(sphereDataConfig);
		SphereData FCn_u(sphereDataConfig);
		SphereData FCn_v(sphereDataConfig);

		ts_l_erk_n_erk.euler_timestep_update_nonlinear(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);

		SphereData R0_h = io_phi;
		SphereData R0_u = io_u;
		SphereData R0_v = io_v;

		SphereData &R1_h = FUn_h;
		SphereData &R1_u = FUn_u;
		SphereData &R1_v = FUn_v;

		SphereData R2_h = FAn_h + FBn_h;
		SphereData R2_u = FAn_u + FBn_u;
		SphereData R2_v = FAn_v + FBn_v;

		SphereData &R3_h = FCn_h;
		SphereData &R3_u = FCn_u;
		SphereData &R3_v = FCn_v;


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
		ts_ups0_rexi.run_timestep(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

		ts_ups1_rexi.run_timestep(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);

		ts_ups2_rexi.run_timestep(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);

		ts_ups3_rexi.run_timestep(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);

		io_phi = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_u = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_v = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
	}
	else
	{
		FatalError("TODO: This order is not implemented, yet!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_rexi_n_etdrk::setup(
		REXI_SimulationVariables &i_rexiSimVars,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size,
		bool i_use_f_sphere
)
{
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	ts_l_erk_n_erk.setup(timestepping_order, timestepping_order2);

	if (timestepping_order != timestepping_order2)
		FatalError("Mismatch of orders, should be equal");

	if (timestepping_order == 0)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);	/* set use_f_sphere to true */
	}
	else if (timestepping_order -= 2)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);	/* set use_f_sphere to true */
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", i_timestep_size, false, true);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", i_timestep_size, false, true);
	}
	else if  (timestepping_order >= 4)
	{
		ts_phi0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size*0.5, false, true);	/* set use_f_sphere to true */
		ts_phi1_rexi.setup(i_rexiSimVars, "phi1", i_timestep_size*0.5, false, true);
		ts_phi2_rexi.setup(i_rexiSimVars, "phi2", i_timestep_size*0.5, false, true);

		// phi0, but with a full time step size
		ts_ups0_rexi.setup(i_rexiSimVars, "phi0", i_timestep_size, false, true);
		ts_ups1_rexi.setup(i_rexiSimVars, "ups1", i_timestep_size, false, true);
		ts_ups2_rexi.setup(i_rexiSimVars, "ups2", i_timestep_size, false, true);
		ts_ups3_rexi.setup(i_rexiSimVars, "ups3", i_timestep_size, false, true);
	}
}


SWE_Sphere_TS_l_rexi_n_etdrk::SWE_Sphere_TS_l_rexi_n_etdrk(
		SimulationVariables &i_simVars,
		SphereOperators &i_op
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



SWE_Sphere_TS_l_rexi_n_etdrk::~SWE_Sphere_TS_l_rexi_n_etdrk()
{
}

