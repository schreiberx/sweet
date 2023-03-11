/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_exp_lc_n_etdrk.hpp"


bool PDESWESphereTS_lg_exp_lc_n_etdrk::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "lg_exp_lc_n_etdrk")
		return true;

	return false;
}

bool PDESWESphereTS_lg_exp_lc_n_etdrk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere->sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	return setup_main(
			io_ops,
			shackExpIntegration,
			shackPDESWETimeDisc->timestepping_order,
			shackPDESWETimeDisc->timestepping_order2,
			shackTimestepControl->current_timestep_size
		);
}


bool PDESWESphereTS_lg_exp_lc_n_etdrk::setup_main(
		sweet::SphereOperators *io_ops,
		sweet::ShackExpIntegration *i_shackExpIntegration,
		int i_timestepping_order,
		int i_timestepping_order2,
		double i_timestep_size
)
{
	ops = io_ops;
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;

	ts_lg_erk_lc_n_erk.setup(ops, i_timestepping_order, -1);

	if (timestepping_order != timestepping_order2)
		SWEETError("Mismatch of orders, should be equal");

	if (timestepping_order == 0 || timestepping_order == 1)
	{
		ts_phi0_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestep_size, false, true, timestepping_order);
	}
	else if (timestepping_order == 2)
	{
		ts_phi0_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestep_size, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestep_size, false, true, timestepping_order);
		ts_phi2_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestep_size, false, true, timestepping_order);
	}
	else if  (timestepping_order == 4)
	{
		ts_phi0_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestep_size*0.5, false, true, timestepping_order);	/* NO Coriolis */
		ts_phi1_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi1", i_timestep_size*0.5, false, true, timestepping_order);
		ts_phi2_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi2", i_timestep_size*0.5, false, true, timestepping_order);

		// phi0, but with a full time step size
		ts_ups0_rexi.setup_variant_50(ops, i_shackExpIntegration, "phi0", i_timestep_size, false, true, timestepping_order);
		ts_ups1_rexi.setup_variant_50(ops, i_shackExpIntegration, "ups1", i_timestep_size, false, true, timestepping_order);
		ts_ups2_rexi.setup_variant_50(ops, i_shackExpIntegration, "ups2", i_timestep_size, false, true, timestepping_order);
		ts_ups3_rexi.setup_variant_50(ops, i_shackExpIntegration, "ups3", i_timestep_size, false, true, timestepping_order);
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}

	return true;
}


std::string PDESWESphereTS_lg_exp_lc_n_etdrk::getIDString()
{
	return "lg_exp_lc_n_etdrk";
}

void PDESWESphereTS_lg_exp_lc_n_etdrk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const sweet::SphereData_Config *sphereDataConfig = io_phi_pert.sphereDataConfig;

	if (timestepping_order == 0 || timestepping_order == 1)
	{
		/*
		 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
		 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);
		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		sweet::SphereData_Spectral phi1_FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi_pert = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		io_vrt = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		io_div = phi0_Un_v + i_fixed_dt*phi1_FUn_v;
	}
	else if (timestepping_order == 2)
	{
		/*
		 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
		 */

		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		sweet::SphereData_Spectral phi1_FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_FUn_v(sphereDataConfig);

		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_FUn_h, phi1_FUn_u, phi1_FUn_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral A_h = phi0_Un_h + i_fixed_dt*phi1_FUn_h;
		sweet::SphereData_Spectral A_u = phi0_Un_u + i_fixed_dt*phi1_FUn_u;
		sweet::SphereData_Spectral A_v = phi0_Un_v + i_fixed_dt*phi1_FUn_v;

		/*
		 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
		 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
		 */

		sweet::SphereData_Spectral FAn_h(sphereDataConfig);
		sweet::SphereData_Spectral FAn_u(sphereDataConfig);
		sweet::SphereData_Spectral FAn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		sweet::SphereData_Spectral phi2_X_h(sphereDataConfig);
		sweet::SphereData_Spectral phi2_X_u(sphereDataConfig);
		sweet::SphereData_Spectral phi2_X_v(sphereDataConfig);

		ts_phi2_rexi.runTimestep(
				FAn_h - FUn_h,
				FAn_u - FUn_u,
				FAn_v - FUn_v,

				phi2_X_h,
				phi2_X_u,
				phi2_X_v,

				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi_pert = A_h + i_fixed_dt*phi2_X_h;
		io_vrt = A_u + i_fixed_dt*phi2_X_u;
		io_div = A_v + i_fixed_dt*phi2_X_v;
	}
	else if (timestepping_order == 4)
	{
		double dt = i_fixed_dt;
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */
		sweet::SphereData_Spectral phi0_Un_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_Un_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral FUn_h(sphereDataConfig);
		sweet::SphereData_Spectral FUn_u(sphereDataConfig);
		sweet::SphereData_Spectral FUn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);



		/*
		 * Some commonly shared buffers
		 */

		sweet::SphereData_Spectral phi1_h(sphereDataConfig);
		sweet::SphereData_Spectral phi1_u(sphereDataConfig);
		sweet::SphereData_Spectral phi1_v(sphereDataConfig);


		/*
		 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
		 */
		ts_phi1_rexi.runTimestep(
				FUn_h, FUn_u, FUn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral A_h = phi0_Un_h + dt_half*phi1_h;
		sweet::SphereData_Spectral A_u = phi0_Un_u + dt_half*phi1_u;
		sweet::SphereData_Spectral A_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
		 */

		sweet::SphereData_Spectral FAn_h(sphereDataConfig);
		sweet::SphereData_Spectral FAn_u(sphereDataConfig);
		sweet::SphereData_Spectral FAn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.runTimestep(
				FAn_h, FAn_u, FAn_v,
				phi1_h, phi1_u, phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral B_h = phi0_Un_h + dt_half*phi1_h;
		sweet::SphereData_Spectral B_u = phi0_Un_u + dt_half*phi1_u;
		sweet::SphereData_Spectral B_v = phi0_Un_v + dt_half*phi1_v;



		/*
		 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
		 */

		sweet::SphereData_Spectral phi0_An_h(sphereDataConfig);
		sweet::SphereData_Spectral phi0_An_u(sphereDataConfig);
		sweet::SphereData_Spectral phi0_An_v(sphereDataConfig);

		ts_phi0_rexi.runTimestep(
				A_h, A_u, A_v,
				phi0_An_h, phi0_An_u, phi0_An_v,
				dt_half,
				i_simulation_timestamp
			);


		sweet::SphereData_Spectral FBn_h(sphereDataConfig);
		sweet::SphereData_Spectral FBn_u(sphereDataConfig);
		sweet::SphereData_Spectral FBn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				B_h, B_u, B_v,
				FBn_h, FBn_u, FBn_v,
				i_simulation_timestamp + dt_half
		);

		ts_phi1_rexi.runTimestep(
				2.0*FBn_h - FUn_h,
				2.0*FBn_u - FUn_u,
				2.0*FBn_v - FUn_v,
				phi1_h,	phi1_u,	phi1_v,
				dt_half,
				i_simulation_timestamp
			);

		sweet::SphereData_Spectral C_h = phi0_An_h + dt_half*phi1_h;
		sweet::SphereData_Spectral C_u = phi0_An_u + dt_half*phi1_u;
		sweet::SphereData_Spectral C_v = phi0_An_v + dt_half*phi1_v;



		/*
		 * R0 - R3
		 */
		sweet::SphereData_Spectral FCn_h(sphereDataConfig);
		sweet::SphereData_Spectral FCn_u(sphereDataConfig);
		sweet::SphereData_Spectral FCn_v(sphereDataConfig);

		ts_lg_erk_lc_n_erk.euler_timestep_update_lc_n(
				C_h, C_u, C_v,
				FCn_h, FCn_u, FCn_v,
				i_simulation_timestamp + dt
		);

		sweet::SphereData_Spectral R0_h = io_phi_pert;
		sweet::SphereData_Spectral R0_u = io_vrt;
		sweet::SphereData_Spectral R0_v = io_div;

		sweet::SphereData_Spectral &R1_h = FUn_h;
		sweet::SphereData_Spectral &R1_u = FUn_u;
		sweet::SphereData_Spectral &R1_v = FUn_v;

		sweet::SphereData_Spectral R2_h = FAn_h + FBn_h;
		sweet::SphereData_Spectral R2_u = FAn_u + FBn_u;
		sweet::SphereData_Spectral R2_v = FAn_v + FBn_v;

		sweet::SphereData_Spectral &R3_h = FCn_h;
		sweet::SphereData_Spectral &R3_u = FCn_u;
		sweet::SphereData_Spectral &R3_v = FCn_v;


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
		ts_ups0_rexi.runTimestep(
				R0_h, R0_u, R0_v,
				dt,		i_simulation_timestamp
			);

		ts_ups1_rexi.runTimestep(
				R1_h, R1_u, R1_v,
				dt,		i_simulation_timestamp
			);

		ts_ups2_rexi.runTimestep(
				R2_h, R2_u, R2_v,
				dt,		i_simulation_timestamp
			);

		ts_ups3_rexi.runTimestep(
				R3_h, R3_u, R3_v,
				dt,		i_simulation_timestamp
			);

		io_phi_pert = R0_h + dt*(R1_h + 2.0*R2_h + R3_h);
		io_vrt = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		io_div = R0_v + dt*(R1_v + 2.0*R2_v + R3_v);
	}
	else
	{
		SWEETError("TODO: This order is not implemented, yet!");
	}
}




PDESWESphereTS_lg_exp_lc_n_etdrk::PDESWESphereTS_lg_exp_lc_n_etdrk()
{
}



PDESWESphereTS_lg_exp_lc_n_etdrk::~PDESWESphereTS_lg_exp_lc_n_etdrk()
{
}

