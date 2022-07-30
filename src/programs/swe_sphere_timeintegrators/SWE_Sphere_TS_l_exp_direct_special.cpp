/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_l_exp_direct_special.hpp"



void SWE_Sphere_TS_l_exp_direct_special::euler_timestep_store_update_lc(
		const SphereData_Spectral &i_phi_pert,
		const SphereData_Spectral &i_vrt,
		const SphereData_Spectral &i_div,
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div,
		double i_simulation_timestamp
)
{
	/*
	 * We have this special function for the lc updates since we need
	 * to initialize the output variables with 0 values and this makes
	 * it more comfortable.
	 */
	o_phi_pert.spectral_set_zero();
	o_vrt.spectral_set_zero();
	o_div.spectral_set_zero();

	timestepping_lc_erk.euler_timestep_update_lc(
			i_phi_pert, i_vrt, i_div,
			o_phi_pert, o_vrt, o_div,
			i_simulation_timestamp
		);
}



void SWE_Sphere_TS_l_exp_direct_special::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	const SphereData_Config *sphereDataConfig = io_phi_pert.sphereDataConfig;

	if (!use_coriolis)
	{
		timestepping_lg_exp_phi0.run_timestep(
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		return;
	}

	if (timestepping_order == 1)
	{
		/*
		 * U_{1} =	\psi_{0}( \Delta t L ) U_{0}
		 * 			+ \Delta t \psi_{1}(\Delta tL) N(U_{0}).
		 */

		SphereData_Spectral phi0_Un_phi_pert(sphereDataConfig), phi0_Un_vrt(sphereDataConfig), phi0_Un_div(sphereDataConfig);
		timestepping_lg_exp_phi0.run_timestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_phi_pert, phi0_Un_vrt, phi0_Un_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_phi_pert(sphereDataConfig), FUn_vrt(sphereDataConfig), FUn_div(sphereDataConfig);
		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_phi_pert, FUn_vrt, FUn_div,
				i_simulation_timestamp
			);

		SphereData_Spectral phi1_FUn_phi_pert(sphereDataConfig), phi1_FUn_vrt(sphereDataConfig), phi1_FUn_div(sphereDataConfig);
		timestepping_lg_exp_phi1.run_timestep(
				FUn_phi_pert, FUn_vrt, FUn_div,
				phi1_FUn_phi_pert, phi1_FUn_vrt, phi1_FUn_div,
				i_fixed_dt,
				i_simulation_timestamp
			);

		io_phi_pert = phi0_Un_phi_pert + i_fixed_dt*phi1_FUn_phi_pert;
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

		timestepping_lg_exp_phi0.run_timestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_h, phi0_Un_u, phi0_Un_v,
				i_fixed_dt,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_h(sphereDataConfig), FUn_u(sphereDataConfig), FUn_v(sphereDataConfig);

		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_h, FUn_u, FUn_v,
				i_simulation_timestamp
		);

		SphereData_Spectral phi1_FUn_h(sphereDataConfig), phi1_FUn_u(sphereDataConfig), phi1_FUn_v(sphereDataConfig);

		timestepping_lg_exp_phi1.run_timestep(
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

		euler_timestep_store_update_lc(
				A_h, A_u, A_v,
				FAn_h, FAn_u, FAn_v,
				i_simulation_timestamp
		);


		SphereData_Spectral phi2_X_h(sphereDataConfig), phi2_X_u(sphereDataConfig), phi2_X_v(sphereDataConfig);

		timestepping_lg_exp_phi2.run_timestep(
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
		return;
	}

	if (timestepping_order == 4)
	{
		double dt = i_fixed_dt;
		double dt_half = dt*0.5;

		/*
		 * Precompute commonly used terms
		 */
		SphereData_Spectral phi0_Un_half_phi_pert(sphereDataConfig), phi0_Un_half_vrt(sphereDataConfig), phi0_Un_half_div(sphereDataConfig);

		timestepping_lg_exp_phi0.run_timestep(
				io_phi_pert, io_vrt, io_div,
				phi0_Un_half_phi_pert, phi0_Un_half_vrt, phi0_Un_half_div,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral FUn_phi_pert(sphereDataConfig), FUn_vrt(sphereDataConfig), FUn_div(sphereDataConfig);
		euler_timestep_store_update_lc(
				io_phi_pert, io_vrt, io_div,
				FUn_phi_pert, FUn_vrt, FUn_div,
				i_simulation_timestamp
		);

		/*
		 * A_{n} = \phi_{0}\left(0.5\Delta tL\right)U_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)F(U_{n},t_{n})
		 */
		SphereData_Spectral phi1_phi_pert(sphereDataConfig), phi1_vrt(sphereDataConfig), phi1_div(sphereDataConfig);
		timestepping_lg_exp_phi1.run_timestep(
				FUn_phi_pert, FUn_vrt, FUn_div,
				phi1_phi_pert, phi1_vrt, phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral A_phi_pert = phi0_Un_half_phi_pert + dt_half*phi1_phi_pert;
		SphereData_Spectral A_vrt = phi0_Un_half_vrt + dt_half*phi1_vrt;
		SphereData_Spectral A_div = phi0_Un_half_div + dt_half*phi1_div;


		/*
		 * B_{n} = \phi_{0}\left(0.5\Delta tL\right)U_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)F(A_{n},t_{n}+0.5\Delta t)
		 */
		SphereData_Spectral FAn_phi_pert(sphereDataConfig), FAn_vrt(sphereDataConfig), FAn_div(sphereDataConfig);

		euler_timestep_store_update_lc(
				A_phi_pert, A_vrt, A_div,
				FAn_phi_pert, FAn_vrt, FAn_div,
				i_simulation_timestamp + dt_half
		);

		timestepping_lg_exp_phi1.run_timestep(
				FAn_phi_pert, FAn_vrt, FAn_div,
				phi1_phi_pert, phi1_vrt, phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral B_phi_pert = phi0_Un_half_phi_pert + dt_half*phi1_phi_pert;
		SphereData_Spectral B_vrt = phi0_Un_half_vrt + dt_half*phi1_vrt;
		SphereData_Spectral B_div = phi0_Un_half_div + dt_half*phi1_div;


		/*
		 * C_{n} = \phi_{0}\left(0.5\Delta tL\right)A_{n}+0.5\Delta t\phi_{1}\left(0.5\Delta tL\right)\left(2F(B_{n},t_{n}+0.5\Delta t)-F(U_{n},t_{n})\right)
		 */
		SphereData_Spectral FBn_phi_pert(sphereDataConfig), FBn_vrt(sphereDataConfig), FBn_div(sphereDataConfig);

		euler_timestep_store_update_lc(
				B_phi_pert, B_vrt, B_div,
				FBn_phi_pert, FBn_vrt, FBn_div,
				i_simulation_timestamp + dt_half
		);

		timestepping_lg_exp_phi1.run_timestep(
				2.0*FBn_phi_pert - FUn_phi_pert,
				2.0*FBn_vrt - FUn_vrt,
				2.0*FBn_div - FUn_div,
				phi1_phi_pert,	phi1_vrt,	phi1_div,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral phi0_An_phi_pert(sphereDataConfig), phi0_An_vrt(sphereDataConfig), phi0_An_div(sphereDataConfig);
		timestepping_lg_exp_phi0.run_timestep(
				A_phi_pert, A_vrt, A_div,
				phi0_An_phi_pert, phi0_An_vrt, phi0_An_div,
				dt_half,
				i_simulation_timestamp
			);

		SphereData_Spectral C_phi_pert = phi0_An_phi_pert + dt_half*phi1_phi_pert;
		SphereData_Spectral C_vrt = phi0_An_vrt + dt_half*phi1_vrt;
		SphereData_Spectral C_div = phi0_An_div + dt_half*phi1_div;



		/*
		 * R0 - R3
		 */
		SphereData_Spectral FCn_phi_pert(sphereDataConfig), FCn_vrt(sphereDataConfig), FCn_div(sphereDataConfig);

		euler_timestep_store_update_lc(
				C_phi_pert, C_vrt, C_div,
				FCn_phi_pert, FCn_vrt, FCn_div,
				i_simulation_timestamp + dt
		);

		SphereData_Spectral R0_phi_pert = io_phi_pert;
		SphereData_Spectral R0_vrt = io_vrt;
		SphereData_Spectral R0_div = io_div;

		SphereData_Spectral &R1_phi_pert = FUn_phi_pert;
		SphereData_Spectral &R1_vrt = FUn_vrt;
		SphereData_Spectral &R1_div = FUn_div;

		SphereData_Spectral R2_phi_pert = FAn_phi_pert + FBn_phi_pert;
		SphereData_Spectral R2_vrt = FAn_vrt + FBn_vrt;
		SphereData_Spectral R2_div = FAn_div + FBn_div;

		SphereData_Spectral &R3_phi_pert = FCn_phi_pert;
		SphereData_Spectral &R3_vrt = FCn_vrt;
		SphereData_Spectral &R3_div = FCn_div;


		/*
		 * U_{n+1} =
		 * 		\psi_{0}(\Delta tL) R_{0}
		 * 			+ \Delta t
		 * 			(
		 * 				  \upsilon_{1}(\Delta tL) R_{1} +
		 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
		 * 				  \upsilon_{3}(\Delta tL) R_{3}
		 * 			)
		 */
		timestepping_lg_exp_phi0.run_timestep(
				R0_phi_pert, R0_vrt, R0_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups1.run_timestep(
				R1_phi_pert, R1_vrt, R1_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups2.run_timestep(
				R2_phi_pert, R2_vrt, R2_div,
				dt,		i_simulation_timestamp
			);

		timestepping_lg_exp_ups3.run_timestep(
				R3_phi_pert, R3_vrt, R3_div,
				dt,		i_simulation_timestamp
			);

		io_phi_pert = R0_phi_pert + dt*(R1_phi_pert + 2.0*R2_phi_pert + R3_phi_pert);
		io_vrt = R0_vrt + dt*(R1_vrt + 2.0*R2_vrt + R3_vrt);
		io_div = R0_div + dt*(R1_div + 2.0*R2_div + R3_div);
		return;
	}

	SWEETError("This time stepping order is not yet supported!");
}



void SWE_Sphere_TS_l_exp_direct_special::setup(
	int i_order,
	bool i_use_coriolis,	///< Include Coriolis term
	const std::string i_function_name
)
{
	timestepping_order = i_order;
	use_coriolis = i_use_coriolis;

	if (use_coriolis)
		timestepping_method = "l_exp_special";
	else
		timestepping_method = "lg_exp_special";


	if (timestepping_order != 1 && timestepping_order != 2 && timestepping_order != 4)
	{
		std::ostringstream ss;
		ss << "Time stepping order " << i_order << " requested";
		SWEETError_nostop(ss.str());
		SWEETError("Only 1st and 2nd order time stepping order supported");
	}

	if (use_coriolis)
	{
		if (i_function_name != "phi0")
		{
			SWEETError("Only phi0 functions supported for the ETDnRK lg/lc scheme");
		}

		timestepping_lg_exp_phi0.setup("phi0");
		timestepping_lg_exp_phi1.setup("phi1");

		if (timestepping_order >= 2)
		{
			timestepping_lg_exp_phi2.setup("phi2");

			if (timestepping_order >= 4)
			{
				timestepping_lg_exp_ups1.setup("ups1");
				timestepping_lg_exp_ups2.setup("ups2");
				timestepping_lg_exp_ups3.setup("ups3");
			}
		}
	}
	else
	{
		timestepping_lg_exp_phi0.setup(i_function_name);
	}
}



void SWE_Sphere_TS_l_exp_direct_special::setup_auto()
{
	if (simVars.disc.timestepping_method == "lg_exp_special")
		setup(simVars.disc.timestepping_order, false, "phi0");
	else if (simVars.disc.timestepping_method == "l_exp_special")
		setup(simVars.disc.timestepping_order, true, "phi0");
}



SWE_Sphere_TS_l_exp_direct_special::SWE_Sphere_TS_l_exp_direct_special(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		timestepping_lg_exp_phi0(i_simVars, i_op),
		timestepping_lg_exp_phi1(i_simVars, i_op),
		timestepping_lg_exp_phi2(i_simVars, i_op),
		//timestepping_lg_exp_phi3(i_simVars, i_op),

		timestepping_lg_exp_ups1(i_simVars, i_op),
		timestepping_lg_exp_ups2(i_simVars, i_op),
		timestepping_lg_exp_ups3(i_simVars, i_op),

		timestepping_lc_erk(i_simVars, i_op)
{
}



SWE_Sphere_TS_l_exp_direct_special::~SWE_Sphere_TS_l_exp_direct_special()
{
}

