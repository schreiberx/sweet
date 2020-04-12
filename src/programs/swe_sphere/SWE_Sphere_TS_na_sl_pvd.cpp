/*
 * SWE_Sphere_TS_na_sl_pvd.cpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2019-10-24: Partly based on plane version
 */

#include "SWE_Sphere_TS_na_sl_pvd.hpp"
#include <sweet/sphere/SphereData_DebugContainer.hpp>





void SWE_Sphere_TS_na_sl_pvd::run_timestep_pert_1st_order(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	FatalError("TODO: Implement na_sl with 1st order");

	const SphereData_Config *sphereDataConfig = io_phi_pert.sphereDataConfig;

	double dt_radius = i_dt/simVars.sim.sphere_radius;

	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_na_sl_pvd: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		U_phi_pert_prev = io_phi_pert;
		U_vort_prev = io_vort;
		U_div_prev = io_div;
	}

	// Output variables
	SphereData_Spectral phi(sphereDataConfig);
	SphereData_Spectral vort(sphereDataConfig);
	SphereData_Spectral div(sphereDataConfig);

	/*
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See Hortal's paper for equation.
	 */

	SphereData_Physical u_lon_prev(sphereDataConfig);
	SphereData_Physical v_lat_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vort_prev, U_div_prev, u_lon_prev, v_lat_prev, false);

	SphereData_Physical u_lon(sphereDataConfig);
	SphereData_Physical v_lat(sphereDataConfig);
	op.vortdiv_to_uv(io_vort, io_div, u_lon, v_lat, false);

	// Departure points and arrival points
	ScalarDataArray pos_lon_d = pos_lon_a;
	ScalarDataArray pos_lat_d = pos_lat_a;

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			dt_radius*u_lon_prev, dt_radius*v_lat_prev,
			dt_radius*u_lon, dt_radius*v_lat,

			pos_lon_a, pos_lat_a,
			pos_lon_d, pos_lat_d,		// OUTPUT

			simVars.disc.timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry, simVars.disc.semi_lagrangian_interpolation_limiter
		);


	/*
	 * Backup fields of previous time step
	 */
	U_phi_pert_prev = io_phi_pert;
	U_vort_prev = io_vort;
	U_div_prev = io_div;

	/*
	 * Phi SL
	 */
	SphereData_Physical phi_D_phys(sphereDataConfig);

	sphereSampler.bicubic_scalar(
			io_phi_pert.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			phi_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	io_phi_pert = SphereData_Spectral(phi_D_phys);

#if 1
	/*
	 * Vort/Div SL
	 */

	SphereData_Physical vort_D_phys(sphereDataConfig);
	SphereData_Physical div_D_phys(sphereDataConfig);

	sphereSampler.bicubic_scalar(
			io_vort.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			vort_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			io_div.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			div_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	io_vort = SphereData_Spectral(vort_D_phys);
	io_div = SphereData_Spectral(div_D_phys);

#else

	/*
	 * U/V SL
	 *
	 * WARNING: We need to use U/V space, since just treating
	 * vorticity and divergence with SL methods doesn't include all terms which
	 * would be included by treating U and V.
	 */

	SphereData_Physical u_D_phys(sphereDataConfig);
	SphereData_Physical v_D_phys(sphereDataConfig);

	sphereSampler.bicubic_scalar(
			u_lon,
			pos_lon_d,
			pos_lat_d,
			u_D_phys,
			true,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			v_lat,
			pos_lon_d,
			pos_lat_d,
			v_D_phys,
			true,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);


	op.uv_to_vortdiv(u_D_phys, v_D_phys, io_vort, io_div);

#endif


}



void SWE_Sphere_TS_na_sl_pvd::run_timestep_pert_2nd_order(
		SphereData_Spectral &io_U_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_U_vort,	///< prognostic variables
		SphereData_Spectral &io_U_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp)
{

	const SphereData_Config *sphereDataConfig = io_U_phi_pert.sphereDataConfig;

	double dt_radius = i_dt/simVars.sim.sphere_radius;

	if (i_dt <= 0)
		FatalError("SWE_Sphere_TS_na_sl_pvd: Only constant time step size allowed (Please set --dt)");

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		U_phi_pert_prev = io_U_phi_pert;
		U_vort_prev = io_U_vort;
		U_div_prev = io_U_div;
	}

	// Output variables
	SphereData_Spectral phi(sphereDataConfig);
	SphereData_Spectral vort(sphereDataConfig);
	SphereData_Spectral div(sphereDataConfig);

	/*
	 * Step 1) SL
	 * Compute Lagrangian trajectories based on SETTLS.
	 * This depends on V(t-\Delta t) and V(t).
	 *
	 * See Hortal's paper for equation.
	 */

	SphereData_Physical u_lon_prev(sphereDataConfig);
	SphereData_Physical v_lat_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vort_prev, U_div_prev, u_lon_prev, v_lat_prev, false);

	SphereData_Physical u_lon(sphereDataConfig);
	SphereData_Physical v_lat(sphereDataConfig);
	op.vortdiv_to_uv(io_U_vort, io_U_div, u_lon, v_lat, false);

	// Departure points and arrival points
	ScalarDataArray pos_lon_d = pos_lon_a;
	ScalarDataArray pos_lat_d = pos_lat_a;

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			dt_radius*u_lon_prev, dt_radius*v_lat_prev,
			dt_radius*u_lon, dt_radius*v_lat,

			pos_lon_a, pos_lat_a,
			pos_lon_d, pos_lat_d,		// OUTPUT

			simVars.disc.timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry, simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical phi_pert_D_phys(sphereDataConfig);
	SphereData_Physical vort_D_phys(sphereDataConfig);
	SphereData_Physical div_D_phys(sphereDataConfig);

	sphereSampler.bicubic_scalar(
			io_U_phi_pert.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			phi_pert_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			io_U_vort.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			vort_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	sphereSampler.bicubic_scalar(
			io_U_div.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			div_D_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	U_phi_pert_prev = io_U_phi_pert;
	U_vort_prev = io_U_vort;
	U_div_prev = io_U_div;

	io_U_phi_pert = SphereData_Spectral(phi_pert_D_phys);
	io_U_vort = SphereData_Spectral(vort_D_phys);
	io_U_div = SphereData_Spectral(div_D_phys);
}



void SWE_Sphere_TS_na_sl_pvd::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp)
{
	if (timestepping_order == 1)
	{
		run_timestep_pert_1st_order(io_phi_pert, io_vort, io_div, i_dt, i_simulation_timestamp);
	}
	else if (timestepping_order == 2)
	{
		run_timestep_pert_2nd_order(io_phi_pert, io_vort, io_div, i_dt, i_simulation_timestamp);
	}
	else
	{
		FatalError("Only orders 1/2 supported (ERRID 098nd89eje)");
	}
}



void SWE_Sphere_TS_na_sl_pvd::setup(
		int i_timestepping_order
)
{
	std::cout << " + SWE_Sphere_TS_na_sl_pvd.setup() called" << std::endl;

	timestepping_order = i_timestepping_order;

	if (timestepping_order != 1 && timestepping_order != 2)
		FatalError("Invalid time stepping order, must be 1 or 2");

	// Setup sampler for future interpolations
	sphereSampler.setup(op.sphereDataConfig);

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig, simVars);

	pos_lon_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);
	pos_lat_a.setup(op.sphereDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position

	pos_lon_a.update_lambda_array_indices([&](int idx, double &io_data)
	{
		int i = idx % op.sphereDataConfig->physical_num_lon;
		//int j = idx / sphereDataConfig->physical_data_size[0];

			io_data = 2.0*M_PI*(double)i/(double)op.sphereDataConfig->physical_num_lon;
			assert(io_data >= 0);
			assert(io_data < 2.0*M_PI);
		}
		);

	pos_lat_a.update_lambda_array_indices([&](int idx, double &io_data)
	{
		//int i = idx % sphereDataConfig->physical_data_size[0];
			int j = idx / op.sphereDataConfig->physical_num_lon;

			io_data = op.sphereDataConfig->lat[j];

			assert(io_data >= -M_PI*0.5);
			assert(io_data <= M_PI*0.5);
		}
	);

}

SWE_Sphere_TS_na_sl_pvd::SWE_Sphere_TS_na_sl_pvd(SimulationVariables &i_simVars, SphereOperators_SphereData &i_op) :
		simVars(i_simVars), op(i_op),

		U_phi_pert_prev(i_op.sphereDataConfig),
		U_vort_prev(i_op.sphereDataConfig),
		U_div_prev(i_op.sphereDataConfig),
		pos_lon_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		pos_lat_a(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posx_d(i_op.sphereDataConfig->physical_array_data_number_of_elements),
		posy_d(i_op.sphereDataConfig->physical_array_data_number_of_elements)
{
}

SWE_Sphere_TS_na_sl_pvd::~SWE_Sphere_TS_na_sl_pvd()
{
}

