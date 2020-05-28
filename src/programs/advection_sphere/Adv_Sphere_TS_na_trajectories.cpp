/*
 * Adv_Sphere_TS_na_sl.cpp
 *
 *  Created on: 29 Mar 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include <benchmarks_sphere/SWESphereBenchmark_nair_lauritzen_sl.hpp>
#include "Adv_Sphere_TS_na_trajectories.hpp"
#include "Adv_Sphere_TS_na_erk.hpp"



void Adv_Sphere_TS_na_trajectories::run_timestep(
		SphereData_Spectral &io_U_phi,		///< prognostic variables
		SphereData_Spectral &io_U_vrt,		///< prognostic variables
		SphereData_Spectral &io_U_div,		///< prognostic variables

		double i_fixed_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const SWESphereBenchmarks *i_sphereBenchmarks,
		SphereData_Physical &io_U_phi_phys
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;
	std::size_t N = sphereDataConfig->physical_array_data_number_of_elements;

	static SphereData_Spectral test;
	static SphereData_Physical test_phys;

	if (i_simulation_timestamp == 0)
	{
		test.setup(sphereDataConfig);
		test = io_U_phi;

		test_phys.setup(sphereDataConfig);
		test_phys = io_U_phi_phys;
	}

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (i_sphereBenchmarks)
	{
		SphereData_Spectral tmp(io_U_phi.sphereDataConfig);
		i_sphereBenchmarks->master->get_reference_state(tmp, io_U_vrt, io_U_div, i_simulation_timestamp);
		i_sphereBenchmarks->master->get_reference_state(tmp, U_vrt_prev, U_div_prev, i_simulation_timestamp - i_fixed_dt);
	}


#if 0
	double K = 64;
	double t = i_simulation_timestamp;
	double dt = i_fixed_dt/K;

	ScalarDataArray pos_lon_D_1(N), pos_lat_D_1(N);
	// compute 2nd order accurate departure points
	i_sphereBenchmarks->compute_departure_3rd_order(
			semiLagrangian.pos_lon_A,
			semiLagrangian.pos_lat_A,
			pos_lon_D_1,
			pos_lat_D_1,
			dt,
			t+dt	// arrival time: current time + dt
		);


	ScalarDataArray pos_lon_D_2(N), pos_lat_D_2(N);
	// compute 2nd order accurate departure points
	i_sphereBenchmarks->compute_departure_3rd_order(
			semiLagrangian.pos_lon_A,
			semiLagrangian.pos_lat_A,
			pos_lon_D_2,
			pos_lat_D_2,
			dt,
			t+dt	// arrival time: current time + dt
		);

	pos_lon_D = pos_lon_D_tmp;
	pos_lat_D = pos_lat_D_tmp;
	t += dt;

#elif 0

	ScalarDataArray pos_lon_D = semiLagrangian.pos_lon_A;
	ScalarDataArray pos_lat_D = semiLagrangian.pos_lat_A;

	double K = 64;
	double t = i_simulation_timestamp;
	double dt = i_fixed_dt/K;

	ScalarDataArray pos_lon_D_tmp(N), pos_lat_D_tmp(N);
	for (int i = 0; i < (int)K; i++)
	{
		// compute 2nd order accurate departure points
		i_sphereBenchmarks->compute_departure_3rd_order(
				pos_lon_D,
				pos_lat_D,
				pos_lon_D_tmp,
				pos_lat_D_tmp,
				dt,
				t+dt	// arrival time: current time + dt
			);

		pos_lon_D = pos_lon_D_tmp;
		pos_lat_D = pos_lat_D_tmp;
		t += dt;
	}

#else

	ScalarDataArray pos_lon_D(N), pos_lat_D(N);

	// compute 2nd order accurate departure points
	i_sphereBenchmarks->master->sl_compute_departure_3rd_order(
			semiLagrangian.pos_lon_A,
			semiLagrangian.pos_lat_A,
			pos_lon_D,
			pos_lat_D,
			i_fixed_dt,
			i_simulation_timestamp+i_fixed_dt	// arrival time: current time + dt
		);

#endif

	// sample phi at departure points
	SphereData_Physical U_phi_phys_D =
		sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi.toPhys(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	io_U_phi = U_phi_phys_D;


	// sample phi at departure points

	U_phi_phys_D =
	sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi_phys,
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	io_U_phi_phys = U_phi_phys_D;

#if 0
	Adv_Sphere_TS_na_erk erk(simVars, op);
	erk.setup(2);

	#if 1
		erk.run_timestep(
				test, io_U_vrt, io_U_div,
				i_fixed_dt,
				i_simulation_timestamp,
				i_sphereBenchmarks,
				test_phys
			);
	#endif


	SphereData_DebugContainer::append((io_U_phi - test), "erk_diff");
	SphereData_DebugContainer::append((io_U_phi_phys - test_phys), "erk_diff phys");
#endif

}



/*
 * Setup
 */
void Adv_Sphere_TS_na_trajectories::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		SWEETError("Only 1st and 2nd order for SL integration supported");
}


Adv_Sphere_TS_na_trajectories::Adv_Sphere_TS_na_trajectories(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)
{
	setup(simVars.disc.timestepping_order);

	semiLagrangian.setup(op.sphereDataConfig);
}



Adv_Sphere_TS_na_trajectories::~Adv_Sphere_TS_na_trajectories()
{
}

