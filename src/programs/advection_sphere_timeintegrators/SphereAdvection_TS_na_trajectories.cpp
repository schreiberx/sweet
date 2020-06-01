/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_trajectories.hpp"
//#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_erk.hpp"

#include "../advection_sphere_benchmarks/BenchmarksSphereAdvection.hpp"
#include "../advection_sphere_benchmarks/BenchmarksSphereAdvection_nair_lauritzen_sl.hpp"



bool SphereAdvection_TS_na_trajectories::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_trajectories";
}

std::string SphereAdvection_TS_na_trajectories::string_id()
{
	return "na_trajectories";
}


void SphereAdvection_TS_na_trajectories::setup_auto()
{
	setup(simVars.disc.timestepping_order);
}


std::string SphereAdvection_TS_na_trajectories::get_help()
{
	std::ostringstream stream;
	stream << " + SphereAdvection_TS_na_trajectories:" << std::endl;
	stream << "    * 'na_trajectories'" << std::endl;

	return stream.str();
}



void SphereAdvection_TS_na_trajectories::run_timestep(
		SphereData_Spectral &io_U_phi,		///< prognostic variables
		SphereData_Spectral &io_U_vrt,		///< prognostic variables
		SphereData_Spectral &io_U_div,		///< prognostic variables

		double i_fixed_dt,					///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const BenchmarksSphereAdvection *i_sphereBenchmarks
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;
	std::size_t N = sphereDataConfig->physical_array_data_number_of_elements;


	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (i_sphereBenchmarks)
	{
		SphereData_Spectral tmp(sphereDataConfig);
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
			io_U_phi.getSphereDataPhysical(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);


}



/*
 * Setup
 */
void SphereAdvection_TS_na_trajectories::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		SWEETError("Only 1st and 2nd order for SL integration supported");
}


SphereAdvection_TS_na_trajectories::SphereAdvection_TS_na_trajectories(
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



SphereAdvection_TS_na_trajectories::~SphereAdvection_TS_na_trajectories()
{
}

