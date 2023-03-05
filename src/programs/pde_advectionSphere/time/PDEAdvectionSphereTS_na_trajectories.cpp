/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDEAdvectionSphereTS_na_trajectories.hpp"
#include "../PDEAdvectionSphereBenchmarksCombined.hpp"
#include "../benchmarks/PDEAdvectionSphereBenchmark_nair_lauritzen_sl.hpp"



bool PDEAdvectionSphereTS_na_trajectories::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_trajectories";
}

std::string PDEAdvectionSphereTS_na_trajectories::string_id()
{
	return "na_trajectories";
}


std::string PDEAdvectionSphereTS_na_trajectories::get_help()
{
	std::ostringstream stream;
	stream << " + SphereAdvection_TS_na_trajectories:" << std::endl;
	stream << "    * 'na_trajectories'" << std::endl;

	return stream.str();
}



void PDEAdvectionSphereTS_na_trajectories::run_timestep(
		sweet::SphereData_Spectral &io_U_phi,		///< prognostic variables
		sweet::SphereData_Physical &io_U_u,
		sweet::SphereData_Physical &io_U_v,

		double i_fixed_dt,
		double i_simulation_timestamp,

		// for varying velocity fields
		const PDEAdvectionSphereBenchmarksCombined *i_sphereBenchmarks
)
{
	const sweet::SphereDataConfig *sphereDataConfig = io_U_phi.sphereDataConfig;
	std::size_t N = sphereDataConfig->physical_array_data_number_of_elements;


	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (shackPDEAdvBenchmark->getVelocities != nullptr)
	{
		shackPDEAdvBenchmark->getVelocities(io_U_u, io_U_v, i_simulation_timestamp, shackPDEAdvBenchmark->getVelocitiesUserData);
	}


	sweet::ScalarDataArray pos_lon_D(N), pos_lat_D(N);

	if (shackPDEAdvBenchmark->callback_slComputeDeparture3rdOrder != nullptr)
	{
		// compute 2nd order accurate departure points
		shackPDEAdvBenchmark->callback_slComputeDeparture3rdOrder(
				shackPDEAdvBenchmark->slComputeDeparture3rdOrderUserData,
				semiLagrangian.pos_lon_A,
				semiLagrangian.pos_lat_A,
				pos_lon_D,
				pos_lat_D,
				i_fixed_dt,
				i_simulation_timestamp+i_fixed_dt	// arrival time: current time + dt
			);
	}


	// sample phi at departure points
	sweet::SphereData_Physical U_phi_phys_D =
		sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi.toPhys(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			shackSemiLagrangian->semi_lagrangian_interpolation_limiter
		);

	io_U_phi = U_phi_phys_D;


	// sample phi at departure points

	U_phi_phys_D =
	sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi.getSphereDataPhysical(),
			pos_lon_D, pos_lat_D,
			false,	// velocity sampling
			false,
			shackSemiLagrangian->semi_lagrangian_interpolation_limiter
		);


}



bool PDEAdvectionSphereTS_na_trajectories::setup(
	sweet::SphereOperators *io_ops
)
{
	PDEAdvectionSphereTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDEAdvTimeDisc->timestepping_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		error.set("Only 1st and 2nd order for SL integration supported");

	semiLagrangian.setup(io_ops->sphereDataConfig, shackSemiLagrangian, timestepping_order);

	sphereSampler.setup(io_ops->sphereDataConfig);
	return true;
}


PDEAdvectionSphereTS_na_trajectories::~PDEAdvectionSphereTS_na_trajectories()
{
}

