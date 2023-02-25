/*
 *  Created on: 29 Mar 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/pde_advectionPlane/time/PDEAdvPlaneTS_na_sl.hpp>


PDEAdvPlaneTS_na_sl::PDEAdvPlaneTS_na_sl(
		sweet::ShackDictionary &io_shackDict,
		sweet::PlaneOperators &i_op
)	:
		op(i_op)
{
	shackTimestepControl = io_shackDict.getAutoRegistration<sweet::ShackTimestepControl>();
	shackPlaneDataOps = io_shackDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	shackTimeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();
	shackBenchmark = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();

	ERROR_CHECK_WITH_RETURN(io_shackDict);

	_setup(shackTimeDisc->timestepping_order);
}



void PDEAdvPlaneTS_na_sl::run_timestep(
		sweet::PlaneData_Spectral &io_phi,		///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	assert(i_dt > 0);

	if (shackBenchmark->getExternalForcesCallback != nullptr)
	{
		shackBenchmark->getExternalForcesCallback(
				1,
				i_simulation_timestamp,
				&io_u,
				shackBenchmark
			);
		shackBenchmark->getExternalForcesCallback(
				2,
				i_simulation_timestamp,
				&io_v,
				shackBenchmark
			);
	}

	if (i_simulation_timestamp == 0)
	{
		prog_u_prev = io_u;
		prog_v_prev = io_v;
	}


	// OUTPUT: position of departure points at t
	sweet::ScalarDataArray posx_d(io_phi.planeDataConfig->physical_array_data_number_of_elements);
	sweet::ScalarDataArray posy_d(io_phi.planeDataConfig->physical_array_data_number_of_elements);

	semiLagrangian.semi_lag_departure_points_settls(
			prog_u_prev.toPhys(), prog_v_prev.toPhys(),
			io_u.toPhys(), io_v.toPhys(),
			posx_a, posy_a,
			i_dt,
			posx_d, posy_d,
			shackPlaneDataOps->plane_domain_size,
			nullptr,
			timestepping_order,

			shackTimeDisc->semi_lagrangian_max_iterations,
			shackTimeDisc->semi_lagrangian_convergence_threshold
	);

	prog_u_prev = io_u;
	prog_v_prev = io_v;

	sweet::PlaneData_Spectral new_prog_phi(io_phi.planeDataConfig);

	if (timestepping_order == 1)
	{
		sampler2D.bilinear_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else if (timestepping_order == 2)
	{
		sampler2D.bicubic_scalar(
				io_phi,
				posx_d,
				posy_d,
				new_prog_phi
		);
	}
	else
	{
		SWEETError("Timestepping order not available");
	}

	io_phi = new_prog_phi;
}



/*
 * Setup
 */
void PDEAdvPlaneTS_na_sl::_setup(
		int i_order	///< order of RK time stepping method
)
{
	prog_u_prev.setup(op.planeDataConfig);
	prog_v_prev.setup(op.planeDataConfig);

	timestepping_order = i_order;

	const sweet::PlaneDataConfig *planeDataConfig = op.planeDataConfig;

	posx_a.setup(planeDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(planeDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position
	posx_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % planeDataConfig->physical_res[0];

			io_data = (double)i*(double)shackPlaneDataOps->plane_domain_size[0]/(double)planeDataConfig->physical_res[0];

			assert(io_data >= 0.0);
			assert(io_data <= shackPlaneDataOps->plane_domain_size[0]);
		}
	);
	posy_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % planeDataConfig->physical_data_size[0];
			int j = idx / (double)planeDataConfig->physical_res[0];

			io_data = (double)j*(double)shackPlaneDataOps->plane_domain_size[1]/(double)planeDataConfig->physical_res[1];

			assert(io_data >= 0.0);
			assert(io_data <= shackPlaneDataOps->plane_domain_size[1]);
		}
	);

	// TODO: Use semiLagrangian.sampler2D
	sampler2D.setup(shackPlaneDataOps->plane_domain_size, planeDataConfig);

	//PXT- This just calls sampler2D.setup, so any reason for having it?
	semiLagrangian.setup(shackPlaneDataOps->plane_domain_size, planeDataConfig);
}


PDEAdvPlaneTS_na_sl::~PDEAdvPlaneTS_na_sl()
{
}
