/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#include "PDESWESphereTS_lg_irk.hpp"

#include <complex>
#include <sweet/core/SWEETError.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"



bool PDESWESphereTS_lg_irk::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	/*
	 * Supported directly in l_irk, not in this class anymore
	 */
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "lg_irk_DEPRECATED")
		return true;

	return false;
}



std::string PDESWESphereTS_lg_irk::getIDString()
{
	return "lg_irk_DEPRECATED";
}



bool PDESWESphereTS_lg_irk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	if (shackPDESWESphere->sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	return setup_main(
			io_ops,
			shackPDESWETimeDisc->timestepping_order,
			shackTimestepControl->current_timestep_size,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter
		);
}


bool PDESWESphereTS_lg_irk::setup(
		sweet::SphereOperators *io_ops,
		int i_timestep_order,
		double i_timestep_size
)
{
	return setup_main(
			io_ops,
			i_timestep_order,
			i_timestep_size,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter
		);
}



bool PDESWESphereTS_lg_irk::setup_main(
		sweet::SphereOperators *io_ops,
		int i_timestep_order,
		double i_timestep_size,
		double i_crank_nicolson_damping_factor
)
{
	ops = io_ops;

	timestepping_order = i_timestep_order;
	timestep_size = i_timestep_size;

	if (i_timestep_order == 1)
	{
		// set this to 1 to ignore it
		crank_nicolson_damping_factor = 1.0;
	}
	else if (i_timestep_order == 2)
	{
		crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;
		lg_erk.setup(ops, 1);
	}
	else
	{
		SWEETError("Only 1st and 2nd order IRK supported so far with this implementation! Use l_cn if you want to have 2nd order Crank-Nicolson!");
	}


	update_coefficients();

	r = shackSphereDataOps->sphere_radius;
	inv_r = 1.0/r;

	gh = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	return true;
}




void PDESWESphereTS_lg_irk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt,		
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	if (std::abs(timestep_size - i_fixed_dt)/std::max(timestep_size, i_fixed_dt) > 1e-10)
	{
		//std::cout << "Warning: Reducing time step size from " << i_fixed_dt << " to " << timestep_size << std::endl;
		timestep_size = i_fixed_dt;

		update_coefficients();
	}

	if (timestepping_order == 2)
	{
		/*
		 * Execute a half ERK time step first for 2nd order
		 */
		lg_erk.runTimestep(io_phi_pert, io_vrt, io_div, i_fixed_dt*(1.0-crank_nicolson_damping_factor), i_simulation_timestamp);
	}

	sweet::SphereData_Spectral phi0 = io_phi_pert;
	sweet::SphereData_Spectral vrt0 = io_vrt;
	sweet::SphereData_Spectral div0 = io_div;

	sweet::SphereData_Spectral phi(ops->sphereDataConfig);
	sweet::SphereData_Spectral vort(ops->sphereDataConfig);
	sweet::SphereData_Spectral div(ops->sphereDataConfig);

	{
		sweet::SphereData_Spectral rhs = gh*div0 + alpha*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_phi_pert = phi*beta;

		rhs = alpha*div0 + ops->laplace(phi0);
		div = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_div = div*beta;

		io_vrt = vrt0;
	}
}



PDESWESphereTS_lg_irk::PDESWESphereTS_lg_irk()
{
}


void PDESWESphereTS_lg_irk::update_coefficients()
{
	alpha = -1.0/timestep_size;
	beta = -1.0/timestep_size;

	{
		/*
		 * Crank-Nicolson method:
		 *
		 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
		 *
		 * with q the CN damping facor with no damping for q=0.5
		 */

		// scale dt by the damping factor to reuse solver structure

		alpha /= crank_nicolson_damping_factor;
		beta /= crank_nicolson_damping_factor;
	}
}



PDESWESphereTS_lg_irk::~PDESWESphereTS_lg_irk()
{
}
