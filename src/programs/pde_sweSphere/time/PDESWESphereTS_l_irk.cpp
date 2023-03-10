/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_irk.hpp"

#include <complex>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/TimeStepSizeChanged.hpp>
#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"


bool PDESWESphereTS_l_irk::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (
			i_timestepping_method == "l_irk"	||
			i_timestepping_method == "lg_irk"	||
			false
	)
	{
		timestepping_method = i_timestepping_method;
		return true;
	}

	return false;
}


std::string PDESWESphereTS_l_irk::getIDString()
{
	return timestepping_method;
}


bool PDESWESphereTS_l_irk::setup_auto(
		sweet::SphereOperators *io_ops
)
{
	if (shackPDESWESphere->sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	if (timestepping_method == "l_irk")
	{
		return setup(
				io_ops,
				timestepping_order,
				shackTimestepControl->current_timestep_size,
				0.5,
				false
			);
	}
	else if (timestepping_method == "lg_irk")
	{
		return setup(
				io_ops,
				timestepping_order,
				shackTimestepControl->current_timestep_size,
				0.5,
				true
			);
	}
	else
	{
		SWEETError("ERROR");
	}
	return false;
}

bool PDESWESphereTS_l_irk::setup(
		sweet::SphereOperators *io_ops,
		int i_timestep_order,
		double i_timestep_size
)
{
	return setup(
			io_ops,
			i_timestep_order,
			i_timestep_size,
			shackPDESWETimeDisc->timestepping_crank_nicolson_filter,
			false
		);
}


bool PDESWESphereTS_l_irk::setup(
		sweet::SphereOperators *io_ops,
		int i_timestepping_order,
		double i_timestep_size,
		double i_crank_nicolson_damping_factor,
		bool i_no_coriolis
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;
	timestep_size = i_timestep_size;
	crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;

	use_f_sphere = shackPDESWESphere->sphere_use_fsphere;
	no_coriolis = i_no_coriolis;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETError("Only 1st and 2nd order supported!");

	clear();

	if (no_coriolis)
	{
		lg_erk = new PDESWESphereTS_lg_erk;
		lg_erk->setup(ops, 1);
	}
	else
	{
		if (use_f_sphere)
		{
			f0 = shackPDESWESphere->sphere_fsphere_f0;
			two_coriolis = 0.0;
		}
		else
		{
			f0 = 0.0;
			two_coriolis = 2.0*shackPDESWESphere->sphere_rotating_coriolis_omega;
		}

		l_erk = new PDESWESphereTS_l_erk;
		l_erk->setup(ops, 1);
	}


	sphere_radius = shackSphereDataOps->sphere_radius;

	update_coefficients(i_timestep_size);
	return true;
}


void PDESWESphereTS_l_irk::runTimestep(
		sweet::SphereData_Spectral &io_phi,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (TimeStepSizeChanged::is_changed(timestep_size, i_fixed_dt, true))
		update_coefficients(i_fixed_dt);


	/*
	 * Crank-Nicolson method:
	 *
	 * (U(t+1) - q dt F(U(t+1))) = (U(t) + q dt F(U(t)))
	 *
	 * with q the CN damping facor with no damping for q=0.5
	 */
	if (timestepping_order == 2)
	{
		/*
		 * Explicit Euler
		 */
		if (no_coriolis)
			lg_erk->runTimestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
		else
			l_erk->runTimestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
	}


	double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	if (no_coriolis)
	{
		sweet::SphereData_Spectral rhs = io_div + ops->implicit_L(io_phi, dt_implicit);
		sweet::SphereData_Spectral div1 = ops->implicit_helmholtz(rhs, -gh0*dt_implicit*dt_implicit, sphere_radius);
		sweet::SphereData_Spectral phi1 = io_phi - dt_implicit*gh0*div1;
		sweet::SphereData_Spectral vrt1 = io_vrt;

		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
	else
	{
		double dt_two_omega = dt_implicit*2.0*shackPDESWESphere->sphere_rotating_coriolis_omega;

		sweet::SphereData_Spectral rhs = io_div + ops->implicit_FJinv(io_vrt, dt_two_omega) + ops->implicit_L(io_phi, dt_implicit);
		sweet::SphereData_Spectral div1 = sphSolverDiv.solve(rhs);

		sweet::SphereData_Spectral phi1 = io_phi - dt_implicit*gh0*div1;
		sweet::SphereData_Spectral vrt1 = ops->implicit_Jinv(io_vrt - ops->implicit_F(div1, dt_two_omega), dt_two_omega);


		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
}

void PDESWESphereTS_l_irk::solveImplicit(
		sweet::SphereData_Spectral &io_phi,	///< rhs variables
		sweet::SphereData_Spectral &io_vrt,	///< rhs variables
		sweet::SphereData_Spectral &io_div,	///< rhs variables

		double dt
)
{
	update_coefficients(dt);
	double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	if (no_coriolis)
	{
		SWEETError("Not supported yet");
	}
	else
	{
		double dt_two_omega = dt*2.0*shackPDESWESphere->sphere_rotating_coriolis_omega;

		// Implicit update using explicit evaluation for implicit_FJinv and implicit_L (with or without NL update)
		sweet::SphereData_Spectral rhs_div = io_div + ops->implicit_FJinv(io_vrt, dt_two_omega) + ops->implicit_L(io_phi, dt);
		sweet::SphereData_Spectral div1 = sphSolverDiv.solve(rhs_div);

		// Update for phi using implicit update for div
		sweet::SphereData_Spectral phi1 = io_phi - dt*gh0*div1;

		// Decoupled implicit update, using the implicit update for div 
		sweet::SphereData_Spectral vrt1 = ops->implicit_Jinv(io_vrt - ops->implicit_F(div1, dt_two_omega), dt_two_omega);

		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
}



void PDESWESphereTS_l_irk::update_coefficients(double i_timestep_size)
{
	timestep_size = i_timestep_size;

	double 	gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;

	if (timestepping_order == 1)
	{
		dt_explicit = -666.0;  // Yeah !!!
		dt_implicit = timestep_size;
	}
	else
	{
		dt_explicit = timestep_size*(1.0-crank_nicolson_damping_factor);
		dt_implicit = timestep_size*crank_nicolson_damping_factor;
	}

	if (!no_coriolis)
	{
		if (!use_f_sphere)
		{
			double dt_two_omega = dt_implicit*two_coriolis;

			sphSolverDiv.setup(sphereDataConfig, 4);
			sphSolverDiv.solver_component_implicit_J(dt_two_omega);
			sphSolverDiv.solver_component_implicit_FJinvF(dt_two_omega);
			sphSolverDiv.solver_component_implicit_L(gh0*dt_implicit, dt_implicit, sphere_radius);
		}
	}
}



PDESWESphereTS_l_irk::PDESWESphereTS_l_irk()
{
}


void PDESWESphereTS_l_irk::clear()
{
	if (lg_erk != nullptr)
	{
		delete lg_erk;
		lg_erk = nullptr;
	}

	if (l_erk != nullptr)
	{
		delete l_erk;
		l_erk = nullptr;
	}
}




PDESWESphereTS_l_irk::~PDESWESphereTS_l_irk()
{
	clear();
}
