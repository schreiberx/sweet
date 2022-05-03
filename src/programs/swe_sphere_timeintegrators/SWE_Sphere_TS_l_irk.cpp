/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_irk.hpp"

#include <complex>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/TimeStepSizeChanged.hpp>
#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"


bool SWE_Sphere_TS_l_irk::implements_timestepping_method(const std::string &i_timestepping_method
#if SWEET_PARAREAL
									,
									int &i_timestepping_order,
									int &i_timestepping_order2
#endif
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
#if SWEET_PARAREAL
	timestepping_order = i_timestepping_order;
	timestepping_order2 = i_timestepping_order2;
#endif
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


std::string SWE_Sphere_TS_l_irk::string_id()
{
	return timestepping_method;
}


void SWE_Sphere_TS_l_irk::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	if (timestepping_method == "l_irk")
	{
		setup(timestepping_order,
				simVars.timecontrol.current_timestep_size,
				0.5,
				false
			);
	}
	else if (timestepping_method == "lg_irk")
	{
		setup(timestepping_order,
				simVars.timecontrol.current_timestep_size,
				0.5,
				true
			);
	}
	else
	{
		SWEETError("ERROR");
	}
}


void SWE_Sphere_TS_l_irk::run_timestep(
		SphereData_Spectral &io_phi,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

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

	SphereData_Spectral o_phi_t(sphereDataConfig);
	SphereData_Spectral o_vrt_t(sphereDataConfig);
	SphereData_Spectral o_div_t(sphereDataConfig);


	if (timestepping_order == 2)
	{
		/*
		 * Explicit Euler
		 */
		if (no_coriolis)
			lg_erk->run_timestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
		else
			l_erk->run_timestep(io_phi, io_vrt, io_div, dt_explicit, i_simulation_timestamp);
	}


	double gh0 = simVars.sim.gravitation*simVars.sim.h0;

	if (no_coriolis)
	{
		SphereData_Spectral rhs = io_div + ops.implicit_L(io_phi, dt_implicit);
		SphereData_Spectral div1 = ops.implicit_helmholtz(rhs, -gh0*dt_implicit*dt_implicit, sphere_radius);
		SphereData_Spectral phi1 = io_phi - dt_implicit*gh0*div1;
		SphereData_Spectral vrt1 = io_vrt;

		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
	else
	{
		double dt_two_omega = dt_implicit*2.0*simVars.sim.sphere_rotating_coriolis_omega;

		SphereData_Spectral rhs = io_div + ops.implicit_FJinv(io_vrt, dt_two_omega) + ops.implicit_L(io_phi, dt_implicit);
		SphereData_Spectral div1 = sphSolverDiv.solve(rhs);

		SphereData_Spectral phi1 = io_phi - dt_implicit*gh0*div1;
		SphereData_Spectral vrt1 = ops.implicit_Jinv(io_vrt - ops.implicit_F(div1, dt_two_omega), dt_two_omega);


		io_phi = phi1;
		io_vrt = vrt1;
		io_div = div1;
	}
}



void SWE_Sphere_TS_l_irk::update_coefficients(double i_timestep_size)
{
	timestep_size = i_timestep_size;

	double 	gh0 = simVars.sim.gravitation*simVars.sim.h0;

	if (timestepping_order == 1)
	{
		dt_explicit = -666.0;
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



SWE_Sphere_TS_l_irk::SWE_Sphere_TS_l_irk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
	simVars(i_simVars),
	ops(i_op),
	sphereDataConfig(ops.sphereDataConfig)
{
}


void SWE_Sphere_TS_l_irk::setup(
		int i_timestepping_order,
		double i_timestep_size,
		double i_crank_nicolson_damping_factor,
		bool i_no_coriolis
)
{
	timestepping_order = i_timestepping_order;
	timestep_size = i_timestep_size;
	crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;

	use_f_sphere = simVars.sim.sphere_use_fsphere;
	no_coriolis = i_no_coriolis;

	if (timestepping_order != 1 && timestepping_order != 2)
		SWEETError("Only 1st and 2nd order supported!");

	free();

	if (no_coriolis)
	{
		lg_erk = new SWE_Sphere_TS_lg_erk(simVars, ops);
		lg_erk->setup(1);
	}
	else
	{
		if (use_f_sphere)
		{
			f0 = simVars.sim.sphere_fsphere_f0;
			two_coriolis = 0.0;
		}
		else
		{
			f0 = 0.0;
			two_coriolis = 2.0*simVars.sim.sphere_rotating_coriolis_omega;
		}

		l_erk = new SWE_Sphere_TS_l_erk(simVars, ops);
		l_erk->setup(1);
	}


	sphere_radius = simVars.sim.sphere_radius;

	update_coefficients(i_timestep_size);
}



void SWE_Sphere_TS_l_irk::free()
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



void SWE_Sphere_TS_l_irk::setup(
		int i_timestep_order,
		double i_timestep_size
)
{
	setup(i_timestep_order, i_timestep_size, simVars.disc.timestepping_crank_nicolson_filter, false);
}



SWE_Sphere_TS_l_irk::~SWE_Sphere_TS_l_irk()
{
	free();
}
