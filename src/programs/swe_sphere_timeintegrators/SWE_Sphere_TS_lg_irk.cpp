/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */


#include "SWE_Sphere_TS_lg_irk.hpp"

#include <complex>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"



bool SWE_Sphere_TS_lg_irk::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
	if (i_timestepping_method == "lg_irk_DEPRECATED")
		return true;

	return false;
}



std::string SWE_Sphere_TS_lg_irk::string_id()
{
	return "lg_irk_DEPRECATED";
}



void SWE_Sphere_TS_lg_irk::setup_auto()
{
	if (simVars.sim.sphere_use_fsphere)
		SWEETError("TODO: Not yet supported");

	setup(
			timestepping_order,
			simVars.timecontrol.current_timestep_size,
			simVars.disc.timestepping_crank_nicolson_filter
		);
}



void SWE_Sphere_TS_lg_irk::run_timestep(
		SphereData_Spectral &io_phi_pert,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

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
		lg_erk->run_timestep(io_phi_pert, io_vrt, io_div, i_fixed_dt*(1.0-crank_nicolson_damping_factor), i_simulation_timestamp);
	}

	SphereData_Spectral phi0 = io_phi_pert;
	SphereData_Spectral vrt0 = io_vrt;
	SphereData_Spectral div0 = io_div;

	SphereData_Spectral phi(sphereDataConfig);
	SphereData_Spectral vort(sphereDataConfig);
	SphereData_Spectral div(sphereDataConfig);

	{
		SphereData_Spectral rhs = gh*div0 + alpha*phi0;
		phi = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_phi_pert = phi*beta;

		rhs = alpha*div0 + op.laplace(phi0);
		div = rhs.spectral_solve_helmholtz(alpha*alpha, -gh, r);
		io_div = div*beta;

		io_vrt = vrt0;
	}
}



SWE_Sphere_TS_lg_irk::SWE_Sphere_TS_lg_irk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
	simVars(i_simVars),
	op(i_op),
	sphereDataConfig(op.sphereDataConfig)
{
}


void SWE_Sphere_TS_lg_irk::update_coefficients()
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



void SWE_Sphere_TS_lg_irk::setup(
		int i_timestep_order,
		double i_timestep_size
)
{
	setup(i_timestep_order, i_timestep_size, simVars.disc.timestepping_crank_nicolson_filter);
}



void SWE_Sphere_TS_lg_irk::setup(
		int i_timestep_order,
		double i_timestep_size,
		double i_crank_nicolson_damping_factor
)
{
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
		lg_erk = new SWE_Sphere_TS_lg_erk(simVars, op);
		lg_erk->setup(1);
	}
	else
	{
		SWEETError("Only 1st and 2nd order IRK supported so far with this implementation! Use l_cn if you want to have 2nd order Crank-Nicolson!");
	}


	update_coefficients();

	r = simVars.sim.sphere_radius;
	inv_r = 1.0/r;

	gh = simVars.sim.gravitation*simVars.sim.h0;


}



SWE_Sphere_TS_lg_irk::~SWE_Sphere_TS_lg_irk()
{
	delete lg_erk;
	lg_erk = nullptr;
}
