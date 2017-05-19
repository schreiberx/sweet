/*
 * Burgers_Plane.cpp
 *
 *  Created on: 15 May 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */
#include "Burgers_Plane.hpp"

#include <sweet/sweetmath.hpp>


#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>


#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>
#include <sweet/FatalError.hpp>


/**
 * Solve viscous Burgers' Equation with IMEX
 */
bool Burgers_Plane::run_timestep_imex(
		PlaneData &io_u,
		PlaneData &io_v,

		double& o_dt,			///< return time step size for the computed time step
		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters,
		PlaneDataConfig* planeDataConfig,

		double i_max_simulation_time
)
{
	if (i_parameters.misc.verbosity > 2)
		std::cout << "Burgers_Plane::run_timestep_imex()" << std::endl;

	PlaneData u=io_u;
	PlaneData v=io_v;

	// Modify timestep to final time if necessary
	double& t = o_dt;
	if (i_parameters.timecontrol.current_simulation_time+i_timestep_size < i_max_simulation_time)
		t = i_timestep_size;
	else
		t = i_max_simulation_time-i_parameters.timecontrol.current_simulation_time;

	// Initialize and set timestep dependent source for manufactured solution
	PlaneData f(planeDataConfig);
	PlaneData ff(planeDataConfig);
	BurgersValidationBenchmarks::set_source(i_parameters.timecontrol.current_simulation_time,i_parameters,i_parameters.disc.use_staggering,f);
	BurgersValidationBenchmarks::set_source(i_parameters.timecontrol.current_simulation_time+0.5*t,i_parameters,i_parameters.disc.use_staggering,ff);
	f.request_data_spectral();
	ff.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	rhs_u += - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + 0.5*t*f;
	rhs_v += - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

	if (i_parameters.disc.use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		lhs = ((-t)*i_parameters.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
        PlaneData u1 = rhs_u.spectral_div_element_wise(lhs);
        PlaneData v1 = rhs_v.spectral_div_element_wise(lhs);

        io_u = u + t*i_parameters.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
              - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
        io_v = v + t*i_parameters.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
              - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

	} else { //Jacobi
		FatalError("NOT available");
	}

	return true;
}


/**
 * Solve viscous Burgers' equation with explicit time stepping
 */
bool Burgers_Plane::run_timestep_explicit_ts(
		PlaneData &io_u,
		PlaneData &io_v,

		double& o_dt,			///< return time step size for the computed time step
		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters,
		PlaneDataConfig* planeDataConfig
)
{
	if (i_parameters.misc.verbosity > 2)
	std::cout << "run_timestep_imex()" << std::endl;

	PlaneData u=io_u;
	PlaneData v=io_v;

	// Modify timestep to final time if necessary
	double& t = o_dt;
	if (i_parameters.timecontrol.current_simulation_time+i_timestep_size < i_parameters.timecontrol.max_simulation_time)
		t = i_timestep_size;
	else
		t = i_parameters.timecontrol.max_simulation_time-i_parameters.timecontrol.current_simulation_time;

	// Initialize and set timestep dependent source for manufactured solution
	PlaneData f(planeDataConfig);
	PlaneData ff(planeDataConfig);

	BurgersValidationBenchmarks::set_source(i_parameters.timecontrol.current_simulation_time,i_parameters,i_parameters.disc.use_staggering,f);
	BurgersValidationBenchmarks::set_source(i_parameters.timecontrol.current_simulation_time+0.5*t,i_parameters,i_parameters.disc.use_staggering,ff);

	f.request_data_spectral();
	ff.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	if (i_parameters.disc.use_spectral_basis_diffs) //spectral
	{
		PlaneData u1 = u + t*i_parameters.sim.viscosity*(op.diff2_c_x(u)+op.diff2_c_y(u))
					   - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + f*t;
		PlaneData v1 = v + t*i_parameters.sim.viscosity*(op.diff2_c_x(v)+op.diff2_c_y(v))
					   - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

		io_u = u + t*i_parameters.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
			  - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
		io_v = v + t*i_parameters.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
			  - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));

	} else { //Jacobi
		FatalError("NOT available");
	}

	return true;
}


/**
 * Solve viscous Burgers' equation with implicit time stepping
 */
bool Burgers_Plane::run_timestep_implicit_ts(
		PlaneData &io_u,
		PlaneData &io_v,

		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_parameters
)
{
	FatalError("NOT available");
	return true;
}


/**
 * Solve  SWE with the novel Semi-Lagrangian Exponential Integrator
 *  SL-REXI
 *
 *  See documentation for details
 *
 */
bool Burgers_Plane::run_timestep_sl(
	PlaneData &io_u,  ///< Current and past fields
	PlaneData &io_v,
	PlaneData &io_u_prev,
	PlaneData &io_v_prev,

	ScalarDataArray &i_posx_a, //Arrival point positions in x and y (this is basically the grid)
	ScalarDataArray &i_posy_a,

	double& o_dt,			///< return time step size for the computed time step
	double i_timestep_size,	///< timestep size

	const SimulationVariables &i_simVars, ///< Parameters for simulation
	PlaneDataConfig* planeDataConfig,

	PlaneOperators &op,     ///< Operator class
	PlaneDataSampler &sampler2D, ///< Interpolation class
	SemiLagrangian &semiLagrangian,  ///< Semi-Lag class
	double* i_stag_displacement,
	double* i_stag_u,
	double* i_stag_v
)
{
	o_dt = i_timestep_size;
	//Padding for last time step
	if (i_simVars.timecontrol.current_simulation_time+o_dt > i_simVars.timecontrol.max_simulation_time)
		o_dt = i_simVars.timecontrol.max_simulation_time-i_simVars.timecontrol.current_simulation_time;


	//Departure points and arrival points
	ScalarDataArray posx_d(io_u.planeDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_u.planeDataConfig->physical_array_data_number_of_elements);

	//Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			io_u_prev, io_v_prev,
			io_u, io_v,
			i_posx_a, i_posy_a,
			o_dt,
			posx_d, posy_d,
			i_stag_displacement
			);

	// Save old velocities
	io_u_prev = io_u;
	io_v_prev = io_v;

	//Now interpolate to the the departure points
	//Departure points are set for physical space
	io_u = sampler2D.bicubic_scalar(
			io_u,
			posx_d,
			posy_d,
			i_stag_u[0],
			i_stag_u[1]
	);

	io_v = sampler2D.bicubic_scalar(
			io_v,
			posx_d,
			posy_d,
			i_stag_v[0],
			i_stag_v[1]
	);

	PlaneData u=io_u;
	PlaneData v=io_v;

	// Initialize and set timestep dependent source for manufactured solution
	PlaneData f(planeDataConfig);

	BurgersValidationBenchmarks::set_source(i_simVars.timecontrol.current_simulation_time+o_dt,i_simVars,i_simVars.disc.use_staggering,f);

	f.request_data_spectral();

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	rhs_u += o_dt*f;

	if (i_simVars.disc.use_spectral_basis_diffs) //spectral
	{
		PlaneData lhs = u;

		lhs = ((-o_dt)*i_simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		io_u = rhs_u.spectral_div_element_wise(lhs);
		io_v = rhs_v.spectral_div_element_wise(lhs);

	} else { //Jacobi
		FatalError("NOT available");
	}

	return true;
}



/**
 * This method computes the analytical solution based on the given initial values.
 */
void Burgers_Plane::run_timestep_direct_solution(
		PlaneData &io_h,
		PlaneData &io_u,
		PlaneData &io_v,

		double i_timestep_size,	///< timestep size

		PlaneOperators &op,
		const SimulationVariables &i_simVars
)
{
	FatalError("NOT available");
}
