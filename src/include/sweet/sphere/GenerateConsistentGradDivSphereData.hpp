/*
 * GenerateConsistentGradDivSphereData.hpp
 *
 *  Created on: 4 Nov 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_GENERATECONSISTENTGRADDIVSPHEREDATA_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_GENERATECONSISTENTGRADDIVSPHEREDATA_HPP_


#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>



/**
 * This class generates data which is consistent with SPH div/grad operators.
 *
 * One of the main problems by using SPH is, that the data has to fulfill
 * certain consistenty requirements.
 *
 * E.g. the GRAD operator should be computed on scalar fields only.
 * However, applying the DIV operator to a scalar field might result in
 * non-sense results.
 *
 * Also, just setting up data in an arbitrary way (e.g. using a Gaussian bump)
 * might lead to such inconsistent data.
 *
 * This class generates a consistent scalar and velocity field.
 *
 * This is done by initializing the (geopotential) scalar field and executing
 * a few time steps which generate a consistent velocity scheme.
 */
class GenerateConsistentGradDivSphereData
{
private:
	SimulationVariables &simVars;

	SphereDataConfig *sphereDataConfig;
	SphereOperators &op;

	std::string prefix_string;

public:
	SphereData prog_h, prog_u, prog_v;

	double center_lon, center_lat;


private:
	// Main routine for method to be used in case of finite differences
	void p_run_euler_timestep_update(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = i_fixed_dt;

		assert(simVars.sim.earth_radius > 0);

		if (!simVars.misc.use_nonlinear_equations)
		{
			if (!simVars.misc.sphere_use_robert_functions)
			{
				// linear equations
				o_h_t = -(op.div_lon(i_u)+op.div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += 2.0*simVars.sim.coriolis_omega*op.mu(i_v);
					o_v_t -= 2.0*simVars.sim.coriolis_omega*op.mu(i_u);
				}
			}
			else
			{
				// use Robert functions for velocity
				// linear equations
				o_h_t = -(op.robert_div_lon(i_u)+op.robert_div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t = -op.robert_grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.robert_grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += 2.0*simVars.sim.coriolis_omega*op.mu(i_v);
					o_v_t -= 2.0*simVars.sim.coriolis_omega*op.mu(i_u);
				}
			}
		}
		else
		{
			if (!simVars.misc.sphere_use_robert_functions)
			{
				/*
				 * Height
				 */
				// non-linear equations
				o_h_t = -(op.div_lon(i_h*i_u)+op.div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);

				/*
				 * Velocity
				 */
				// linear terms
				o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += 2.0*simVars.sim.coriolis_omega*op.mu(i_v);
					o_v_t -= 2.0*simVars.sim.coriolis_omega*op.mu(i_u);
				}

				// non-linear terms
				o_u_t -= (i_u*op.grad_lon(i_u) + i_v*op.grad_lat(i_u))*(1.0/simVars.sim.earth_radius);
				o_v_t -= (i_u*op.grad_lon(i_v) + i_v*op.grad_lat(i_v))*(1.0/simVars.sim.earth_radius);
			}
			else
			{
				/*
				 * Height
				 */
				// non-linear equations
				o_h_t = -(op.robert_div_lon(i_h*i_u)+op.robert_div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);

				/*
				 * Velocity
				 */
				// linear terms
				o_u_t = -op.robert_grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.robert_grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					o_u_t += 2.0*simVars.sim.coriolis_omega*op.mu(i_v);
					o_v_t -= 2.0*simVars.sim.coriolis_omega*op.mu(i_u);
				}

				// non-linear terms
				o_u_t -= (i_u*op.robert_grad_lon(i_u) + i_v*op.robert_grad_lat(i_u))*(1.0/simVars.sim.earth_radius);
				o_v_t -= (i_u*op.robert_grad_lon(i_v) + i_v*op.robert_grad_lat(i_v))*(1.0/simVars.sim.earth_radius);

			}
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			double scalar = simVars.sim.viscosity/(simVars.sim.earth_radius*simVars.sim.earth_radius);

			o_h_t += op.laplace(i_h)*scalar;
			o_u_t += op.laplace(i_u)*scalar;
			o_v_t += op.laplace(i_v)*scalar;
		}
	}



private:
	void p_run_timesteps(
			double i_timestep_size,
			double i_max_simulation_time
	)
	{
		// maximum simulation time
		double max_simulation_time = i_max_simulation_time;

		double current_timestep_size = i_timestep_size;

		double current_simulation_time = 0;

		std::cout << "using timestep size: " << current_timestep_size << std::endl;
		std::cout << "max sim time: " << max_simulation_time << std::endl;

		SphereDataTimesteppingExplicitRK timestepping;

		std::cout << max_simulation_time << ", " << current_simulation_time << std::endl;
		while (max_simulation_time > 0.0 && max_simulation_time > current_simulation_time)
		{
			double o_dt;

			timestepping.run_rk_timestep(
					this,
					&GenerateConsistentGradDivSphereData::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					prog_h, prog_u, prog_v,
					o_dt,
					current_timestep_size,
					4,
					current_simulation_time,
					max_simulation_time
				);

			// advance time step and provide information to parameters
			current_timestep_size = o_dt;
			current_simulation_time += o_dt;

			std::cout << "." << std::flush;
		}
	}


private:
	void p_setup_initial_conditions_gaussian(
			double i_center_lat = M_PI/3,
			double i_center_lon = M_PI/3
	)
	{
		double exp_fac = 10.0;

		double center_lat = i_center_lat;
		double center_lon = i_center_lon;

		auto initial_condition_h = [&](double lon, double mu, double &o_data)
		{
			// https://en.wikipedia.org/wiki/Great-circle_distance
			// d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2))
			// exp(-pow(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lon1-lon2)), 2)*A)

			double phi1 = asin(mu);
			double phi2 = center_lat;
			double lambda1 = lon;
			double lambda2 = center_lon;

			double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

			o_data = exp(-d*d*exp_fac)*0.1*simVars.sim.h0 + simVars.sim.h0;
		};

		prog_h.physical_update_lambda_gaussian_grid(initial_condition_h);
		prog_u.physical_set_zero();
		prog_v.physical_set_zero();
	}



public:
	void generate()
	{
		bool load_and_read_from_file = false;

		char *envvar = getenv("SWEET_HACK_123_BLARG");

		if (envvar != nullptr)
		{
			std::cout << envvar << std::endl;
			exit(1);
			load_and_read_from_file = atoi(envvar);
		}

		if (load_and_read_from_file)
		{
			std::cout << "Using SWEET_HACK_123_BLARG hack... loading data from cached files" << std::endl;

			bool ok = true;
			ok &= prog_h.physical_file_load(prefix_string+"_final_prog_h.csv");
			ok &= prog_u.physical_file_load(prefix_string+"_final_prog_u.csv");
			ok &= prog_v.physical_file_load(prefix_string+"_final_prog_v.csv");

			if (ok)
				return;

			std::cerr << "Failed to load file from at least one file... regenerating data" << std::endl;;
		}

		p_setup_initial_conditions_gaussian(center_lon, center_lat);

		if (load_and_read_from_file)
		{
			prog_h.physical_file_write(prefix_string+"_initial_prog_h.csv");
			prog_u.physical_file_write(prefix_string+"_initial_prog_u.csv");
			prog_v.physical_file_write(prefix_string+"_initial_prog_v.csv");
		}

		// NOTE: This is only a very rough approximation of the time step size
		double timestep_size = simVars.sim.earth_radius/(
				std::max(simVars.disc.res_physical[0], simVars.disc.res_physical[1])*
					std::max(simVars.sim.gravitation, std::max(simVars.sim.h0, simVars.sim.coriolis_omega))
		);

		// run 200 time steps
		double simtime = timestep_size * 100;

		p_run_timesteps(timestep_size, simtime);

		std::cout << std::endl;

		// rescale all variables to 1.0
		prog_h = ((prog_h-simVars.sim.h0).physical_rescale_to_max_abs(1.0) + simVars.sim.h0);
		prog_u = prog_u.physical_rescale_to_max_abs(1.0);
		prog_v = prog_v.physical_rescale_to_max_abs(1.0);

		if (load_and_read_from_file)
		{
			prog_h.physical_file_write(prefix_string+"_final_prog_h.csv");
			prog_u.physical_file_write(prefix_string+"_final_prog_u.csv");
			prog_v.physical_file_write(prefix_string+"_final_prog_v.csv");
		}
	}

	GenerateConsistentGradDivSphereData(
			SimulationVariables &i_simVars,

			SphereDataConfig *i_sphereDataConfig,
			SphereOperators &i_op,

			double i_center_lon =  M_PI/3,
			double i_center_lat = M_PI/3,
			std::string i_prefix_string=""
	)	:
		simVars(i_simVars),
		sphereDataConfig(i_sphereDataConfig),
		op(i_op)
	{
		prog_h.setup(sphereDataConfig);
		prog_u.setup(sphereDataConfig);
		prog_v.setup(sphereDataConfig);

		center_lon = i_center_lon;
		center_lat = i_center_lat;

		if (i_prefix_string != "")
			prefix_string = i_prefix_string;
		else
			prefix_string = std::string("gen_divgrad_data_")+sphereDataConfig->getUniqueIDString();
	}
};



#endif /* SRC_INCLUDE_SWEET_SPHERE_GENERATECONSISTENTGRADDIVSPHEREDATA_HPP_ */
