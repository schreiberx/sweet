/*
 * GenerateConsistentGradDivSphereData.hpp
 *
 *  Created on: 4 Nov 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SPHERE_GENERATECONSISTENTGRADDIVSPHEREDATA_HPP_
#define SRC_INCLUDE_SWEET_SPHERE_GENERATECONSISTENTGRADDIVSPHEREDATA_HPP_


#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>
#include "../../../programs/swe_sphere/SWE_Sphere_TS_l_erk.hpp"



/**
 * This class generates data which is consistent with SPH div/grad operators.
 *
 * One of the main problems by using SPH is, that the data has to fulfill
 * certain consistency requirements.
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

		timestepping.resetAndSetup(sphereDataConfig, 4);

		SWE_Sphere_TS_l_erk l_erk(simVars, op);
		l_erk.setup(4);

		SphereDataPhysical prog_u_g = prog_u.getSphereDataPhysical();
		SphereDataPhysical prog_v_g = prog_v.getSphereDataPhysical();

		SphereData phi = prog_h * simVars.sim.gravitation;
		SphereData vort(sphereDataConfig);
		SphereData div(sphereDataConfig);

		op.robert_uv_to_vortdiv(prog_u_g, prog_v_g, vort, div);

		std::cout << max_simulation_time << ", " << current_simulation_time << std::endl;
		while (max_simulation_time > 0.0 && max_simulation_time > current_simulation_time)
		{
			l_erk.run_timestep(phi, vort, div, current_timestep_size, current_simulation_time);

			// advance time step and provide information to parameters
			current_simulation_time += current_timestep_size;

			std::cout << "." << std::flush;
		}

		prog_h = phi/simVars.sim.gravitation;
		op.robert_vortdiv_to_uv(vort, div, prog_u_g,prog_v_g);
		prog_u = prog_u_g;
		prog_v = prog_v_g;
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
		)*0.1;

		// run 200 time steps
		double simtime = timestep_size * 10;

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
