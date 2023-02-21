/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/advection_sphere_timeintegrators/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/advection_sphere_benchmarks/
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */


#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include "../programs/advection_sphere_benchmarks/BenchmarksSphereAdvection.hpp"
#include "../programs/advection_sphere_timeintegrators/SphereAdvectionTimeSteppers.hpp"
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>



// Sphere data config
SphereData_Config sphereDataConfigInstance;
SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;


SimulationVariables simVars;


class SimulationInstance
{
public:
	std::vector<SphereData_Spectral*> prognostic_variables;
	std::vector<SphereData_Spectral*> prognostic_variables_t0;

	//SphereData_Spectral prog_vort, prog_div;
	SphereData_Physical velocity_field_u, velocity_field_v;

	SphereAdvectionTimeSteppers timeSteppers;

	SphereOperators_SphereData op;

	bool time_varying_fields;

	/*
	 * LMax error to h0
	 */
	double max_error_h0;
	/*
	 * RMS error to h0
	 */
	double rms_error_h0;



	BenchmarksSphereAdvection sphereBenchmarks;



public:
	SimulationInstance()	:
		velocity_field_u(sphereDataConfig),
		velocity_field_v(sphereDataConfig),

		op(sphereDataConfig, &(simVars.sim)),
		time_varying_fields(false)
	{
		reset();
	}


	void alloc_prognostic_variables(std::size_t i_size)
	{
		prognostic_variables.resize(i_size);
		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			prognostic_variables[i] = new SphereData_Spectral(sphereDataConfig);

		prognostic_variables_t0.resize(i_size);
		for (std::size_t i = 0; i < prognostic_variables_t0.size(); i++)
			prognostic_variables_t0[i] = new SphereData_Spectral(sphereDataConfig);
	}


	void free_prognostic_variables()
	{
		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			delete prognostic_variables[i];
		prognostic_variables.clear();

		for (std::size_t i = 0; i < prognostic_variables_t0.size(); i++)
			delete prognostic_variables_t0[i];
		prognostic_variables_t0.clear();
	}


	void reset()
	{
		simVars.reset();

		sphereBenchmarks.setup(simVars, op);

		int num_field_variables = sphereBenchmarks.master->get_num_prognostic_fields();

		free_prognostic_variables();
		alloc_prognostic_variables(num_field_variables);

		sphereBenchmarks.master->get_initial_state(
				prognostic_variables, velocity_field_u, velocity_field_v
			);

		for (std::size_t i = 0; i < prognostic_variables.size(); i++)
			*prognostic_variables_t0[i] = *prognostic_variables[i];

		time_varying_fields = sphereBenchmarks.master->has_time_varying_state();

		// setup sphereDataconfig instance again
		sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans, simVars.misc.verbosity, simVars.parallelization.num_threads_space);

		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);
	}



	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;


		timeSteppers.master->run_timestep(
				prognostic_variables, velocity_field_u, velocity_field_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time,
				(time_varying_fields ? &sphereBenchmarks : nullptr)
			);

		double dt = simVars.timecontrol.current_timestep_size;

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.misc.verbosity >= 10)
			std::cout << simVars.timecontrol.current_timestep_nr << ": " << simVars.timecontrol.current_simulation_time/(60*60*24.0) << std::endl;

#if 0
		max_error_h0 = (prog_phi_pert_t0-prog_phi_pert).toPhys().physical_reduce_max_abs();
		rms_error_h0 = (prog_phi_pert_t0-prog_phi_pert).toPhys().physical_reduce_rms();
#else
		max_error_h0 = -1;
		rms_error_h0 = -1;
#endif
	}


	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		double diff = std::abs(simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time);

		if (	simVars.timecontrol.max_simulation_time != -1 &&
			(
					simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time	||
					diff/simVars.timecontrol.max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
			)
		)
			return true;

		return false;
	}


};



int main(int i_argc, char *i_argv[])
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		std::cout << std::endl;
		return -1;
	}

	simVars.outputConfig();

	int initial_spectral_modes = simVars.disc.space_res_spectral[0];

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	int max_modes = 256;

	if (simVars.disc.timestepping_order == 1)
		max_modes = 512;
	else if (simVars.disc.timestepping_order == 2)
		max_modes = 256;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		simVars.timecontrol.current_timestep_size *= 0.5;
		simVars.timecontrol.setup_timestep_size = simVars.timecontrol.current_timestep_size;

		if (simVars.disc.timestepping_method == "na_sl")
		{
			simVars.disc.space_res_spectral[0] = i;
			simVars.disc.space_res_spectral[1] = i;

			simVars.disc.space_res_physical[0] = 2*i;
			simVars.disc.space_res_physical[1] = i;
		}
		else
		{
			simVars.disc.space_res_spectral[0] = initial_spectral_modes;
			simVars.disc.space_res_spectral[1] = initial_spectral_modes;

			simVars.disc.space_res_physical[0] = 0;
			simVars.disc.space_res_physical[1] = 0;
		}


		sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans, simVars.misc.verbosity, simVars.parallelization.num_threads_space);

		if (1)
		{
			std::cout << "Checking for right velocity" << std::endl;
			SphereTimestepping_SemiLagrangian sl(simVars, sphereDataConfig);


			/*
			 * Convert to Cartesian velocity space
			 */
			ScalarDataArray V_lon_D(sphereDataConfig->physical_array_data_number_of_elements);
			ScalarDataArray V_lat_D(sphereDataConfig->physical_array_data_number_of_elements);
			V_lon_D.set_all(1.0);
			V_lat_D.set_all(2.0);

			ScalarDataArray V_x_D, V_y_D, V_z_D;
			SWEETVectorMath::velocity_latlon_to_cartesian__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_lon_D,
					V_lat_D,
					V_x_D,
					V_y_D,
					V_z_D
				);

			ScalarDataArray V_lon_tmp, V_lat_tmp;
			SWEETVectorMath::velocity_cartesian_to_latlon__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_x_D,
					V_y_D,
					V_z_D,
					V_lon_tmp, V_lat_tmp
			);

			double err_lon = (V_lon_D - V_lon_tmp).reduce_maxAbs();
			double err_lat = (V_lat_D - V_lat_tmp).reduce_maxAbs();

			if (err_lon > 1e-10)
			{
				std::cerr << "Error: " << err_lon << std::endl;
				SWEETError("Error lon too high!");
			}

			if (err_lat > 1e-10)
			{
				std::cerr << "Error: " << err_lat << std::endl;
				SWEETError("Error lat too high!");
			}
		}

		std::cout << "Testing with " << sphereDataConfigInstance.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << simVars.timecontrol.current_timestep_size << std::endl;

		SimulationInstance simulation;


		{
			while (!simulation.should_quit())
				simulation.run_timestep();

			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << simulation.max_error_h0 << std::endl;
			std::cout << "RMS error: " << simulation.rms_error_h0 << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / simulation.max_error_h0;
				std::cout << "Convergence: " << conv << std::endl;

				if (conv*1.1 < std::pow(2.0, (double)simVars.disc.timestepping_order))
					SWEETError("Convergence not given!");
			}

			if (simulation.max_error_h0  > 1e10)
				SWEETError("Lmax error exceeded threshold!");

			prev_max_error = simulation.max_error_h0;

			std::cout << "*********************************************" << std::endl;
		}
	}

	return 0;
}
