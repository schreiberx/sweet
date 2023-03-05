/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMEINTEGRATORS_SWE_SPHERE_NORMALMODEANALYSIS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMEINTEGRATORS_SWE_SPHERE_NORMALMODEANALYSIS_HPP_

#include <sweet/core/sphere/SphereData_Config.hpp>



class NormalModeAnalysisSphere
{

public:
	template <typename TCallbackClass>
	static
	void normal_mode_analysis(
			sweet::SphereData_Spectral &io_prog_phi,
			sweet::SphereData_Spectral &io_prog_vort,
			sweet::SphereData_Spectral &io_prog_div,

			sweet::ShackDictionary &i_shackDict,

			TCallbackClass *i_class,
			void(TCallbackClass::* const i_run_timestep_method)(void)
	)
	{
		const sweet::SphereDataConfig *sphereDataConfig = io_prog_phi.sphereDataConfig;

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		char buffer_real[1024];
		const char* filename = i_shackDict.iodata.output_file_name.c_str();

		if (i_shackDict.timecontrol.max_timesteps_nr > 0)
			sprintf(buffer_real, filename, "normal_modes_physical", i_shackDict.timecontrol.current_timestep_size*i_shackDict.timecontrol.max_timesteps_nr*i_shackDict.iodata.output_time_scale);
		else
			sprintf(buffer_real, filename, "normal_modes_physical", i_shackDict.timecontrol.current_timestep_size*i_shackDict.iodata.output_time_scale);

		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		sweet::SphereData_Spectral* prog[3] = {&io_prog_phi, &io_prog_vort, &io_prog_div};

		int max_prog_id = 3;
		prog[0]->spectral_set_zero();
		prog[1]->spectral_set_zero();
		prog[2]->spectral_set_zero();

		int leapfrog_start_num_timesteps = 16;
		double leapfrog_start_timesteps_size;
		double leapfrog_end_timestep_size;
		double leapfrog_original_timestep_size;

		if (i_shackDict.disc.timestepping_method.find("_lf") != std::string::npos)
		{
			SWEETError("TODO: Get this Leapfrog running");

			std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
			std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
			std::cout << "Therefore, we also halve the time step size here" << std::endl;

			leapfrog_original_timestep_size = i_shackDict.timecontrol.current_timestep_size;
			// do 2 leapfrog time steps -> half time step size
			i_shackDict.timecontrol.current_timestep_size *= 0.5;

			leapfrog_start_timesteps_size = i_shackDict.timecontrol.current_timestep_size/(double)leapfrog_start_num_timesteps;
			leapfrog_end_timestep_size = i_shackDict.timecontrol.current_timestep_size;
		}

		int num_timesteps = 1;
		if (i_shackDict.misc.normal_mode_analysis_generation >= 10)
		{
			if (i_shackDict.timecontrol.max_timesteps_nr > 0)
				num_timesteps = i_shackDict.timecontrol.max_timesteps_nr;
		}

		if (i_shackDict.timecontrol.max_simulation_time > 0)
			file << "# t " << i_shackDict.timecontrol.max_simulation_time << std::endl;
		else
			file << "# t " << (num_timesteps*i_shackDict.timecontrol.current_timestep_size) << std::endl;

		file << "# g " << i_shackDict.sim.gravitation << std::endl;
		file << "# h " << i_shackDict.sim.h0 << std::endl;
		file << "# r " << i_shackDict.sim.sphere_radius << std::endl;
		file << "# f " << i_shackDict.sim.sphere_rotating_coriolis_omega << std::endl;

		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			std::cout << "normal mode analysis for prog " << outer_prog_id << std::endl;

			if (i_shackDict.misc.normal_mode_analysis_generation == 1 || i_shackDict.misc.normal_mode_analysis_generation == 11)
			{
				SWEETError("Not supported anymore");
#if 0
				// iterate over physical space
				for (int outer_i = 0; outer_i < sphereDataConfig->physical_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					i_shackDict.timecontrol.current_timestep_nr = 0;
					i_shackDict.timecontrol.current_simulation_time = 0;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_set_zero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] = 1;

					/*
					 * RUN timestep
					 */

					// In case of a multi-step scheme, reset it!
					if (i_shackDict.disc.timestepping_method.find("_lf") != std::string::npos)
					{
						SWEETError("TODO 01943934");
						//spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, i_shackDict.disc.timestepping_order, i_shackDict.disc.leapfrog_robert_asselin_filter);

						i_shackDict.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

						for (int i = 0; i < leapfrog_start_num_timesteps; i++)
							(i_class->*i_run_timestep_method)();

						i_shackDict.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

						(i_class->*i_run_timestep_method)();

						i_shackDict.timecontrol.current_timestep_size = leapfrog_original_timestep_size;
					}
					else
					{
						(i_class->*i_run_timestep_method)();
					}

					for (int i = 1; i < num_timesteps; i++)
					{
						(i_class->*i_run_timestep_method)();
					}


					if (i_shackDict.misc.normal_mode_analysis_generation == 1)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/i_shackDict.timecontrol.current_timestep_size);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_physical();
						for (int k = 0; k < sphereDataConfig->physical_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->physical_space_data[k];
							if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->physical_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
#endif
			}
			else if (i_shackDict.misc.normal_mode_analysis_generation == 2 || i_shackDict.misc.normal_mode_analysis_generation == 12)
			{

				// iterate over physical space
				for (int outer_i = 0; outer_i < sphereDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					for (int imag_i = 0; imag_i < 2; imag_i++)
					{
						// reset time control
						i_shackDict.timecontrol.current_timestep_nr = 0;
						i_shackDict.timecontrol.current_simulation_time = 0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->spectral_set_zero();

						// activate mode
						if (imag_i)
							prog[outer_prog_id]->spectral_space_data[outer_i].imag(1);
						else
							prog[outer_prog_id]->spectral_space_data[outer_i].real(1);

						// In case of a multi-step scheme, reset it!
						if (i_shackDict.disc.timestepping_method.find("_lf") != std::string::npos)
						{
							SWEETError("TODO 01943934");
							//spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, i_shackDict.disc.timestepping_order, i_shackDict.disc.leapfrog_robert_asselin_filter);

							i_shackDict.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

							for (int i = 0; i < leapfrog_start_num_timesteps; i++)
								(i_class->*i_run_timestep_method)();

							i_shackDict.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

							(i_class->*i_run_timestep_method)();

							i_shackDict.timecontrol.current_timestep_size = leapfrog_original_timestep_size;
						}
						else
						{
							(i_class->*i_run_timestep_method)();
						}

						for (int i = 1; i < num_timesteps; i++)
							(i_class->*i_run_timestep_method)();



						if (i_shackDict.misc.normal_mode_analysis_generation == 2)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog[inner_prog_id]->operator*=(1.0/i_shackDict.timecontrol.current_timestep_size);
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].imag();
								if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->spectral_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
			}
			else if (i_shackDict.misc.normal_mode_analysis_generation == 3 || i_shackDict.misc.normal_mode_analysis_generation == 13)
			{
				// iterate over spectral space
				for (int outer_i = 0; outer_i < sphereDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					i_shackDict.timecontrol.current_timestep_nr = 0;
					i_shackDict.timecontrol.current_simulation_time = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_set_zero();

					// activate mode via real coefficient
					prog[outer_prog_id]->spectral_space_data[outer_i].real(1);


					// In case of a multi-step scheme, reset it!
					if (i_shackDict.disc.timestepping_method.find("_lf") != std::string::npos)
					{
						SWEETError("TODO 01839471");
						//spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, i_shackDict.disc.timestepping_order, i_shackDict.disc.leapfrog_robert_asselin_filter);

						i_shackDict.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

						for (int i = 0; i < leapfrog_start_num_timesteps; i++)
							(i_class->*i_run_timestep_method)();

						i_shackDict.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

						(i_class->*i_run_timestep_method)();

						i_shackDict.timecontrol.current_timestep_size = leapfrog_original_timestep_size;

						if (num_timesteps > 1)
							SWEETError("Doesn't make sense because the previous time step is half the time step size in advance");
					}
					else
					{
						(i_class->*i_run_timestep_method)();
					}

					for (int i = 1; i < num_timesteps; i++)
					{
						(i_class->*i_run_timestep_method)();
					}

					if (i_shackDict.misc.normal_mode_analysis_generation == 3)
					{
						/*
						 * Compute
						 *    1/dt * (U(t+1) - U(t))
						 * for linearization
						 */
						prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/i_shackDict.timecontrol.current_timestep_size);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->spectral_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
		}
	}
};



#endif
