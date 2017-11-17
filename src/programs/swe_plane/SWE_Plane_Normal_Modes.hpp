/*
 * SWE_Plane_TimeSteppers.hpp
 *
 *  Created on: 17 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *      based on previous implementation by Martin Schreiber in swe_plane.cpp
 *
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_
#define SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>


/**
 * SWE Plane normal mode
 */
class SWE_Plane_Normal_Modes
{
public:

	bool linear_only = false;


	SWE_Plane_Normal_Modes()
	{
	}

	void normal_mode_analysis();

	void setup(	)
	{
	};

	~SWE_Plane_Normal_Modes()
	{

	}
};

	void SWE_Plane_Normal_Modes::normal_mode_analysis()
	{
		// dummy time step to get time step size
		if (simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Normal mode analysis requires setting fixed time step size");

		//run_timestep();

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		const char* filename;
		char buffer_real[1024];

		if (simVars.misc.output_file_name_prefix == "")
			filename = "output_%s_normalmodes.csv";
		else
			filename = simVars.misc.output_file_name_prefix.c_str();


		sprintf(buffer_real, filename, "normal_modes_physical", simVars.timecontrol.current_timestep_size*simVars.misc.output_time_scale);
		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		PlaneData* prog[3] = {&prog_h_pert, &prog_u, &prog_v};

		/*
		 * Maximum number of prognostic variables
		 *
		 * Advection e.g. has only one
		 */
		int max_prog_id = -1;
		if (simVars.pde.id == 0 || simVars.pde.id == 1)
		{
			max_prog_id = 3;
			prog_h_pert.physical_set_zero();
			prog_u.physical_set_zero();
			prog_v.physical_set_zero();
		}
		else
		{
			max_prog_id = 1;
			prog_h_pert.physical_set_zero();
		}

#if 0
		if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
		{
			FatalError("Not yet tested and supported");
			std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
			std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
			std::cout << "Therefore, we also halve the time step size here" << std::endl;

			simVars.timecontrol.current_timestep_size = 0.5*simVars.sim.CFL;
			simVars.sim.CFL = -simVars.timecontrol.current_timestep_size;
		}
#endif

		int num_timesteps = 1;
		if (simVars.disc.normal_mode_analysis_generation >= 10)
		{
			if (simVars.timecontrol.max_timesteps_nr > 0)
				num_timesteps = simVars.timecontrol.max_timesteps_nr;
		}

		if (simVars.timecontrol.max_simulation_time > 0)
			file << "# t " << simVars.timecontrol.max_simulation_time << std::endl;
		else
			file << "# t " << (num_timesteps*(-simVars.sim.CFL)) << std::endl;

		file << "# g " << simVars.sim.gravitation << std::endl;
		file << "# h " << simVars.sim.h0 << std::endl;
		file << "# r " << simVars.sim.earth_radius << std::endl;
		file << "# f " << simVars.sim.coriolis_omega << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
		file << "# specnummodes " << specmodes << std::endl;
		file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
		file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

		file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
		file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
		file << "# normalmodegeneration " << simVars.disc.normal_mode_analysis_generation << std::endl;
		file << "# antialiasing ";

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		file << 1;
#else
		file << 0;
#endif

		file << std::endl;


		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			if (simVars.disc.normal_mode_analysis_generation == 1 || simVars.disc.normal_mode_analysis_generation == 11)
			{
				// iterate over physical space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->physical_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					simVars.timecontrol.current_timestep_nr = 0;
					simVars.timecontrol.current_simulation_time = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->physical_set_zero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] = 1;

					/*
					 * RUN timestep
					 */

					run_timestep();

					if (simVars.disc.normal_mode_analysis_generation == 1)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							(*prog[inner_prog_id]) /= simVars.timecontrol.current_timestep_size;
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_physical();
						for (std::size_t k = 0; k < planeDataConfig->physical_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->physical_space_data[k];
							if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->physical_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
#if 1
			else if (simVars.disc.normal_mode_analysis_generation == 3 || simVars.disc.normal_mode_analysis_generation == 13)
			{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
				FatalError("Only available with if plane spectral space is activated during compile time!");
#else

				// iterate over spectral space
				for (int r = 0; r < 2; r++)
				{

					for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
					{
						for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
						{
							// reset time control
							simVars.timecontrol.current_timestep_nr = 0;
							simVars.timecontrol.current_simulation_time = 0;

							std::cout << "." << std::flush;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog[inner_prog_id]->spectral_set_zero();

							// activate mode via real coefficient
							prog[outer_prog_id]->p_spectral_set(j, i, 1.0);

							/*
							 * RUN timestep
							 */
							run_timestep();


							if (simVars.disc.normal_mode_analysis_generation == 3)
							{
								/*
								 * compute
								 * 1/dt * (U(t+1) - U(t))
								 */
								prog[outer_prog_id]->request_data_spectral();

								std::complex<double> val = prog[outer_prog_id]->p_spectral_get(j, i);
								val = val - 1.0;
								prog[outer_prog_id]->p_spectral_set(j, i, val);

								for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
									(*prog[inner_prog_id]) /= simVars.timecontrol.current_timestep_size;
							}


							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							{
								prog[inner_prog_id]->request_data_spectral();

								/*
								 * REAL
								 */

								for (int r = 0; r < 2; r++)
								{
									for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
									{
										for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
										{
											file << prog[inner_prog_id]->p_spectral_get(j, i).real();
											file << "\t";
										}
									}
								}
							}


							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							{
								/*
								 * IMAG
								 */
								int c = 0;
								for (int r = 0; r < 2; r++)
								{
									for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
									{
										for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
										{
											file << prog[inner_prog_id]->p_spectral_get(j, i).imag();

											if (inner_prog_id != max_prog_id-1 || c != specmodes-1)
												file << "\t";
											else
												file << std::endl;

											c++;
										}
									}
								}
							}
						}
					}
				}
#endif
			}
#else
			else if (simVars.disc.normal_mode_analysis_generation == 3 || simVars.disc.normal_mode_analysis_generation == 13)
			{
				PlaneDataComplex t1(planeDataConfig);
				PlaneDataComplex t2(planeDataConfig);
				PlaneDataComplex t3(planeDataConfig);
				PlaneDataComplex* prog_cplx[3] = {&t1, &t2, &t3};

				// iterate over spectral space
				for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_complex_array_data_number_of_elements; outer_i++)
				{
					// reset time control
					simVars.timecontrol.current_timestep_nr = 0;
					simVars.timecontrol.current_simulation_time = 0;

					std::cout << "." << std::flush;

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog_cplx[inner_prog_id]->spectral_set_zero();

					// activate mode via real coefficient
					prog_cplx[outer_prog_id]->request_data_spectral();
					prog_cplx[outer_prog_id]->spectral_space_data[outer_i].real(1);

					// convert PlaneDataComplex to PlaneData
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						*prog[inner_prog_id] = Convert_PlaneDataComplex_To_PlaneData::physical_convert(*prog_cplx[inner_prog_id]);
						prog[inner_prog_id]->spectral_zeroAliasingModes();
					}


					/*
					 * RUN timestep
					 */
					run_timestep();

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->spectral_zeroAliasingModes();
#warning "update this physical_convert maybe to spectral_convert"

						*prog_cplx[inner_prog_id] = Convert_PlaneData_To_PlaneDataComplex::physical_convert(*prog[inner_prog_id]);

						prog_cplx[inner_prog_id]->request_data_spectral();
					}

					if (simVars.disc.normal_mode_analysis_generation == 3)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog_cplx[outer_prog_id]->request_data_spectral();
						prog_cplx[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog_cplx[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);
					}


					// convert PlaneDataComplex to PlaneData
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog_cplx[inner_prog_id]->request_data_spectral();

						/*
						 * REAL
						 */
						for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
						{
							file << prog_cplx[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}

						/*
						 * IMAG
						 */
						for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
						{
							file << prog_cplx[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != planeDataConfig->spectral_complex_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
#endif
		}
	}



#endif /* SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_ */
