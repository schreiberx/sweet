/*
 * SWE_Plane_Normal_Modes.hpp
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
#include <functional>
#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif
/**
 * SWE Plane normal mode
 */
class SWE_Plane_Normal_Modes
{
public:
	template <typename TCallbackClass>
	static
	void normal_mode_analysis(
			PlaneData &io_prog_h_pert, // h: surface height (perturbation)
			PlaneData &io_prog_u, // u: velocity in x-direction
			PlaneData &io_prog_v, // v: velocity in y-direction
			SimulationVariables &i_simVars, // Simulation variables
			TCallbackClass *i_class,
			void(TCallbackClass::* const i_run_timestep_method)(void)
	)
	{

		const PlaneDataConfig *planeDataConfig = io_prog_h_pert.planeDataConfig;

		// dummy time step to get time step size
		if (i_simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Normal mode analysis requires setting fixed time step size");

		/*
		 *
		 * Mode-wise normal mode analysis
		 *
		 *
		 */

		if (i_simVars.disc.normal_mode_analysis_generation == 4)
		{
#if !SWEET_EIGEN
			FatalError("SWE_Plane_Normal_Modes: Cannot test this without Eigen library. Please compile with --eigen=enable");
#endif
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
			FatalError("SWE_Plane_Normal_Modes: This test was build for linear or linearized models, so please compile without dealising --plane-spectral-dealiasing=disable.");
#endif


			/*
			 * Setup all output files
			 */
			const char* filename; //general filename
			char buffer_real[1024];

			if (i_simVars.misc.output_file_name_prefix == "")
				filename = "output_%s_t%020.8f.csv";
			else
				filename = i_simVars.misc.output_file_name_prefix.c_str();

			sprintf(buffer_real, filename, "normal_modes_plane", i_simVars.timecontrol.current_timestep_size*i_simVars.misc.output_time_scale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to files of the form '" << buffer_real << "'" << std::endl;

			//Positive inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_plane_igpos", i_simVars.timecontrol.current_timestep_size*i_simVars.misc.output_time_scale);
			std::ofstream file_igpos(buffer_real, std::ios_base::trunc);

			//Negative inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_plane_igneg", i_simVars.timecontrol.current_timestep_size*i_simVars.misc.output_time_scale);
			std::ofstream file_igneg(buffer_real, std::ios_base::trunc);

			//Geostrophic modes
			sprintf(buffer_real, filename, "normal_modes_plane_geo", i_simVars.timecontrol.current_timestep_size*i_simVars.misc.output_time_scale);
			std::ofstream file_geo(buffer_real, std::ios_base::trunc);

			//std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);
			file_igpos << std::setprecision(20);
			file_igneg << std::setprecision(20);
			file_geo << std::setprecision(20);

			file << "# dt " << i_simVars.timecontrol.current_timestep_size << std::endl;
			file << "# g " << i_simVars.sim.gravitation << std::endl;
			file << "# h " << i_simVars.sim.h0 << std::endl;
			file << "# r " << i_simVars.sim.earth_radius << std::endl;
			file << "# f " << i_simVars.sim.f0 << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
			file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

			file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
			file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
			file << "# normalmodegeneration " << i_simVars.disc.normal_mode_analysis_generation << std::endl;
			file << "# antialiasing ";
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
			file << 1;
#else
			file << 0;
#endif
			file << std::endl;

			PlaneData* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			int max_prog_id = 3;
			//The basic state is with zero in all variables
			// The only non zero variable in the basic state is the total height
			//    for which the constant is added within run_timestep()
			io_prog_h_pert.physical_set_zero();
			io_prog_u.physical_set_zero();
			io_prog_v.physical_set_zero();

			//int num_timesteps = 1;

			// Timestep and perturbation
			double dt = i_simVars.timecontrol.current_timestep_size;
			double eps = dt;

			//Matrix representing discrete linear operator in spectral space
			Eigen::MatrixXcf A(3,3) ;
			//Eigen solver
			Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
			//Final eigenvalues
			std::complex<double> eval[3];

			//For each spectral mode
			//for (int r = 0; r < 2; r++) //only required to get the symmetric half of the spectrum
			//{
			int r = 0;

			for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
			{
				std::cout << "." << std::flush;
				for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					//This is the mode to be analysed
					//std::cout << "Mode (i,j)= (" << i << " , " << j <<")"<< std::endl;


					for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
					{

						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->spectral_set_zero();

						// activate mode via real coefficient
						prog[outer_prog_id]->p_spectral_set(j, i, 1.0);
						//Activate the symetric couterpart of the mode (only needed if j>0 )
						if (j > 0)
							prog[outer_prog_id]->p_spectral_set(planeDataConfig->spectral_data_size[1]-j, i, 1.0);

						/*
						 * RUN timestep
						 */
						prog[outer_prog_id]->request_data_physical();
						(i_class->*i_run_timestep_method)();

						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_spectral();

						std::complex<double> val = prog[outer_prog_id]->p_spectral_get(j, i);
						val = val - 1.0; //subtract U(0) from mode
						prog[outer_prog_id]->p_spectral_set(j, i, val);

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							(*prog[inner_prog_id]) /= eps;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							A(inner_prog_id,outer_prog_id)=prog[inner_prog_id]->p_spectral_get(j, i);;
						}

					}

					//std::cout << "Lik matrix" << std::endl;
					//std::cout << A << std::endl;

					//std::cout<<"Normal modes" << std::endl;
					ces.compute(A);
					for(int i=0; i<3; i++)
					{
						eval[i]=ces.eigenvalues()[i];
						//std::cout << "Eigenvalue "<< i << " : " << eval[i].real() <<" "<<eval[i].imag() << std::endl;

					}
					/* We will try to separate the modes in 3 types:
					 * -positive inertia-gravity (imag>f) - we will adopt coriolis f to test as if > zero, since the exact freq is sqrt(f^2+cK*K)
					 * -negative inertia-gravity (imag<-f)
					 * -negative inertia-gravity (imag aprox 0) - we will fit all other modes here
					 */
					for(int i=0; i<3; i++)
					{
						if(eval[i].imag() > 0.5 * i_simVars.sim.f0)
						{
							//std::cout<< "IG pos mode: " << eval[i].imag() << std::endl;
							//file_igpos << eval[i].imag();
							file_igpos << eval[i];
							file_igpos << "\t";
						}
						if(eval[i].imag() < - 0.5 * i_simVars.sim.f0)
						{
							//std::cout<< "IG neg mode: " << eval[i].imag() << std::endl;
							//file_igneg << eval[i].imag();
							file_igneg << eval[i];
							file_igneg << "\t";
						}
						if(eval[i].imag() >= - 0.5 * i_simVars.sim.f0 && eval[i].imag() <=  0.5 * i_simVars.sim.f0 )
						{
							//std::cout<< "IG geo mode: " << eval[i].imag() << std::endl;
							//file_geo << eval[i].imag();
							file_geo << eval[i];
							file_geo << "\t";
						}
					}
					//std::cout<<"-------------------------" << std::endl;
				}
				file_igpos << std::endl;
				file_igneg << std::endl;
				file_geo << std::endl;
			}

			//}
			//std::cout<<"-------------------------" << std::endl;
			//FatalError("still needs work...");
		}
		/*
		 * Do a normal mode analysis using perturbation, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		else
		{

			//run_timestep();
			const char* filename;
			char buffer_real[1024];

			if (i_simVars.misc.output_file_name_prefix == "")
				filename = "output_%s_normalmodes.csv";
			else
				filename = i_simVars.misc.output_file_name_prefix.c_str();


			sprintf(buffer_real, filename, "normal_modes_physical", i_simVars.timecontrol.current_timestep_size*i_simVars.misc.output_time_scale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

			std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);

			PlaneData* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			/*
			 * Maximum number of prognostic variables
			 *
			 * Advection e.g. has only one
			 */
			int max_prog_id = -1;
			if (i_simVars.pde.id == 0 || i_simVars.pde.id == 1)
			{
				max_prog_id = 3;
				io_prog_h_pert.physical_set_zero();
				io_prog_u.physical_set_zero();
				io_prog_v.physical_set_zero();
			}
			else
			{
				max_prog_id = 1;
				io_prog_h_pert.physical_set_zero();
			}

#if 0
			if (i_simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
			{
				FatalError("Not yet tested and supported");
				std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
				std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
				std::cout << "Therefore, we also halve the time step size here" << std::endl;

				i_simVars.timecontrol.current_timestep_size = 0.5*i_simVars.sim.CFL;
				i_simVars.sim.CFL = -i_simVars.timecontrol.current_timestep_size;
			}
#endif

			int num_timesteps = 1;
			if (i_simVars.disc.normal_mode_analysis_generation >= 10)
			{
				if (i_simVars.timecontrol.max_timesteps_nr > 0)
					num_timesteps = i_simVars.timecontrol.max_timesteps_nr;
			}

			if (i_simVars.timecontrol.max_simulation_time > 0)
				file << "# t " << i_simVars.timecontrol.max_simulation_time << std::endl;
			else
				file << "# t " << (num_timesteps*(-i_simVars.sim.CFL)) << std::endl;

			file << "# g " << i_simVars.sim.gravitation << std::endl;
			file << "# h " << i_simVars.sim.h0 << std::endl;
			file << "# r " << i_simVars.sim.earth_radius << std::endl;
			file << "# f " << i_simVars.sim.coriolis_omega << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
			file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

			file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
			file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
			file << "# normalmodegeneration " << i_simVars.disc.normal_mode_analysis_generation << std::endl;
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
				if (i_simVars.disc.normal_mode_analysis_generation == 1 || i_simVars.disc.normal_mode_analysis_generation == 11)
				{
					// iterate over physical space
					for (std::size_t outer_i = 0; outer_i < planeDataConfig->physical_array_data_number_of_elements; outer_i++)
					{
						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

						std::cout << "." << std::flush;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->physical_set_zero();

						// activate mode
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] = 1;

						/*
						 * RUN timestep
						 */

						(i_class->*i_run_timestep_method)();

						if (i_simVars.disc.normal_mode_analysis_generation == 1)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog[outer_prog_id]->request_data_physical();
							prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								(*prog[inner_prog_id]) /= i_simVars.timecontrol.current_timestep_size;
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
				else if (i_simVars.disc.normal_mode_analysis_generation == 3 || i_simVars.disc.normal_mode_analysis_generation == 13)
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
								i_simVars.timecontrol.current_timestep_nr = 0;
								i_simVars.timecontrol.current_simulation_time = 0;

								std::cout << "." << std::flush;

								for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
									prog[inner_prog_id]->spectral_set_zero();

								// activate mode via real coefficient
								prog[outer_prog_id]->p_spectral_set(j, i, 1.0);

								/*
								 * RUN timestep
								 */
								(i_class->*i_run_timestep_method)();


								if (i_simVars.disc.normal_mode_analysis_generation == 3)
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
										(*prog[inner_prog_id]) /= i_simVars.timecontrol.current_timestep_size;
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
				else if (i_simVars.disc.normal_mode_analysis_generation == 3 || i_simVars.disc.normal_mode_analysis_generation == 13)
				{
					PlaneDataComplex t1(planeDataConfig);
					PlaneDataComplex t2(planeDataConfig);
					PlaneDataComplex t3(planeDataConfig);
					PlaneDataComplex* prog_cplx[3] = {&t1, &t2, &t3};

					// iterate over spectral space
					for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_complex_array_data_number_of_elements; outer_i++)
					{
						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

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
						(i_class->*i_run_timestep_method)();

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							prog[inner_prog_id]->spectral_zeroAliasingModes();
#warning "update this physical_convert maybe to spectral_convert"

							*prog_cplx[inner_prog_id] = Convert_PlaneData_To_PlaneDataComplex::physical_convert(*prog[inner_prog_id]);

							prog_cplx[inner_prog_id]->request_data_spectral();
						}

						if (i_simVars.disc.normal_mode_analysis_generation == 3)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog_cplx[outer_prog_id]->request_data_spectral();
							prog_cplx[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog_cplx[inner_prog_id]->operator*=(1.0/i_simVars.timecontrol.current_timestep_size);
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
	}


	~SWE_Plane_Normal_Modes()
	{

	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_ */
