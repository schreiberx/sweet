/*
 * PararealController.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_


#include <parareal/Parareal_ConsolePrefix.hpp>
#include <parareal/Parareal_SimulationInstance_GenericData.hpp>
#include <parareal/Parareal_SimulationVariables.hpp>

#include <parareal/Parareal_GenericData.hpp>
#include <parareal/Parareal_GenericData_Scalar.hpp>
#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

/**
 * This class takes over the control and
 * calls methods offered via PararealSimulation.
 *
 * \param t_SimulationInstance	class which implements the Parareal_SimulationInstance interfaces
 */
////template <class t_SimulationInstance>
template <class t_tsmType, int N>
class Parareal_Controller_Serial_GenericData
{
	/**
	 * Array with instantiations of PararealSimulations
	 */
////	t_SimulationInstance *simulationInstances = nullptr;

	/**
	 * Pointers to interfaces of simulationInstances
	 * This helps to clearly separate between the allocation of the simulation classes and the parareal interfaces.
	 */
	std::vector<Parareal_SimulationInstance_GenericData<t_tsmType, N>*> parareal_simulationInstances = {};

	SimulationVariables* simVars;

	PlaneDataConfig* planeDataConfig = nullptr;
	SphereData_Config* sphereDataConfig = nullptr;

	// Operators
	PlaneOperators op_plane;
	SphereOperators_SphereData op_sphere;
	SphereOperators_SphereData op_sphere_nodealiasing;

	// Geometry (scalar, plane, sphere)
	std::string geometry = "";

	// Model (ODE, Burgers, SWE)
	std::string model = "";

	t_tsmType* timeSteppersFine = nullptr;
	t_tsmType* timeSteppersCoarse = nullptr;


	/**
	 * Pointer to parareal simulation variables.
	 * These variables are used as a singleton
	 */
	Parareal_SimulationVariables *pVars;

	/**
	 * Class which helps prefixing console output
	 */
	Parareal_ConsolePrefix CONSOLEPREFIX;




public:

	// Scalar
	Parareal_Controller_Serial_GenericData(SimulationVariables i_simVars, std::string i_geometry, std::string i_model,
						t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse):
		simVars(&i_simVars), geometry(i_geometry), model(i_model), timeSteppersFine(i_timeSteppersFine), timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

	// Plane
	Parareal_Controller_Serial_GenericData(SimulationVariables& i_simVars, PlaneDataConfig* i_planeDataConfig, PlaneOperators &i_op_plane, std::string i_geometry, std::string i_model,
						t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse):
		simVars(&i_simVars), planeDataConfig(i_planeDataConfig), op_plane(i_op_plane), geometry(i_geometry), model(i_model), timeSteppersFine(i_timeSteppersFine), timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

	// Sphere
	Parareal_Controller_Serial_GenericData(SimulationVariables& i_simVars, SphereData_Config* i_sphereDataConfig, SphereOperators_SphereData &i_op_sphere, SphereOperators_SphereData &i_op_sphere_nodealiasing, std::string i_geometry, std::string i_model,
						t_tsmType* i_timeSteppersFine, t_tsmType* i_timeSteppersCoarse):
		simVars(&i_simVars), sphereDataConfig(i_sphereDataConfig), op_sphere(i_op_sphere), op_sphere_nodealiasing(i_op_sphere_nodealiasing), geometry(i_geometry), model(i_model), timeSteppersFine(i_timeSteppersFine), timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};


///	Parareal_Controller_Serial()
///	{
///	}


	~Parareal_Controller_Serial_GenericData()
	{
		cleanup();
	}


	inline
	void CONSOLEPREFIX_start(
			const char *i_prefix
	)
	{
//		if (pVars->verbosity > 0)
			CONSOLEPREFIX.start(i_prefix);
	}


	inline
	void CONSOLEPREFIX_start(int i_prefix)
	{
//		if (pVars->verbosity > 0)
			CONSOLEPREFIX.start(i_prefix);
	}

	inline
	void CONSOLEPREFIX_end()
	{
		if (pVars->verbosity > 0)
			CONSOLEPREFIX.end();
	}


	void cleanup()
	{
		for (typename std::vector<Parareal_SimulationInstance_GenericData<t_tsmType, N>*>::iterator it = this->parareal_simulationInstances.begin();
															it != this->parareal_simulationInstances.end();
															it++)
			if (*it)
				delete *it;
	}


	void setup(
			Parareal_SimulationVariables *i_pararealSimVars
	)
	{

		cleanup();

		pVars = i_pararealSimVars;

		if (!pVars->enabled)
			return;

		if (pVars->coarse_slices <= 0)
		{
			std::cerr << "Invalid number of coarse slices" << std::endl;
			exit(1);
		}

		if (pVars->max_simulation_time <= 0)
		{
			std::cerr << "Invalid simulation time" << std::endl;
			exit(1);
		}

		// allocate raw simulation instances
		////simulationInstances = new t_SimulationInstance[pVars->coarse_slices];

		CONSOLEPREFIX.start("[MAIN] ");
		std::cout << "Resetting simulation instances" << std::endl;

		// convert to pararealsimulationInstances to get Parareal interfaces
		for (int k = 0; k < pVars->coarse_slices; k++)
		{
			CONSOLEPREFIX_start(k);
		
			parareal_simulationInstances.push_back(new Parareal_SimulationInstance_GenericData<t_tsmType, N>);
		///	parareal_simulationInstances[k] = &(Parareal_SimulationInstance&)(simulationInstances[k]);
			if ( ! this->geometry.compare("scalar") )
				parareal_simulationInstances[k]->setup(this->simVars,
								       this->geometry,
								       this->model,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);
			else if ( !this->geometry.compare("plane") )
				parareal_simulationInstances[k]->setup(this->simVars,
								       this->planeDataConfig,
								       this->geometry,
								       this->model,
								       &this->op_plane,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);
			else if ( !this->geometry.compare("sphere") )
				parareal_simulationInstances[k]->setup(this->simVars,
								       this->sphereDataConfig,
								       this->geometry,
								       this->model,
								       &this->op_sphere,
								       &this->op_sphere_nodealiasing,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);
			else
				SWEETError("Unknown geometry");
		}


		CONSOLEPREFIX_start("[MAIN] ");
		std::cout << "Setup time frames" << std::endl;

		/*
		 * SETUP time frame
		 */
		// size of coarse time step

		double time_slice_size = pVars->max_simulation_time / pVars->coarse_slices;

		if (pVars->coarse_timestep_size < 0)
			pVars->coarse_timestep_size = time_slice_size;

		// if time slices are not homogeneous, this should be called by each parareal_simulationInstance
		//for (int k = 0; k < pVars->coarse_slices; k++)
		//{
			CONSOLEPREFIX_start(0);
			parareal_simulationInstances[0]->sim_check_timesteps(time_slice_size);
		//}

		CONSOLEPREFIX_start(0);
		parareal_simulationInstances[0]->sim_set_timeframe(0, time_slice_size);

		for (int k = 1; k < pVars->coarse_slices-1; k++)
		{
			CONSOLEPREFIX_start(k);
			parareal_simulationInstances[k]->sim_set_timeframe(time_slice_size*k, time_slice_size*(k+1));
		}

		CONSOLEPREFIX_start(pVars->coarse_slices-1);
		parareal_simulationInstances[pVars->coarse_slices-1]->sim_set_timeframe(i_pararealSimVars->max_simulation_time-time_slice_size, i_pararealSimVars->max_simulation_time);


		/*
		 * Setup first simulation instance
		 */
		CONSOLEPREFIX_start(0);
		parareal_simulationInstances[0]->sim_setup_initial_data();
		CONSOLEPREFIX_end();

		CONSOLEPREFIX_start("[MAIN] ");
		std::cout << "Finished setup parareal" << std::endl;
		CONSOLEPREFIX_end();
	}

	void run()
	{


#if SWEET_DEBUG
		// DEBUG: full fine simulation
		CONSOLEPREFIX_start(0);
		parareal_simulationInstances[0]->run_timestep_fine();
		*(parareal_simulationInstances[0]->parareal_data_fine_exact_debug) = *(parareal_simulationInstances[0]->parareal_data_fine);
		for (int i = 1; i < pVars->coarse_slices; i++)
		{
			CONSOLEPREFIX_start(i - 1);
			Parareal_GenericData &tmp = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_fine();
			Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_fine_previous_timestep(); // SL

			// use coarse time step output data as initial data of next coarse time step
			CONSOLEPREFIX_start(i);
			parareal_simulationInstances[i]->sim_set_data(tmp);
			parareal_simulationInstances[i]->sim_set_data_fine_previous_time_slice(tmp2); // SL

			// run coarse time step
			parareal_simulationInstances[i]->run_timestep_fine();

			*(parareal_simulationInstances[i]->parareal_data_fine_exact_debug) = *(parareal_simulationInstances[i]->parareal_data_fine);
		}
#endif


		// Store initial solution:
		parareal_simulationInstances[0]->output_data_file(
				0,
				0,
				true
			);

		CONSOLEPREFIX_start("[MAIN] ");
		std::cout << "Initial propagation" << std::endl;

		/**
		 * Initial propagation
		 */
		CONSOLEPREFIX_start(0);
		parareal_simulationInstances[0]->run_timestep_coarse();
		*(parareal_simulationInstances[0]->parareal_data_output) = *(parareal_simulationInstances[0]->parareal_data_coarse);
		for (int i = 1; i < pVars->coarse_slices; i++)
		{
			CONSOLEPREFIX_start(i - 1);
			Parareal_GenericData &tmp = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_coarse();
			Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_coarse_previous_timestep(); // SL

			// use coarse time step output data as initial data of next coarse time step
			CONSOLEPREFIX_start(i);
			parareal_simulationInstances[i]->sim_set_data(tmp);
			parareal_simulationInstances[i]->sim_set_data_coarse_previous_time_slice(tmp2); // SL

			// run coarse time step
			parareal_simulationInstances[i]->run_timestep_coarse();

			*(parareal_simulationInstances[i]->parareal_data_output) = *(parareal_simulationInstances[i]->parareal_data_coarse);

		}

		// Store initial propagation:
		if (pVars->store_iterations)
			for (int i = 0; i < pVars->coarse_slices; i++)
				parareal_simulationInstances[i]->output_data_file(
						0,  // 0-th iteration
						i
					);
		// Store initial error relative to reference solution
		if (pVars->load_ref_csv_files)
			for (int i = 0; i < pVars->coarse_slices; i++)
				parareal_simulationInstances[i]->store_parareal_error(
						0,
						i,
						pVars->path_ref_csv_files,
						"ref");
		// Store initial error relative to fine solution
		if (pVars->load_fine_csv_files)
			for (int i = 0; i < pVars->coarse_slices; i++)
				parareal_simulationInstances[i]->store_parareal_error(
						0,
						i,
						pVars->path_fine_csv_files,
						"fine");


		/**
		 * We run as much Parareal iterations as there are coarse slices
		 */
//		int start_slice = 0;

		int k = 0;
		for (; k < pVars->coarse_slices; k++)
		{
			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Iteration Nr. " << k << std::endl;
			/*
			 * All the following loops should start with 0.
			 * For debugging reasons, we leave it here at 0
			 */


			/**
			 * Fine time stepping
			 */
			for (int i = k; i < pVars->coarse_slices; i++)
			{
				CONSOLEPREFIX_start(i);

				if (i > 0)
				{
				    // SL:
				    Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_fine_previous_timestep();
				    parareal_simulationInstances[i]->sim_set_data_fine_previous_time_slice(tmp2);
				}

				parareal_simulationInstances[i]->run_timestep_fine();
			}

			/**
			 * Compute difference between coarse and fine solution
			 */
			for (int i = k; i < pVars->coarse_slices; i++)
			{
				CONSOLEPREFIX_start(i);
				parareal_simulationInstances[i]->compute_difference();
			}

			/**
			 * 1) Coarse time stepping
			 * 2) Compute output + convergence check
			 * 3) Forward to next frame
			 */
			double max_convergence = -2;
			for (int i = k; i < pVars->coarse_slices; i++)
			{
				CONSOLEPREFIX_start(i);

				if (i > 0)
				{
				    // SL:
				    Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_coarse_previous_timestep();
				    parareal_simulationInstances[i]->sim_set_data_coarse_previous_time_slice(tmp2);
				}

				parareal_simulationInstances[i]->run_timestep_coarse();

				// compute convergence
				double convergence = parareal_simulationInstances[i]->compute_output_data(true);
				std::cout << "                        iteration " << k << ", time slice " << i << ", convergence: " << convergence << std::endl;
				if (max_convergence != -1)
					max_convergence = (convergence==-1)?(convergence):(std::max(max_convergence,convergence));

				if (pVars->store_iterations)
					parareal_simulationInstances[i]->output_data_file(
							k + 1,
							i
						);
				if (pVars->load_ref_csv_files)
					parareal_simulationInstances[i]->store_parareal_error(
							k + 1,
							i,
							pVars->path_ref_csv_files,
							"ref");
				if (pVars->load_fine_csv_files)
					parareal_simulationInstances[i]->store_parareal_error(
							k + 1,
							i,
							pVars->path_fine_csv_files,
							"fine");


				CONSOLEPREFIX.start(i);
				parareal_simulationInstances[i]->output_data_console(
						k + 1,
						i
					);

				// last coarse time step slice?
				if (i == pVars->coarse_slices-1)
				{
					// convergence check activated?
					if (pVars->convergence_error_threshold >= 0)
					{
						if (max_convergence >= 0)
						{
							// convergence given?
							if (max_convergence < pVars->convergence_error_threshold)
							{
								CONSOLEPREFIX_start("[MAIN] ");
								std::cout << "Convergence reached at iteration " << k << " with convergence value " << max_convergence << std::endl;
								goto converged;
							}
						}
					}
				}


				// forward to next time slice if it exists
				if (i < pVars->coarse_slices-1)
				{
					CONSOLEPREFIX_start(i);
					Parareal_GenericData &tmp = parareal_simulationInstances[i]->get_reference_to_output_data();

					CONSOLEPREFIX_start(i+1);
					parareal_simulationInstances[i+1]->sim_set_data(tmp);
				}

				parareal_simulationInstances[i]->check_for_nan_parareal();

#if SWEET_DEBUG
				// fine serial solution should be retrieved
				if (i == k)
				{
					std::cout << "Comparing parareal to fine solution at iteration " << i << " and end of timeframe " << i << ". Solutions should be identical." << std::endl;
					parareal_simulationInstances[i]->compare_to_fine_exact();
					std::cout << "Comparison OK!" << std::endl;
				}
#endif

			}
		}

converged:

		CONSOLEPREFIX_end();
	}
};





#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_ */
