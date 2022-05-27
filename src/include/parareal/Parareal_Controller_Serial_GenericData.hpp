/*
 * PararealController.hpp
 *
 *  Created on: 11 Apr 2016
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_
#define SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_

// Checking if geometry and model have been correctly defined
#if SWEET_PARAREAL_SCALAR
	#if SWEET_PARAREAL_PLANE || SWEET_PARAREAL_SPHERE
		#error "More than one geometry has been defined for parareal"
	#endif
#elif SWEET_PARAREAL_PLANE
	#if SWEET_PARAREAL_SCALAR || SWEET_PARAREAL_SPHERE
		#error "More than one geometry has been defined for parareal"
	#endif
	#if (!SWEET_PARAREAL_PLANE_SWE) && (!SWEET_PARAREAL_PLANE_BURGERS)
		#error "No model has been defined for parareal on the plane"
	#endif
	#if SWEET_PARAREAL_PLANE_SWE && SWEET_PARAREAL_PLANE_BURGERS
		#error "More than one model has been defined for parareal on the plane"
	#endif
#elif SWEET_PARAREAL_SPHERE
	#if SWEET_PARAREAL_SCALAR || SWEET_PARAREAL_PLANE
		#error "More than one geometry has been defined for parareal"
	#endif
#else
	#error "No geometry has been defined for parareal"
#endif


#include <parareal/Parareal_ConsolePrefix.hpp>
#include <parareal/Parareal_SimulationInstance_GenericData.hpp>
#include <parareal/Parareal_SimulationVariables.hpp>

#include <parareal/Parareal_GenericData.hpp>

#if SWEET_PARAREAL_SCALAR
#include <parareal/Parareal_GenericData_Scalar.hpp>

#elif SWEET_PARAREAL_PLANE
#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#elif SWEET_PARAREAL_SPHERE
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <map>

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
	 * Pointers to interfaces of simulationInstances
	 * This helps to clearly separate between the allocation of the simulation classes and the parareal interfaces.
	 */
	std::vector<Parareal_SimulationInstance_GenericData<t_tsmType, N>*> parareal_simulationInstances = {};

	SimulationVariables* simVars;


	// Operators and DataConfig
#if SWEET_PARAREAL_PLANE
	PlaneOperators op_plane;
	PlaneDataConfig* planeDataConfig = nullptr;
#elif SWEET_PARAREAL_SPHERE
	SphereOperators_SphereData op_sphere;
	SphereOperators_SphereData op_sphere_nodealiasing;
	SphereData_Config* sphereDataConfig = nullptr;
#endif


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


// MPI
	int mpi_nprocs = 1;
	int mpi_rank = 0;
	std::vector<int> slices_for_proc = {};
	std::vector<int> proc_for_slices = {};
	std::map<int, int> global_to_local_slice;



public:

#if SWEET_PARAREAL==2

	// every communication is made between mpi_proc = 0 and mpi_proc > 0
	void communicate_solution( Parareal_GenericData* io_data, int source_rank, int dest_rank)
	{

		bool send = false;
		bool recv = false;

		////if (end_timeslice == slices_per_proc.back())
		////	send = true;
		////else if (end_timeslice == slices_per_proc[0] - 1)
		////	recv = true;
		if (mpi_rank == source_rank)
			send = true;
		else if (mpi_rank == dest_rank)
			recv = true;
		else
			return;

	#if SWEET_PARAREAL_SCALAR
		double u;
	#elif SWEET_PARAREAL_PLANE_SWE
		PlaneData_Spectral h(planeDataConfig);
		PlaneData_Spectral u(planeDataConfig);
		PlaneData_Spectral v(planeDataConfig);
	#elif SWEET_PARAREAL_PLANE_BURGERS
		PlaneData_Spectral u(planeDataConfig);
		PlaneData_Spectral v(planeDataConfig);
	#elif SWEET_PARAREAL_SPHERE
		SphereData_Spectral phi(sphereDataConfig);
		SphereData_Spectral vrt(sphereDataConfig);
		SphereData_Spectral div(sphereDataConfig);
	#endif

		int size = io_data.size();
		int tag = 1000 * end_timeslice;

		if (send)
		{

	#if SWEET_PARAREAL_SCALAR
			parareal_simulationInstances[0]->GenericData_Scalar_to_dataArrays(io_data, u);
			MPI_Send(u, size, MPI_DOUBLE, dest_rank, tag + 0, MPI_COMM_WORLD);
	#elif SWEET_PARAREAL_PLANE_SWE
			parareal_simulationInstances[0]->GenericData_PlaneData_Spectral_to_dataArrays(io_data, h, u, v);
			MPI_Send(h.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Send(u.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 1, MPI_COMM_WORLD);
			MPI_Send(v.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 2, MPI_COMM_WORLD);
	#elif SWEET_PARAREAL_PLANE_BURGERS
			parareal_simulationInstances[0]->GenericData_PlaneData_Spectral_to_dataArrays(io_data, u, v);
			MPI_Send(u.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Send(v.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 1, MPI_COMM_WORLD);
	#elif SWEET_PARAREAL_SPHERE
			parareal_simulationInstances[0]->GenericData_SphereData_Spectral_to_dataArrays(io_data, phi, vrt, div);
			MPI_Send(phi.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Send(vrt.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 1, MPI_COMM_WORLD);
			MPI_Send(div.spectral_space_data, size, MPI_DOUBLE, dest_rank, tag + 2, MPI_COMM_WORLD);
	#endif

		}
		else if (recv)
		{

	#if SWEET_PARAREAL_SCALAR
			MPI_Recv(u, size, MPI_DOUBLE, source_rank, tag + 0, MPI_COMM_WORLD);
			parareal_simulationInstances[0]->GenericData_Scalar_to_dataArrays(io_data, u);
	#elif SWEET_PARAREAL_PLANE_SWE
			MPI_Recv(h.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Recv(u.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 1, MPI_COMM_WORLD);
			MPI_Recv(v.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 2, MPI_COMM_WORLD);
			parareal_simulationInstances[0]->data_arrays_to_GenericData_PlaneData_Spectral(io_data, h, u, v);
	#elif SWEET_PARAREAL_PLANE_BURGERS
			MPI_Recv(u.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Recv(v.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 1, MPI_COMM_WORLD);
			parareal_simulationInstances[0]->dataArrays_to_GenericData_PlaneData_Spectral(io_data, u, v);
	#elif SWEET_PARAREAL_SPHERE
			MPI_Recv(phi.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 0, MPI_COMM_WORLD);
			MPI_Recv(vrt.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 1, MPI_COMM_WORLD);
			MPI_Recv(div.spectral_space_data, size, MPI_DOUBLE, source_rank, tag + 2, MPI_COMM_WORLD);
			parareal_simulationInstances[0]->dataArrays_to_GenericData_SphereData_Spectral(io_data, phi, vrt, div);
	#endif


		}

	}


#endif


#if SWEET_PARAREAL_SCALAR
	// Scalar
	Parareal_Controller_Serial_GenericData(SimulationVariables* i_simVars,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

#elif SWEET_PARAREAL_PLANE
	// Plane
	Parareal_Controller_Serial_GenericData(SimulationVariables* i_simVars,
						PlaneDataConfig* i_planeDataConfig,
						PlaneOperators &i_op_plane,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		planeDataConfig(i_planeDataConfig),
		op_plane(i_op_plane),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};

#elif SWEET_PARAREAL_SPHERE
	// Sphere
	Parareal_Controller_Serial_GenericData(SimulationVariables* i_simVars,
						SphereData_Config* i_sphereDataConfig,
						SphereOperators_SphereData &i_op_sphere,
						SphereOperators_SphereData &i_op_sphere_nodealiasing,
						t_tsmType* i_timeSteppersFine,
						t_tsmType* i_timeSteppersCoarse):
		simVars(i_simVars),
		sphereDataConfig(i_sphereDataConfig),
		op_sphere(i_op_sphere),
		op_sphere_nodealiasing(i_op_sphere_nodealiasing),
		timeSteppersFine(i_timeSteppersFine),
		timeSteppersCoarse(i_timeSteppersCoarse)
	{
	};
#endif

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
			//Parareal_SimulationVariables *i_pararealSimVars
	)
	{

		cleanup();

#if SWEET_PARAREAL==2
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
#endif


		if (mpi_rank == 0)
		{
			pVars = &this->simVars->parareal;

			if (!pVars->enabled)
				return;

			if (pVars->coarse_slices <= 0)
			{
				std::cerr << "Invalid number of coarse slices" << std::endl;
				exit(1);
			}

			if (pVars->coarse_slices % mpi_nprocs != 0)
			{
				SWEETError("Number of coarse slices must be a multiple integer of the number of MPI processes");
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

		}


		// Distribute slices for each MPI proc
		// mpi_rank = 0 contains all processes to compute all serial parts
		int slices_per_proc = pVars->coarse_slices / mpi_nprocs;

		if (mpi_rank == 0)
			for (int k = 0; k < pVars->coarse_slices; k++)
				global_to_local_slice.emplace(std::make_pair(k, k)); // slice k is the i-th slice treated by proc mpi_rank

		int i = 0;
		for (int k = mpi_rank * slices_per_proc; k < (mpi_rank + 1) * slices_per_proc; k++)
		{
			slices_for_proc.push_back(k); // slice k is treated by proc mpi_rank
			if (mpi_rank > 0)
				global_to_local_slice.emplace(std::make_pair(k, i)); // slice k is the i-th slice treated by proc mpi_rank
			i++;
		}
		if (mpi_rank > 0)
			for (int k = 0; k < pVars->coarse_slices; k++)
				if (k < mpi_rank * slices_per_proc || k >= (mpi_rank + 1) * slices_per_proc)
					global_to_local_slice.emplace(std::make_pair(k, -1)); // slice k is not treated by proc mpi_rank



		//// Distribute MPI proc for each slice
		proc_for_slices = std::vector<int>(pVars->coarse_slices);
		if (mpi_rank == 0)
			for (int i = 0; i < mpi_nprocs; i++)
				for (int k = i * slices_per_proc; k < (i + 1) * slices_per_proc; k++)
					proc_for_slices[k] = i;

#if SWEET_PARAREAL==2
		MPI_Bcast(&proc_for_slices[0], proc_for_slices.size() , MPI_INT, 0, MPI_COMM_WORLD );
#endif


		// convert to pararealsimulationInstances to get Parareal interfaces
		for (int k = 0; k < pVars->coarse_slices; k++)
		{

			// except for mpi_rank = 0, only create instances for slices treated by mpi_rank
			if (mpi_rank > 0)
				if (k < slices_for_proc[0] || k > slices_for_proc.back())
					continue;

			CONSOLEPREFIX_start(k);

			int local_k = global_to_local_slice.at(k);
			parareal_simulationInstances.push_back(new Parareal_SimulationInstance_GenericData<t_tsmType, N>);
#if SWEET_PARAREAL_SCALAR
				parareal_simulationInstances[local_k]->setup(this->simVars,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);

#elif SWEET_PARAREAL_PLANE
				parareal_simulationInstances[local_k]->setup(this->simVars,
								       this->planeDataConfig,
								       &this->op_plane,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);

#elif SWEET_PARAREAL_SPHERE
				parareal_simulationInstances[local_k]->setup(this->simVars,
								       this->sphereDataConfig,
								       &this->op_sphere,
								       &this->op_sphere_nodealiasing,
								       this->timeSteppersFine,
								       this->timeSteppersCoarse);
#endif

		}


		if (mpi_rank == 0)
		{
			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Setup time frames" << std::endl;
		}

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

		////for (int k = 0; k < pVars->coarse_slices; k++)
		for (int k = slices_for_proc[0]; k <= slices_for_proc.back(); k++)
		{
			CONSOLEPREFIX_start(k);
			int local_k = global_to_local_slice.at(k);
			parareal_simulationInstances[local_k]->sim_set_timeframe(time_slice_size*k, time_slice_size*(k+1));
		}

		/*
		 * Setup first simulation instance
		 */
		if (mpi_rank == 0)
		{
			CONSOLEPREFIX_start(0);
			parareal_simulationInstances[0]->sim_setup_initial_data();
			CONSOLEPREFIX_end();

			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Finished setup parareal" << std::endl;
			CONSOLEPREFIX_end();
		}
	}



	void run()
	{

#if SWEET_DEBUG
		// DEBUG: full fine simulation
		if (mpi_rank == 0)
		{
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

		}
#endif


		if (mpi_rank == 0)
		{
			// Store initial solution:
			parareal_simulationInstances[0]->output_data_file(
					0,
					0,
					true
				);

			CONSOLEPREFIX_start("[MAIN] ");
			std::cout << "Initial propagation" << std::endl;
		}

		/**
		 * Initial propagation
		 */
		if (mpi_rank == 0)
		{
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
		}

		/**
		 * We run as much Parareal iterations as there are coarse slices
		 */
//		int start_slice = 0;

		int k = 0;
		for (; k < pVars->coarse_slices; k++)
		{


			if (mpi_rank == 0)
			{
				CONSOLEPREFIX_start("[MAIN] ");
				std::cout << "Iteration Nr. " << k << std::endl;
			}
			/*
			 * All the following loops should start with 0.
			 * For debugging reasons, we leave it here at 0
			 */


			/**
			 * Fine time stepping (in parallel)
			 */
			//for (int i = k; i < pVars->coarse_slices; i++)
			for (int i = slices_for_proc[0]; i <= slices_for_proc.back(); i++)
			{

				// solution already converged at time slices i < k
				if (i < k)
					continue;

				int working_rank = proc_for_slices[i];
				int local_slice = global_to_local_slice.at(i);

				if (mpi_rank == working_rank)
					CONSOLEPREFIX_start(i);

				if (mpi_rank == working_rank || mpi_rank == 0)
				{

					// For SL
					if (i > 0)
					{
						int local_prev_slice = global_to_local_slice.at(i-1);

						Parareal_GenericData* tmp2;
						// data taken from mpi_rank itself
						if (local_prev_slice >= 0)
							tmp2 = &parareal_simulationInstances[local_prev_slice]->get_reference_to_data_timestep_fine_previous_timestep();
							//Parareal_GenericData &tmp2 = parareal_simulationInstances[i-1]->get_reference_to_data_timestep_fine_previous_timestep();
#if SWEET_PARAREAL==2
						else
							this->communicate_solution(&tmp2, 0, mpi_rank);
#endif

						if (mpi_rank == working_rank)
							parareal_simulationInstances[local_slice]->sim_set_data_fine_previous_time_slice(*tmp2);
					}

					if (mpi_rank == working_rank)
						parareal_simulationInstances[local_slice]->run_timestep_fine();
				}
			}

			/**
			 * Compute difference between coarse and fine solution
			 */
			///for (int i = k; i < pVars->coarse_slices; i++)
			for (int i = slices_for_proc[0]; i <= slices_for_proc.back(); i++)
			{
				CONSOLEPREFIX_start(i);
				int local_slice = global_to_local_slice.at(i);
				parareal_simulationInstances[local_slice]->compute_difference();
#if SWEET_PARAREAL==2
				// communicate from each proc to mpi_rank = 0
				Parareal_GenericData* tmp2;
				int working_rank = proc_for_slices[i];
				if (mpi_rank == working_rank || mpi_rank == 0)
				{
					tmp2 = parareal_simulationInstances[local_slice]->get_reference_to_data_diff();
					this->communicate_solution(tmp2, mpi_rank, 0);
					if (mpi_rank == 0)
						parareal_simulationInstances[local_slice]->sim_set_data_diff(tmp2);
				}
#endif
			}

			/**
			 * 1) Coarse time stepping (serial)
			 * 2) Compute output + convergence check (serial)
			 * 3) Forward to next frame (serial)
			 */
			if (mpi_rank == 0)
			{
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
		}

converged:

		CONSOLEPREFIX_end();
	}
};





#endif /* SRC_INCLUDE_PARAREAL_PARAREAL_CONTROLLER_SERIAL_GENERICDATA_HPP_ */
