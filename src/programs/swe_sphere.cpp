/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 * MULE_SCONS_OPTIONS: --fortran-source=enable --sphere-spectral-space=enable
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif


#include <stdexcept>

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData_Physical.hpp>
	#include <sweet/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
	#include <sweet/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>
#endif

#include "swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereOperators_SphereDataComplex.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>

#include <sweet/Stopwatch.hpp>
#include <sweet/SWEETError.hpp>

#include "swe_sphere_timeintegrators/SWE_Sphere_TimeSteppers.hpp"
#include "swe_sphere_timeintegrators/SWE_Sphere_NormalModeAnalysis.hpp"

#include <sweet/SimulationBenchmarkTiming.hpp>
#include <sweet/sphere/SphereData_DebugContainer.hpp>

#if SWEET_PARAREAL
	#include <parareal/Parareal.hpp>
#endif


SimulationVariables simVars;

// Plane data config
SphereData_Config sphereDataConfigInstance;
SphereData_Config sphereDataConfigInstance_nodealiasing;
SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;
SphereData_Config *sphereDataConfig_nodealiasing = &sphereDataConfigInstance_nodealiasing;


#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif



/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class SimulationInstance
{
public:
	SphereOperators_SphereData op;
	SphereOperators_SphereData op_nodealiasing;

	SWE_Sphere_TimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// Implementation of different time steppers
	SWE_Sphere_TimeSteppers timeSteppersCoarse;
#endif


	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	SphereData_Spectral prog_phi_pert;
	SphereData_Spectral prog_vrt;
	SphereData_Spectral prog_div;

	Stopwatch stopwatch;

#if SWEET_GUI
	PlaneData_Physical viz_plane_data;
#endif

	int render_primitive_id = 1;

	SphereHelpers_Diagnostics sphereDiagnostics;

#if SWEET_MPI
	int mpi_rank;
#endif

	// was the output of the time step already done for this simulation state?
	double timestep_last_output_simtime;

	BenchmarksSphereSWE sphereBenchmarks;

public:
	SimulationInstance()	:
		op(sphereDataConfig, &(simVars.sim)),
		op_nodealiasing(sphereDataConfig_nodealiasing, &(simVars.sim)),
		prog_phi_pert(sphereDataConfig),
		prog_vrt(sphereDataConfig),
		prog_div(sphereDataConfig),

#if SWEET_GUI
		viz_plane_data(planeDataConfig),
#endif
		sphereDiagnostics(
				sphereDataConfig,
				simVars,
				simVars.misc.verbosity
		)

	{
#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif


		reset();

	}



	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		sphereDiagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
				op,
				prog_phi_pert,
				prog_vrt,
				prog_div,
				simVars
		);
	}



	void reset()
	{
		SimulationBenchmarkTimings::getInstance().main_setup.start();

		simVars.reset();
		simVars.iodata.output_time_scale = 1.0/(60.0*60.0);

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		simVars.iodata.output_next_sim_seconds = 0;

		if (simVars.timecontrol.current_timestep_size <= 0)
			SWEETError("Only fixed time step size supported");

		if (simVars.benchmark.setup_dealiased)
		{
			//std::cout << "A" << std::endl;
			//std::cout << sphereDataConfig->getConfigInformationString() << std::endl;
			//exit(1);
			// use dealiased physical space for setup
			sphereBenchmarks.setup(simVars, op);
			sphereBenchmarks.master->get_initial_state(prog_phi_pert, prog_vrt, prog_div);
		}
		else
		{
			//std::cout << "B" << std::endl;
			//std::cout << sphereDataConfig_nodealiasing->getConfigInformationString() << std::endl;
			//exit(1);
			// this is not the default since noone uses it
			// use reduced physical space for setup to avoid spurious modes
			SphereData_Spectral prog_phi_pert_nodealiasing(sphereDataConfig_nodealiasing);
			SphereData_Spectral prog_vrt_nodealiasing(sphereDataConfig_nodealiasing);
			SphereData_Spectral prog_div_nodealiasing(sphereDataConfig_nodealiasing);

			sphereBenchmarks.setup(simVars, op_nodealiasing);
			sphereBenchmarks.master->get_initial_state(prog_phi_pert_nodealiasing, prog_vrt_nodealiasing, prog_div_nodealiasing);

			prog_phi_pert.load_nodealiasing(prog_phi_pert_nodealiasing);
			prog_vrt.load_nodealiasing(prog_vrt_nodealiasing);
			prog_div.load_nodealiasing(prog_div_nodealiasing);
		}

		/*
		 * SETUP time steppers
		 */
		timeSteppers.setup(simVars.disc.timestepping_method,
				op, simVars);

		std::cout << "[MULE] timestepper_string_id: " << timeSteppers.master->string_id() << std::endl;

		update_diagnostics();

		simVars.diag.backup_reference();

		SimulationBenchmarkTimings::getInstance().main_setup.stop();

		// start at one second in the past to ensure output at t=0
		timestep_last_output_simtime = simVars.timecontrol.current_simulation_time-1.0;

		/*
		 * Output configuration here to ensure that updated variables are included in this output
		 */
#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			simVars.outputConfig();
		}


		/*
		 * Output data for the first time step as well if output of datafiels is requested
		 */
		if (simVars.iodata.output_each_sim_seconds >= 0)
			timestep_do_output();

		stopwatch.start();
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv_spec_evol(
			const SphereData_Spectral &i_sphereData,
			const char* i_name		///< name of output variable
	)
	{
		char buffer[1024];
		std::string phase = "_arg";
		std::string ampl = "_amp"; 

		const char* filename_template_ampl = "output_spec_ampl_%s.txt"; //.c_str();
		const char* filename_template_arg = "output_spec_arg_%s.txt"; //.c_str();
		int reduce_mode_factor = 4;

		sprintf(buffer, filename_template_arg, i_name);
		i_sphereData.spectrum_phase_file_write_line(buffer, 
			i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale,
			20, 10e-20, reduce_mode_factor);

		sprintf(buffer, filename_template_ampl, i_name);
		i_sphereData.spectrum_abs_file_write_line(buffer, 
			i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale,
			20, 10e-20, reduce_mode_factor);

		return buffer;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const SphereData_Spectral &i_sphereData,
			const char* i_name,		///< name of output variable
			bool i_phi_shifted = false
	)
	{
		char buffer[1024];

		// create copy
		SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin(
			const SphereData_Spectral &i_sphereData,
			const char* i_name
	)
	{
		char buffer[1024];

		SphereData_Spectral sphereData(i_sphereData);
		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return;
#endif

		if (simVars.iodata.output_file_name.length() == 0)
			return;


		std::cout << "Writing output files at simulation time: " << simVars.timecontrol.current_simulation_time << " secs" << std::endl;

		if (simVars.iodata.output_file_mode == "csv")
		{
			std::string output_filename;

			SphereData_Spectral h = prog_phi_pert*(1.0/simVars.sim.gravitation);
			h += simVars.sim.h0;

			output_filename = write_file_csv(h, "prog_h");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;

			output_filename = write_file_csv(prog_phi_pert, "prog_phi_pert");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << prog_phi_pert.toPhys().physical_reduce_min() << ", max: " << prog_phi_pert.toPhys().physical_reduce_max() << ")" << std::endl;

			SphereData_Physical u(sphereDataConfig);
			SphereData_Physical v(sphereDataConfig);

			op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);

			output_filename = write_file_csv(u, "prog_u");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(v, "prog_v");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(prog_vrt, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(prog_div, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			SphereData_Spectral potvrt = (prog_phi_pert/simVars.sim.gravitation)*prog_vrt;

			output_filename = write_file_csv(potvrt, "prog_potvrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;
		}
		else if (simVars.iodata.output_file_mode == "bin")
		{
			std::string output_filename;

			{
				output_filename = write_file_bin(prog_phi_pert, "prog_phi_pert");
				output_reference_filenames = output_filename;
				SphereData_Physical prog_phys = prog_phi_pert.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(prog_vrt, "prog_vrt");
				output_reference_filenames += ";"+output_filename;
				SphereData_Physical prog_phys = prog_vrt.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(prog_div, "prog_div");
				output_reference_filenames += ";"+output_filename;
				SphereData_Physical prog_phys = prog_div.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}
		}
		else if (simVars.iodata.output_file_mode == "csv_spec_evol"){

			std::string output_filename;

			{ 
				/*
				* Spectral kinetic energy and potential enstrophy calculation and output
				*
				* Details in Jakob-Chien, Ruediger, James J. Hack, and David L. Williamson. 
				* "Spectral transform solutions to the shallow water test set." Journal of Computational Physics 119, no. 1 (1995): 164-187.
				*/
				// Kinetic energy is given in spectral space as
				// KE per mode = a^2/((n(n+1)))*(vrt*conj(vrt))+a^2/((n(n+1)))*(div*conj(div))
				// r = a/(sqrt(n(n+1))) (root_laplace)
				// KE per mode = (r*vrt*conj(r*vrt))+(r*div*conj(r*div))
				SphereData_Spectral rlap_vrt = op.inv_root_laplace(prog_vrt); 
				SphereData_Spectral rlap_div = op.inv_root_laplace(prog_div); 
				SphereData_Spectral kin_en = rlap_vrt + rlap_div ;

				output_filename = write_file_csv_spec_evol(kin_en, "kin_en"); 
				std::cout << " + " << output_filename << " (Total Kin Energy : " << 0.25*kin_en.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// For Barotropic vort eq: See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// Kinetic energy is given in spectral space as
				// Vortical energy per mode is (0.5 n*(n+1) / a^2) *psi*conj(psi) in spectral space
				//SphereData_Spectral psi = op.inv_laplace(prog_vrt); // 
				// multiply psi by sqrt( n * (n+1))/a (apply root laplacian)
				//SphereData_Spectral psi_root = op.root_laplace(psi);
				//output_filename = write_file_csv_spec_evol(psi_root*std::sqrt(0.5), "spec_energy"); 
				//std::cout << " + " << output_filename << " (Kinetic energy : " << (0.5)*psi_root.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// enstrophy per mode is 0.5 vrt*conj(vrt) in spectral space
				// Total enstrophy is the sum of these (counting twice modes with m>0 and once when m=0)
				output_filename = write_file_csv_spec_evol(prog_vrt, "enstrophy"); 
				std::cout << " + " << output_filename << " (Total Enstrophy : " << prog_vrt.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

			}
		}
		else
		{
			SWEETError("Unknown output file mode '"+simVars.iodata.output_file_mode+"'");
		}
	}



	void timestep_do_output()
	{
		if (simVars.misc.compute_errors)
		{
			/*
			 * Check for stationary solutions
			 */
			if (
					simVars.benchmark.benchmark_name != "williamson2"					&&
					simVars.benchmark.benchmark_name != "williamson2_linear"			&&
					simVars.benchmark.benchmark_name != "galewsky_nobump"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_linear"	&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_1"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_2"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_4"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_8"			&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_16"		&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_32"		&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_64"		&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_128"		&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_256"		&&
					simVars.benchmark.benchmark_name != "geostrophic_balance_512"
			)
			{
				std::cout << "Benchmark name: " << simVars.benchmark.benchmark_name << std::endl;
				SWEETError("Analytical solution not available for this benchmark");
			}

			SphereData_Spectral anal_solution_phi_pert(sphereDataConfig);
			SphereData_Spectral anal_solution_vrt(sphereDataConfig);
			SphereData_Spectral anal_solution_div(sphereDataConfig);

			sphereBenchmarks.setup(simVars, op);
			sphereBenchmarks.master->get_initial_state(anal_solution_phi_pert, anal_solution_vrt, anal_solution_div);

			/*
			 * Compute difference
			 */
			SphereData_Spectral diff_phi = prog_phi_pert - anal_solution_phi_pert;
			SphereData_Spectral diff_vrt = prog_vrt - anal_solution_vrt;
			SphereData_Spectral diff_div = prog_div - anal_solution_div;

#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				double error_phi = diff_phi.toPhys().physical_reduce_max_abs();
				double error_vrt = diff_vrt.toPhys().physical_reduce_max_abs();
				double error_div = diff_div.toPhys().physical_reduce_max_abs();

				
				std::ios init(NULL);
				init.copyfmt(std::cout);
				std::cout << "[MULE] errors." << std::setw(8) << std::setfill('0') << simVars.timecontrol.current_timestep_nr << ": ";
				std::cout.copyfmt(init);

				std::cout << "simtime=" << simVars.timecontrol.current_simulation_time;
				std::cout << "\terror_linf_phi=" << error_phi;
				std::cout << "\terror_linf_vrt=" << error_vrt;
				std::cout << "\terror_linf_div=" << error_div;
				std::cout << std::endl;
			}
		}

		write_file_output();

		update_diagnostics();

		if (simVars.misc.verbosity > 1)
		{

#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				update_diagnostics();

				// Print header
				if (simVars.timecontrol.current_timestep_nr == 0)
				{
					std::cout << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
					std::cout << std::endl;
				}

				// Print simulation time, energy and pot enstrophy
				std::cout << simVars.timecontrol.current_simulation_time << "\t";
				std::cout << simVars.diag.total_mass << "\t";
				std::cout << simVars.diag.potential_energy << "\t";
				std::cout << simVars.diag.kinetic_energy << "\t";
				std::cout << simVars.diag.total_energy << "\t";
				std::cout << simVars.diag.total_potential_enstrophy << "\t";

				std::cout << (simVars.diag.total_mass-simVars.diag.ref_total_mass)/simVars.diag.total_mass << "\t";
				std::cout << (simVars.diag.potential_energy-simVars.diag.ref_potential_energy)/simVars.diag.potential_energy << "\t";
				std::cout << (simVars.diag.kinetic_energy-simVars.diag.ref_kinetic_energy)/simVars.diag.kinetic_energy << "\t";
				std::cout << (simVars.diag.total_energy-simVars.diag.total_energy)/simVars.diag.total_energy << "\t";
				std::cout << (simVars.diag.total_potential_enstrophy-simVars.diag.total_potential_enstrophy)/simVars.diag.total_potential_enstrophy << std::endl;

				static double start_tot_energy = -1;
				if (start_tot_energy == -1)
					start_tot_energy = simVars.diag.total_energy;
			}
		}


		if (simVars.misc.verbosity > 0)
		{
#if SWEET_MPI
			if (mpi_rank == 0)
#endif
				std::cout << "prog_phi min/max:\t" << prog_phi_pert.toPhys().physical_reduce_min() << ", " << prog_phi_pert.toPhys().physical_reduce_max() << std::endl;
		}

		if (simVars.iodata.output_each_sim_seconds > 0)
			while (simVars.iodata.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.iodata.output_next_sim_seconds += simVars.iodata.output_each_sim_seconds;
	}




public:
	bool timestep_check_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return false;
#endif
		if (simVars.misc.gui_enabled)
			update_diagnostics();

		if (simVars.misc.verbosity > 0)
			std::cout << "." << std::flush;

		if (simVars.iodata.output_each_sim_seconds < 0)
			return false;

		if (simVars.timecontrol.current_simulation_time == timestep_last_output_simtime)
			return false;

		timestep_last_output_simtime = simVars.timecontrol.current_simulation_time;

		if (simVars.timecontrol.current_simulation_time < simVars.timecontrol.max_simulation_time - simVars.iodata.output_each_sim_seconds*1e-10)
		{
			if (simVars.iodata.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
				return false;
		}

		if (simVars.misc.verbosity > 0)
			std::cout << std::endl;

		timestep_do_output();

		return true;
	}



public:
	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		if (simVars.timecontrol.max_wallclock_time >= 0)
		{
			double t = stopwatch.getTimeSinceStart();
			if (simVars.timecontrol.max_wallclock_time <= t)
			{
				std::cout << "[MULE] max_wallclock_time: exceeded (" << t << ")" << std::endl;
				return true;
			}
		}

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



	bool detect_instability()
	{
		if (prog_phi_pert.spectral_is_first_nan_or_inf())
		{
			std::cout << "Infinity value detected" << std::endl;
			std::cerr << "Infinity value detected" << std::endl;
			return true;
		}

		return false;
	}



	void run_timestep()
	{
#if SWEET_GUI
		if (simVars.misc.gui_enabled && simVars.misc.normal_mode_analysis_generation == 0)
			timestep_check_output();
#endif

		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_phi_pert, prog_vrt, prog_div,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);


		// advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_check_output();
#endif
	}



	void normalmode_analysis()
	{
		NormalModeAnalysisSphere::normal_mode_analysis(
				prog_phi_pert,
				prog_vrt,
				prog_div,
				simVars,
				this,
				&SimulationInstance::run_timestep
			);
	}



#if SWEET_GUI

	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				run_timestep();
	}


	int max_viz_types = 9;


	void vis_get_vis_data_array(
			const PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		// request rendering of sphere
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		if (simVars.misc.vis_id < 0)
		{
			int n = -simVars.misc.vis_id-1;
			if (n <  (int)SphereData_DebugContainer().size())
			{
				SphereData_DebugContainer::DataContainer &d = SphereData_DebugContainer().container_data()[n];
				if (d.is_spectral)
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(d.data_spectral, planeDataConfig);
				else
					viz_plane_data = Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(d.data_physical, planeDataConfig);

				*o_dataArray = &viz_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
		}

		int id = simVars.misc.vis_id % max_viz_types;

		switch (id)
		{
			default:

			case 0:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_phi_pert), planeDataConfig);
				break;

			case 1:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_vrt), planeDataConfig);
				break;

			case 2:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_div), planeDataConfig);
				break;

			case 3:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(simVars.sim.h0 + SphereData_Spectral(prog_phi_pert)/simVars.sim.gravitation, planeDataConfig);
				break;

			case 4:
			{
				SphereData_Physical u(prog_vrt.sphereDataConfig);
				SphereData_Physical v(prog_vrt.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(u, planeDataConfig);
				break;
			}

			case 5:
			{
				SphereData_Physical u(prog_vrt.sphereDataConfig);
				SphereData_Physical v(prog_vrt.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(v, planeDataConfig);
				break;
			}

			case 6:
			case 7:
			case 8:
			{
				SphereData_Spectral anal_solution_phi_pert(sphereDataConfig);
				SphereData_Spectral anal_solution_vrt(sphereDataConfig);
				SphereData_Spectral anal_solution_div(sphereDataConfig);

				sphereBenchmarks.setup(simVars, op);
				sphereBenchmarks.master->get_initial_state(anal_solution_phi_pert, anal_solution_vrt, anal_solution_div);

				switch (id)
				{
				case 6:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_phi_pert - anal_solution_phi_pert, planeDataConfig);
					break;

				case 7:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_vrt - anal_solution_vrt, planeDataConfig);
					break;

				case 8:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_div - anal_solution_div, planeDataConfig);
					break;
				}
			}
		}


		double viz_min = viz_plane_data.physical_reduce_min();
		double viz_max = viz_plane_data.physical_reduce_max();

		viz_max = std::max(std::abs(viz_max), std::abs(viz_min));
		viz_min = -viz_max;

		*o_viz_min = viz_min;
		*o_viz_max = viz_max;


		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		std::string description = "";


		bool found = false;
		if (simVars.misc.vis_id < 0)
		{
			int n = -simVars.misc.vis_id-1;

			if (n <  (int)SphereData_DebugContainer().size())
			{
				description = std::string("DEBUG_")+SphereData_DebugContainer().container_data()[n].description;
				found = true;
			}
		}

		int id = simVars.misc.vis_id % max_viz_types;

		if (!found)
		{
			switch (id)
			{
			default:
			case 0:
				description = "phi_pert";
				break;

			case 1:
				description = "vrt";
				break;

			case 2:
				description = "div";
				break;

			case 3:
				description = "h";
				break;

			case 4:
				description = "u";
				break;

			case 5:
				description = "v";
				break;

			case 6:
				description = "phi diff t0";
				break;

			case 7:
				description = "vrt diff t0";
				break;

			case 8:
				description = "div diff t0";
				break;
			}
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
#if SWEET_MPI
				"Rank %i,"
				","
#endif
				"Visualization %i: %s,"
				"MaxVal: %.12e,"
				"MinVal: %.12e,"
				","
				"Time: %f secs,"
				"Time: %f hours,"
				"Time: %f days,"
				"timestep nr.: %i,"
				"timestep size: %f,"
				","
				"TMass: %.12e,"
				"TEnergy: %.12e,"
				"PotEnstrophy: %.12e,"
				","
				"Colorscale: lowest [Blue... green ... red] highest"
				,
#if SWEET_MPI
				mpi_rank,
#endif
				simVars.misc.vis_id,
				description.c_str(),
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min(),

				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0),
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,

				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy

		);

		return title_string;
	}



	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
			break;

		case 'c':
			write_file_output();
			break;
		}
	}
#endif


};


int main_real(int i_argc, char *i_argv[])
{

	// Time counter
	SimulationBenchmarkTimings::getInstance().main.start();

#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING_SPACE
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
		{
				std::cerr << "MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye..." << std::endl;
				exit(-1);
		}
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	// input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			nullptr
	};

#if SWEET_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.printOptions();
#endif
		return -1;
	}

	if (simVars.misc.verbosity > 3)
		std::cout << " + setup SH sphere transformations..." << std::endl;

	sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans, simVars.misc.verbosity);

	int res_physical_nodealias[2] = {
			2*simVars.disc.space_res_spectral[0],
			simVars.disc.space_res_spectral[1]
		};

	if (simVars.misc.verbosity > 3)
		std::cout << " + setup SH sphere transformations (nodealiasing)..." << std::endl;

	sphereDataConfigInstance_nodealiasing.setupAuto(res_physical_nodealias, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans, simVars.misc.verbosity);


#if SWEET_GUI
	if (simVars.misc.verbosity > 3)
		std::cout << " + setup FFT plane transformations..." << std::endl;

	planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);

	if (simVars.misc.verbosity > 3)
		std::cout << " + setup finished" << std::endl;

#if SWEET_MPI
	std::cout << "Helo from MPI rank: " << mpi_rank << std::endl;

	// only start simulation and time stepping for first rank
	if (mpi_rank > 0)
	{
		/*
		 * Deactivate all output for ranks larger than the current one
		 */
		simVars.misc.verbosity = 0;
		simVars.iodata.output_each_sim_seconds = -1;
	}
#endif

	{
#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			std::cout << "SPH config string: " << sphereDataConfigInstance.getConfigInformationString() << std::endl;
		}

#if SWEET_PARAREAL
		if (simVars.parareal.enabled)
		{

			simVars.iodata.output_time_scale = 1.0/(60.0*60.0);

			//SphereOperators op(sphereDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);
			SphereOperators_SphereData op(sphereDataConfig, &(simVars.sim));
			SphereOperators_SphereData op_nodealiasing(sphereDataConfig_nodealiasing, &(simVars.sim));

			SWE_Sphere_TimeSteppers* timeSteppersFine = new SWE_Sphere_TimeSteppers;
			SWE_Sphere_TimeSteppers* timeSteppersCoarse = new SWE_Sphere_TimeSteppers;

			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller<SWE_Sphere_TimeSteppers, 3> parareal_Controller(&simVars,
												sphereDataConfig,
												op,
												op_nodealiasing,
												timeSteppersFine,
												timeSteppersCoarse);

			// setup controller. This initializes several simulation instances
			parareal_Controller.setup();

			// execute the simulation
			parareal_Controller.run();

			delete timeSteppersFine;
			delete timeSteppersCoarse;
		}
		else
#endif

#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationSWE);
			delete simulationSWE;
		}
		else
#endif
		{
			SimulationInstance *simulationSWE = new SimulationInstance;

			if (simVars.misc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normalmode_analysis();
			}
			else
			{
				// Do first output before starting timer
				simulationSWE->timestep_check_output();
#if SWEET_MPI
				// Start counting time
				if (mpi_rank == 0)
				{
					std::cout << "********************************************************************************" << std::endl;
					std::cout << "Parallel performance information: MPI barrier & timer starts here" << std::endl;
					std::cout << "********************************************************************************" << std::endl;
				}
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				SimulationBenchmarkTimings::getInstance().main_timestepping.start();

				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Test for some output to be done
					simulationSWE->timestep_check_output();

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simVars.misc.instability_checks)
					{
#if SWEET_MPI
						if (mpi_rank == 0)
#endif
						{
							if (simulationSWE->detect_instability())
							{
								std::cout << "INSTABILITY DETECTED" << std::endl;
								std::cerr << "INSTABILITY DETECTED" << std::endl;
								// IMPORANT: EXIT IN CASE OF INSTABILITIES
								exit(1);
								break;
							}
						}
					}
				}

				// Stop counting time
				SimulationBenchmarkTimings::getInstance().main_timestepping.stop();

#if SWEET_MPI
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				if (simVars.misc.verbosity > 0)
					std::cout << std::endl;
#if SWEET_MPI
				// Start counting time
				if (mpi_rank == 0)
				{
					std::cout << "********************************************************************************" << std::endl;
					std::cout << "Parallel performance information: timer stopped here" << std::endl;
					std::cout << "********************************************************************************" << std::endl;
				}
#endif

				// Do some output after the time loop
				simulationSWE->timestep_check_output();
			}

#if SWEET_MPI
			// Start counting time
			if (mpi_rank == 0)
#endif
			{
				if (simVars.iodata.output_file_name.size() > 0)
					std::cout << "[MULE] reference_filenames: " << simulationSWE->output_reference_filenames << std::endl;
			}

			std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;

			delete simulationSWE;
		}

		SimulationBenchmarkTimings::getInstance().main.stop();
	}


#if SWEET_MPI
	if (mpi_rank == 0)
#endif
	{
		// End of run output results
		std::cout << std::endl;
		SimulationBenchmarkTimings::getInstance().output();

		std::cout << "***************************************************" << std::endl;
		std::cout << "* Other timing information (direct)" << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "[MULE] simVars.timecontrol.current_timestep_nr: " << simVars.timecontrol.current_timestep_nr << std::endl;
		std::cout << "[MULE] simVars.timecontrol.current_timestep_size: " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "* Other timing information (derived)" << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "[MULE] simulation_benchmark_timings.time_per_time_step (secs/ts): " << SimulationBenchmarkTimings::getInstance().main_timestepping()/(double)simVars.timecontrol.current_timestep_nr << std::endl;
	}

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}




int main(int i_argc, char *i_argv[])
{
	try
	{
		return main_real(i_argc, i_argv);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}

}

