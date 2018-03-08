/*
 * swe_sphere.cpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData.hpp>
	#include <sweet/Convert_SphereData_To_PlaneData.hpp>
#endif

#include <benchmarks_sphere/SphereBenchmarksCombined.hpp>

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataPhysical.hpp>
#include <sweet/sphere/SphereDiagnostics.hpp>


#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>

#include "swe_sphere/SWE_Sphere_TimeSteppers.hpp"
#include "swe_sphere/SWE_Sphere_NormalModeAnalysis.hpp"



SimulationVariables simVars;

// Plane data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig sphereDataConfigInstance_nodealiasing;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig_nodealiasing = &sphereDataConfigInstance_nodealiasing;


#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif



/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */
//bool param_compute_error = false;


class SimulationInstance
{
public:
	SphereOperators op;
	SphereOperators op_nodealiasing;

	SWE_Sphere_TimeSteppers timeSteppers;


	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	SphereData prog_phi;
	SphereData prog_vort;
	SphereData prog_div;


	REXI_Terry<> rexi;

#if SWEET_GUI
	PlaneData viz_plane_data;
#endif

	int render_primitive_id = 1;

	SphereDiagnostics sphereDiagnostics;

#if SWEET_MPI
	int mpi_rank;
#endif


public:
	SimulationInstance()	:
		op(sphereDataConfig, simVars.sim.earth_radius),
		op_nodealiasing(sphereDataConfig_nodealiasing, simVars.sim.earth_radius),
		prog_phi(sphereDataConfig),
		prog_vort(sphereDataConfig),
		prog_div(sphereDataConfig),

#if SWEET_GUI
		viz_plane_data(planeDataConfig),
#endif
		sphereDiagnostics(sphereDataConfig, simVars)
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

		sphereDiagnostics.update_phi_vort_div_2_mass_energy_enstrophy(
				op,
				prog_phi,
				prog_vort,
				prog_div,
				simVars
		);
	}



	void reset()
	{
		simVars.reset();

		if (simVars.misc.sphere_use_robert_functions != 1)
			FatalError("Only Robert formulation allowed");

		// one month runtime
		if (simVars.timecontrol.max_simulation_time == -1)
		{
			if (simVars.setup.benchmark_scenario_id == 10)
			{
				simVars.timecontrol.max_simulation_time = 31*60*60*24;

				// 200 h
				simVars.timecontrol.max_simulation_time = 200*60*60;
			}
		}

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		simVars.misc.output_next_sim_seconds = 0;


		if (simVars.sim.CFL < 0)
			simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

#if 0
		if (simVars.timecontrol.current_timestep_size <= 0)
		{
			// TRY to guess optimal time step size

			// time step size
			if (sphereDataConfig->physical_num_lat < 256)
				simVars.timecontrol.current_timestep_size = 0.002*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
			else
				simVars.timecontrol.current_timestep_size = 0.001*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
		}
#endif

		if (simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Only fixed time step size supported");


		if (simVars.setup.benchmark_setup_dealiased)
		{
			// use dealiased physical space for setup
			SphereBenchmarksCombined::setupInitialConditions(prog_phi, prog_vort, prog_div, simVars, op);
		}
		else
		{
			// this is the default
			// use reduced physical space for setup to avoid spurious modes
			SphereData prog_phi_nodealiasing(sphereDataConfig_nodealiasing);
			SphereData prog_vort_nodealiasing(sphereDataConfig_nodealiasing);
			SphereData prog_div_nodealiasing(sphereDataConfig_nodealiasing);

			SphereBenchmarksCombined::setupInitialConditions(prog_phi_nodealiasing, prog_vort_nodealiasing, prog_div_nodealiasing, simVars, op_nodealiasing);

			prog_phi.load_nodealiasing(prog_phi_nodealiasing);
			prog_vort.load_nodealiasing(prog_vort_nodealiasing);
			prog_div.load_nodealiasing(prog_div_nodealiasing);
		}


#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			simVars.outputConfig();

#if 0
			std::cout << std::endl;
			std::cout << "LOCAL PARAMETERS:" << std::endl;
			std::cout << " + param_compute_error: " << param_compute_error << std::endl;
			std::cout << std::endl;
#endif
		}

		/*
		 * SETUP time steppers
		 */
		timeSteppers.setup(simVars.disc.timestepping_method, op, simVars);

		update_diagnostics();

		simVars.diag.backup_reference();
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const SphereData &i_sphereData,
			const char* i_name,		///< name of output variable
			bool i_phi_shifted
		)
	{
		char buffer[1024];

		// create copy
		SphereData sphereData(i_sphereData);

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;
	}



	void write_file_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return;
#endif

		if (simVars.misc.output_file_name_prefix.length() == 0)
			return;

		std::string output_filename;

		std::cout << "Simulation time: " << simVars.timecontrol.current_simulation_time << std::endl;

		SphereData h = prog_phi*(1.0/simVars.sim.gravitation);
		output_filename = write_file(h, "prog_h", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << " (min: " << SphereData(h).physical_reduce_min() << ", max: " << SphereData(h).physical_reduce_max() << ")" << std::endl;

		SphereDataPhysical u(sphereDataConfig);
		SphereDataPhysical v(sphereDataConfig);

		op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);

		output_filename = write_file(u, "prog_u", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;

		output_filename = write_file(v, "prog_v", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;

		output_filename = write_file(prog_vort, "prog_vort", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;

		SphereData potvort = (prog_phi/simVars.sim.gravitation)*prog_vort;

		output_filename = write_file(potvort, "prog_potvort", simVars.setup.benchmark_scenario_id == 0);
		std::cout << output_filename << std::endl;
	}



	void timestep_do_output()
	{
		if (simVars.misc.compute_errors)
		{
			if (
					simVars.setup.benchmark_scenario_name != "geostrophic_balance"		&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_1"		&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_2"		&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_4"		&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_8"		&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_16"	&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_32"	&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_64"	&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_128"	&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_256"	&&
					simVars.setup.benchmark_scenario_name != "geostrophic_balance_512"
			)
			{
				if (
						simVars.setup.benchmark_scenario_id != 10	 &&
						simVars.setup.benchmark_scenario_id != 11	 &&
						simVars.setup.benchmark_scenario_id != 101
				)
				{
					std::cout << "Benchamrk name: " << simVars.setup.benchmark_scenario_name << std::endl;
					std::cout << "Benchmark scenario id: " << simVars.setup.benchmark_scenario_id << std::endl;
					FatalError("Analytical solution not available for this benchmark");
				}
			}

			SphereData anal_solution_h(sphereDataConfig);
			SphereData anal_solution_u(sphereDataConfig);
			SphereData anal_solution_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(anal_solution_h, anal_solution_u, anal_solution_v, simVars, op);

			SphereDataPhysical anal_solution_hg = anal_solution_h.getSphereDataPhysical();
			SphereDataPhysical anal_solution_ug = anal_solution_u.getSphereDataPhysical();
			SphereDataPhysical anal_solution_vg = anal_solution_v.getSphereDataPhysical();

			double error_h = -1;
			double error_u = -1;
			double error_v = -1;
			double error_vort = -1;
			double error_div = -1;


			SphereData h = prog_phi*(1.0/simVars.sim.gravitation);
			SphereDataPhysical hg = h.getSphereDataPhysical();

			SphereData tmp_vort(sphereDataConfig);
			SphereData tmp_div(sphereDataConfig);


			/*
			 * Convert analytical correct velocities to vort/div
			 */
			if (simVars.misc.sphere_use_robert_functions)
				op.robert_uv_to_vortdiv(anal_solution_ug, anal_solution_vg, tmp_vort, tmp_div);
			else
				op.uv_to_vortdiv(anal_solution_ug, anal_solution_vg, tmp_vort, tmp_div);


			/*
			 * Compute difference
			 */
			SphereData diff_vort = prog_vort - tmp_vort;
			SphereData diff_div = prog_div - tmp_div;


			/*
			 * Convert back to u-v
			 */
			SphereDataPhysical diff_u(sphereDataConfig);
			SphereDataPhysical diff_v(sphereDataConfig);
			if (simVars.misc.sphere_use_robert_functions)
				op.robert_vortdiv_to_uv(diff_vort, diff_div, diff_u, diff_v);
			else
				op.vortdiv_to_uv(diff_vort, diff_div, diff_u, diff_v);

			error_h = hg.physical_reduce_max_abs(anal_solution_hg);
			error_u = diff_u.physical_reduce_max_abs();
			error_v = diff_v.physical_reduce_max_abs();
			error_vort = diff_vort.physical_reduce_max_abs();
			error_div = diff_div.physical_reduce_max_abs();
#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				std::cerr << "error time, h, u, v, vort, div:\t" << simVars.timecontrol.current_simulation_time << "\t" << error_h << "\t" << error_u << "\t" << error_v << "\t" << error_vort << "\t" << error_div << std::endl;
			}
		}

		write_file_output();

		// output line break
		std::cout << std::endl;

		if (simVars.misc.verbosity > 1)
		{
			update_diagnostics();

#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				// Print header
				if (simVars.timecontrol.current_timestep_nr == 0)
				{
					std::cerr << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
					std::cerr << std::endl;
				}

				// Print simulation time, energy and pot enstrophy
				std::cerr << simVars.timecontrol.current_simulation_time << "\t";
				std::cerr << simVars.diag.total_mass << "\t";
				std::cerr << simVars.diag.potential_energy << "\t";
				std::cerr << simVars.diag.kinetic_energy << "\t";
				std::cerr << simVars.diag.total_energy << "\t";
				std::cerr << simVars.diag.total_potential_enstrophy << "\t";

				std::cerr << (simVars.diag.total_mass-simVars.diag.ref_total_mass)/simVars.diag.total_mass << "\t";
				std::cerr << (simVars.diag.potential_energy-simVars.diag.ref_potential_energy)/simVars.diag.potential_energy << "\t";
				std::cerr << (simVars.diag.kinetic_energy-simVars.diag.ref_kinetic_energy)/simVars.diag.kinetic_energy << "\t";
				std::cerr << (simVars.diag.total_energy-simVars.diag.total_energy)/simVars.diag.total_energy << "\t";
				std::cerr << (simVars.diag.total_potential_enstrophy-simVars.diag.total_potential_enstrophy)/simVars.diag.total_potential_enstrophy << std::endl;

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
				std::cout << "prog_phi min/max:\t" << SphereData(prog_phi).physical_reduce_min() << ", " << SphereData(prog_phi).physical_reduce_max() << std::endl;
		}


#if 0
		if (	param_compute_error &&
			//	simVars.pde.use_nonlinear_equations == 0 &&
				simVars.setup.benchmark_scenario_id == 10
		)
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);

			std::cout << "ERRORS - time, RMS(h,u,v), MAXABS(h,u,v):\t";
			std::cout << simVars.timecontrol.current_simulation_time << "\t";

			std::cout << (SphereData(test_h)-SphereData(prog_h)).physical_reduce_rms() << "\t";
			std::cout << (SphereData(test_u)-SphereData(prog_u)).physical_reduce_rms() << "\t";
			std::cout << (SphereData(test_v)-SphereData(prog_v)).physical_reduce_rms() << "\t";

			std::cout << (SphereData(test_h)-SphereData(prog_h)).physical_reduce_max_abs() << "\t";
			std::cout << (SphereData(test_u)-SphereData(prog_u)).physical_reduce_max_abs() << "\t";
			std::cout << (SphereData(test_v)-SphereData(prog_v)).physical_reduce_max_abs() << std::endl;
		}
#endif

		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
	}



public:
	bool timestep_check_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return false;
#endif

		if (simVars.misc.verbosity > 0)
			std::cout << "." << std::flush;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		timestep_do_output();

		return true;
	}



public:
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



	bool detect_instability()
	{
		double max_abs_value = std::abs(simVars.sim.h0)*2.0*simVars.sim.gravitation;

		if (
				SphereData(prog_phi).physical_reduce_max_abs() > max_abs_value &&
				simVars.setup.benchmark_scenario_id != 4
		)
		{
			std::cerr << "Instability detected (max abs value of h > " << max_abs_value << ")" << std::endl;
			return true;
		}

		if (SphereData(prog_phi).physical_isAnyNaNorInf())
		{
			std::cerr << "Inf value detected" << std::endl;
			return true;
		}

		return false;
	}



	void run_timestep()
	{
#if SWEET_GUI
		if (simVars.misc.gui_enabled && simVars.disc.normal_mode_analysis_generation == 0)
			timestep_check_output();
#endif

		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_phi, prog_vort, prog_div,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);


		/*
		 * Add implicit viscosity
		 */
		if (simVars.sim.viscosity != 0)
		{
			double scalar = simVars.sim.viscosity*simVars.timecontrol.current_timestep_size;
			double r = simVars.sim.earth_radius;

			/*
			 * (1-dt*visc*D2)p(t+dt) = p(t)
			 */
			prog_phi = prog_phi.spectral_solve_helmholtz(1.0, -scalar, r);
			prog_vort = prog_vort.spectral_solve_helmholtz(1.0, -scalar, r);
			prog_div = prog_div.spectral_solve_helmholtz(1.0, -scalar, r);
		}

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
				prog_phi,
				prog_vort,
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
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data
	)
	{
		// request rendering of sphere
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		int id = simVars.misc.vis_id % 5;
		switch (id)
		{
		default:
		case 0:
			// USE COPY TO AVOID FORWARD/BACKWARD TRANSFORMATION
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_phi)/simVars.sim.gravitation, planeDataConfig);
			break;

		case 1:
		{
			SphereDataPhysical u(prog_vort.sphereDataConfig);
			SphereDataPhysical v(prog_vort.sphereDataConfig);

			// Don't use Robert, since we're not interested in the Robert formulation here
			op.vortdiv_to_uv(prog_vort, prog_div, u, v);
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(u, planeDataConfig);
			break;
		}

		case 2:
		{
			SphereDataPhysical u(prog_vort.sphereDataConfig);
			SphereDataPhysical v(prog_vort.sphereDataConfig);

			// Don't use Robert, since we're not interested in the Robert formulation here
			op.vortdiv_to_uv(prog_vort, prog_div, u, v);
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(v, planeDataConfig);
			break;
		}

		case 3:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_vort), planeDataConfig);
			break;

		case 4:
			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_div), planeDataConfig);
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		const char* description = "";

		int id = simVars.misc.vis_id % 4;

		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "U";
			break;

		case 2:
			description = "V";
			break;

		case 3:
			description = "vort";
			break;

		case 4:
			description = "div";
			break;
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
#if SWEET_MPI
				"Rank %i - "
#endif
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				mpi_rank,
#endif
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				viz_plane_data.reduce_max(),
				viz_plane_data.reduce_min()
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

#if 0
case 'C':
			// dump data arrays to VTK
			prog_h.file_physical_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h.file_physical_loadData("swe_rexi_dump_h.csv", simVars.setup_velocityformulation_progphiuv.input_data_binary);
			prog_u.file_physical_loadData("swe_rexi_dump_u.csv", simVars.setup_velocityformulation_progphiuv.input_data_binary);
			prog_v.file_physical_loadData("swe_rexi_dump_v.csv", simVars.setup_velocityformulation_progphiuv.input_data_binary);
			break;
#endif
		}
	}
#endif
};



int main(int i_argc, char *i_argv[])
{

#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_SPACE_THREADING
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

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			nullptr
	};


	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.printOptions();
#endif
#if SWEET_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		if (mpi_rank == 0)
#endif
		{
			std::cout << "	--compute-error [0/1]	Output errors (if available, default: 1)" << std::endl;
		}
		return -1;
	}

	sphereDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	int res_physical_nodealias[2] = {
			2*simVars.disc.res_spectral[0],
			simVars.disc.res_spectral[1]
		};

	sphereDataConfigInstance_nodealiasing.setupAuto(res_physical_nodealias, simVars.disc.res_spectral);


#if SWEET_GUI
	planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);



#if SWEET_MPI

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	std::cout << "Helo from MPI rank: " << mpi_rank << std::endl;

	// only start simulation and time stepping for first rank
	if (mpi_rank > 0)
	{
		/*
		 * Deactivate all output for ranks larger than the current one
		 */
		simVars.misc.verbosity = 0;
		simVars.misc.output_each_sim_seconds = -1;
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
			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller_Serial<SimulationInstance> parareal_Controller_Serial;

			// setup controller. This initializes several simulation instances
			parareal_Controller_Serial.setup(&simVars.parareal);

			// execute the simulation
			parareal_Controller_Serial.run();
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


			if (simVars.disc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normalmode_analysis();
			}
			else
			{
				// Do first output before starting timer
				simulationSWE->timestep_check_output();


				//Time counter
				Stopwatch time;

#if SWEET_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				// Start counting time
				if (mpi_rank == 0)
					std::cout << "TIMER RESET" << std::endl;
#endif
				time.reset();


				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Do the output here before the quit check.
					// This assures that in case of writing the data for the 1st and last time step,
					// this is not included in the timings
					simulationSWE->timestep_check_output();

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simVars.misc.stability_checks)
					{
						if (simulationSWE->detect_instability())
						{
							std::cout << "INSTABILITY DETECTED" << std::endl;
							break;
						}
					}
				}


				// Stop counting time
				time.stop();

#if SWEET_MPI
				MPI_Barrier(MPI_COMM_WORLD);
				// Start counting time
				if (mpi_rank == 0)
					std::cout << "TIMER STOP" << std::endl;
#endif

				double wallclock_time = time();

#if SWEET_MPI
				if (mpi_rank == 0)
#endif
				{
					// End of run output results
					std::cout << std::endl;
					std::cout << "***************************************************" << std::endl;
					std::cout << "Wallclock time (seconds): " << wallclock_time << std::endl;
					std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
					std::cout << "Time per time step: " << wallclock_time/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
					std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;
				}

				// Output final time step output!
				simulationSWE->timestep_check_output();
			}


			delete simulationSWE;
		}
	}


#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}
