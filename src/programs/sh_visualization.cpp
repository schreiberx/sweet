/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData_Physical.hpp>
	#include <sweet/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
	#include <sweet/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>
#endif

#include "swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereData_SpectralComplex.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereOperators_SphereDataComplex.hpp>
#include <sweet/SWEETError.hpp>
#include <sweet/sphere/SphereData_DebugContainer.hpp>


SimulationVariables simVars;


/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class SimulationInstance
{
public:
	SphereOperators_SphereData *ops;

	SphereData_Config *sphereDataConfig;
#if SWEET_GUI

	PlaneDataConfig *planeDataConfig;
#endif

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	SphereData_Spectral *prog_phi_pert = nullptr;
	SphereData_Spectral *prog_vrt = nullptr;
	SphereData_Spectral *prog_div = nullptr;

#if SWEET_GUI
	PlaneData_Physical *viz_plane_data = nullptr;
#endif

	int render_primitive_id = 1;

#if SWEET_MPI
	int mpi_rank;
#endif

	// was the output of the time step already done for this simulation state?
	double timestep_last_output_simtime;

	int mode_m = 0;
	int mode_n = 0;

	bool viz_reset = false;

public:
	SimulationInstance()
	{
#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
		reset();
	}

	~SimulationInstance()
	{
		cleanup();
	}

	void setup_mode()
	{
		SphereData_DebugContainer::clear();

		std::cout << std::endl;
		std::cout << "setup_mode() called" << std::endl;
		std::cout << " + m: " << mode_m << std::endl;
		std::cout << " + n: " << mode_n << std::endl;
		std::cout << std::endl;


		if (mode_n < 0 ||  mode_m < 0)
		{
			std::cout << "mode_n < 0 ||  mode_m < 0" << std::endl;
			return;
		}

		if (mode_n > sphereDataConfig->spectral_modes_n_max)
		{
			std::cout << "mode_n > sphereDataConfig->spectral_modes_n_max" << std::endl;
			return;
		}

		if (mode_m > sphereDataConfig->spectral_modes_m_max)
		{
			std::cout << "mode_m > sphereDataConfig->spectral_modes_m_max" << std::endl;
			return;
		}

		if (mode_m > mode_n)
		{
			std::cout << "ERROR: mode_m > mode_n" << std::endl;
			return;
		}


		SphereData_Spectral a = *prog_phi_pert;

		SphereData_DebugContainer::append(a, "a");

		SphereData_Spectral mu_a = ops->mu(a);
		SphereData_DebugContainer::append(mu_a, "mu(a)");

		SphereData_Spectral mug_a = ops->mug*a.toPhys();
		SphereData_DebugContainer::append(mug_a, "mug*a");

		SphereData_Spectral diff = mu_a - mug_a;
		SphereData_DebugContainer::append(diff, "diff");

		std::cout << "a" << std::endl;
		if (sphereDataConfig->spectral_modes_m_max <= 32)
			a.spectral_print(10, 1e-14);
		std::cout << "  maxabs: " << a.toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << std::endl;


		std::cout << "mu_a" << std::endl;
		if (sphereDataConfig->spectral_modes_m_max <= 32)
			mu_a.spectral_print(10, 1e-14);
		std::cout << "  maxabs: " << mu_a.toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << std::endl;

		std::cout << "mug_a" << std::endl;
		if (sphereDataConfig->spectral_modes_m_max <= 32)
			mug_a.spectral_print(10, 1e-14);
		std::cout << "  maxabs: " << mug_a.toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << std::endl;

		std::cout << "diff" << std::endl;
		if (sphereDataConfig->spectral_modes_m_max <= 32)
			diff.spectral_print(10, 1e-14);
		std::cout << "  maxabs: " << diff.toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << std::endl;
	}

	void cleanup()
	{
		delete prog_phi_pert;
		prog_phi_pert = nullptr;

		delete prog_vrt;
		prog_vrt = nullptr;

		delete prog_div;
		prog_div = nullptr;

		delete sphereDataConfig;
		sphereDataConfig = nullptr;

#if SWEET_GUI
		delete viz_plane_data;
		viz_plane_data = nullptr;
#endif

		delete ops;
		ops = nullptr;

		delete sphereDataConfig;
		sphereDataConfig = nullptr;

#if SWEET_GUI
		delete planeDataConfig;
		planeDataConfig = nullptr;
#endif

		SphereData_DebugContainer::clear();
	}


	void reset(bool online_reset = false)
	{
		cleanup();

		if (online_reset)
		{
			simVars.disc.space_res_physical[0] = 0;
			simVars.disc.space_res_physical[1] = 0;

			viz_reset = true;
		}
		else
		{
			simVars.reset();
		}

		sphereDataConfig = new SphereData_Config;
		sphereDataConfig->setupAuto(
				simVars.disc.space_res_physical,
				simVars.disc.space_res_spectral,
				simVars.misc.reuse_spectral_transformation_plans,
				simVars.misc.verbosity,
				simVars.parallelization.num_threads_space
			);

		std::cout << "SPH config string: " << sphereDataConfig->getConfigInformationString() << std::endl;

		ops = new SphereOperators_SphereData(sphereDataConfig, &(simVars.sim));
		prog_phi_pert = new SphereData_Spectral(sphereDataConfig);
		prog_vrt = new SphereData_Spectral(sphereDataConfig);
		prog_div = new SphereData_Spectral(sphereDataConfig);


#if SWEET_GUI
		planeDataConfig = new PlaneDataConfig;
		planeDataConfig->setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);

		viz_plane_data = new PlaneData_Physical(planeDataConfig);
#endif

		simVars.iodata.output_time_scale = 1.0/(60.0*60.0);

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		simVars.iodata.output_next_sim_seconds = 0;

		if (simVars.benchmark.benchmark_name == "m")
		{
			simVars.misc.vis_id = -1;

			setup_mode();

			std::cout << "Use m/M to change mode M" << std::endl;
			std::cout << "Use n/N to change mode N" << std::endl;
		}
		else if (simVars.benchmark.benchmark_name == "a")
		{
			SphereData_DebugContainer::clear();

			SphereData_Spectral one(sphereDataConfig);
			one.spectral_set_zero();
			one += 1;

			SphereData_DebugContainer::append(one, "one");
			SphereData_DebugContainer::append(ops->mug, "fg");

			SphereData_Spectral mu_one = ops->mu(one);
			SphereData_DebugContainer::append(mu_one, "mu(one)");

			SphereData_Spectral mug_one = ops->mug*one.toPhys();
			SphereData_DebugContainer::append(mug_one, "mug*one");

			SphereData_Spectral diff = mu_one - mug_one;
			SphereData_DebugContainer::append(diff, "diff");
		}
		else if (simVars.benchmark.benchmark_name == "b")
		{
			SphereData_DebugContainer::clear();

			SphereData_Spectral mu = ops->mug;

			SphereData_DebugContainer::append(mu, "mu");
			SphereData_DebugContainer::append(ops->mug, "mug");

			SphereData_Spectral mu_mu = ops->mu(mu);
			SphereData_DebugContainer::append(mu_mu, "mu(mu)");

			SphereData_Spectral mug_mu = ops->mug*mu.toPhys();
			SphereData_DebugContainer::append(mug_mu, "mug*mu");

			SphereData_Spectral diff = mu_mu - mug_mu;
			SphereData_DebugContainer::append(diff, "diff");
		}
		else
		{
			BenchmarksSphereSWE sphereBenchmarks;
			sphereBenchmarks.setup(simVars, *ops);
			sphereBenchmarks.master->get_initial_state(*prog_phi_pert, *prog_vrt, *prog_div);

			SphereData_Spectral a = *prog_phi_pert;

			SphereData_DebugContainer::append(a, "a");
			SphereData_DebugContainer::append(ops->mug, "mug");

			SphereData_Spectral mu_a = ops->mu(a);
			SphereData_DebugContainer::append(mu_a, "mu(a)");

			SphereData_Spectral mug_a = ops->mug*a.toPhys();
			SphereData_DebugContainer::append(mug_a, "mug*a");

			SphereData_Spectral diff = mu_a - mug_a;
			SphereData_DebugContainer::append(diff, "diff");
		}

		// start at one second in the past to ensure output at t=0
		timestep_last_output_simtime = simVars.timecontrol.current_simulation_time-1.0;

		/*
		 * Output configuration here to ensure that updated variables are included in this output
		 */
		simVars.outputConfig();
	}



public:
	bool should_quit()
	{
		return false;
	}


	void run_timestep()
	{
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
			bool *o_viz_reset
	)
	{
		*o_viz_reset = viz_reset;
		viz_reset = false;

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
					*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(d.data_spectral, planeDataConfig);
				else
					*viz_plane_data = Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(d.data_physical, planeDataConfig);

				*o_dataArray = viz_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
		}

		int id = simVars.misc.vis_id % max_viz_types;

		switch (id)
		{
			default:
			case 0:
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(*prog_phi_pert), planeDataConfig);
				break;

			case 1:
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(*prog_vrt), planeDataConfig);
				break;

			case 2:
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(*prog_div), planeDataConfig);
				break;

			case 3:
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(simVars.sim.h0 + SphereData_Spectral(*prog_phi_pert)/simVars.sim.gravitation, planeDataConfig);
				break;

			case 4:
			{
				SphereData_Physical u(prog_vrt->sphereDataConfig);
				SphereData_Physical v(prog_vrt->sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				ops->vrtdiv_to_uv(*prog_vrt, *prog_div, u, v);
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(u, planeDataConfig);
				break;
			}

			case 5:
			{
				SphereData_Physical u(prog_vrt->sphereDataConfig);
				SphereData_Physical v(prog_vrt->sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				ops->vrtdiv_to_uv(*prog_vrt, *prog_div, u, v);
				*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(v, planeDataConfig);
				break;
			}

			case 6:
			case 7:
			case 8:
			{
				SphereData_Spectral anal_solution_phi_pert(sphereDataConfig);
				SphereData_Spectral anal_solution_vort(sphereDataConfig);
				SphereData_Spectral anal_solution_div(sphereDataConfig);

				BenchmarksSphereSWE sphereBenchmarks;
				sphereBenchmarks.setup(simVars, *ops);
				sphereBenchmarks.master->get_initial_state(anal_solution_phi_pert, anal_solution_vort, anal_solution_div);

				switch (id)
				{
				case 6:
					*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(*prog_phi_pert - anal_solution_phi_pert, planeDataConfig);
					break;

				case 7:
					*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(*prog_vrt - anal_solution_vort, planeDataConfig);
					break;

				case 8:
					*viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(*prog_div - anal_solution_div, planeDataConfig);
					break;
				}
			}
		}


		double viz_min = viz_plane_data->physical_reduce_min();
		double viz_max = viz_plane_data->physical_reduce_max();

		viz_max = std::max(std::abs(viz_max), std::abs(viz_min));
		viz_min = -viz_max;

		*o_viz_min = viz_min;
		*o_viz_max = viz_max;


		*o_dataArray = viz_plane_data;
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
				description = "vort";
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
				description = "vort diff t0";
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
				"MaxVal: %.6e,"
				"MinVal: %.6e,"
				","
				"Time: %f secs,"
				"Time: %f hours,"
				"Time: %f days,"
				"timestep nr.: %i,"
				"timestep size: %f,"
				","
				"TMass: %.6e,"
				"TEnergy: %.6e,"
				"PotEnstrophy: %.6e,"
				","
				"Colorscale: lowest [Blue... green ... red] highest"
				,
#if SWEET_MPI
				mpi_rank,
#endif
				simVars.misc.vis_id,
				description.c_str(),
				viz_plane_data->reduce_max(),
				viz_plane_data->reduce_min(),

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

		case 'm':
			mode_m++;
			setup_mode();
			break;

		case 'M':
			mode_m--;
			setup_mode();
			break;

		case 'n':
			mode_n++;
			setup_mode();
			break;

		case 'N':
			mode_n--;
			setup_mode();
			break;

		case 't':
			simVars.disc.space_res_spectral[0] *= 2;
			simVars.disc.space_res_spectral[1] *= 2;
			reset(true);
			break;

		case 'T':
			simVars.disc.space_res_spectral[0] /= 2;
			simVars.disc.space_res_spectral[1] /= 2;
			reset(true);
			break;
		}
	}
#endif
};



int main(int i_argc, char *i_argv[])
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr))
		return -1;

	if (simVars.misc.verbosity > 3)
		std::cout << " + setup SH sphere transformations..." << std::endl;


#if SWEET_GUI
	if (simVars.misc.verbosity > 3)
		std::cout << " + setup FFT plane transformations..." << std::endl;

#endif

	if (simVars.misc.verbosity > 3)
		std::cout << " + setup finished" << std::endl;

	{

#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *sim = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(sim);
			delete sim;
		}
		else
#endif
		{
			SimulationInstance *sim = new SimulationInstance;

			{
				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (sim->should_quit())
						break;

					// Main call for timestep run
					sim->run_timestep();
				}


				if (simVars.misc.verbosity > 0)
					std::cout << std::endl;
			}

			std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;

			delete sim;
		}
	}

	return 0;
}
