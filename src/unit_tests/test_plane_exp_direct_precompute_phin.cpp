/*
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_plane_timeintegrators/
 */

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include "../include/sweet/plane/PlaneData_Spectral.hpp"
#include "../include/sweet/plane/PlaneData_Physical.hpp"
#include "../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include "../programs/swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;

class SimulationInstance
{
public:
	PlaneData_Spectral prog_h_pert, prog_u, prog_v;

	// Initial values for comparison with analytical solution
	PlaneData_Spectral t0_prog_h_pert, t0_prog_u, t0_prog_v;

	// Forcings
	PlaneData_Spectral force_h_pert, force_u, force_v;

	PlaneDataGridMapping gridMapping;

	SWE_Plane_TimeSteppers timeSteppers;

	//Swe benchmarks
	SWEPlaneBenchmarksCombined swePlaneBenchmarks;

	PlaneOperators op;

#if SWEET_GUI
	PlaneData_Spectral viz_plane_data;

	int render_primitive_id = 0;
#endif

	SWEPlaneBenchmarksCombined planeBenchmarkCombined;

public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h_pert(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

#if SWEET_GUI
		vis(planeDataConfig),
#endif

		t0_prog_h_pert(planeDataConfig),
		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

		force_h_pert(planeDataConfig),
		force_u(planeDataConfig),
		force_v(planeDataConfig),

		// Initialises operators
		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		// Calls initialisation of the run (e.g. sets u, v, h)
		reset();
	}

	virtual ~SimulationInstance()
	{
	}

	void reset()
	{
		SimulationBenchmarkTimings::getInstance().main_setup.start();

		simVars.reset();

		if (simVars.benchmark.benchmark_name == "")
		{
			std::cout << "Benchmark scenario not selected (option --benchmark-name [string])" << std::endl;
			swePlaneBenchmarks.printBenchmarkInformation();
			SWEETError("Benchmark name not given");
		}

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		PlaneData_Physical prog_h_pert_phys(planeDataConfig);
		PlaneData_Physical prog_u_phys(planeDataConfig);
		PlaneData_Physical prog_v_phys(planeDataConfig);

		// set to some values for first touch NUMA policy (HPC stuff)
		prog_h_pert_phys.physical_set_all_value(simVars.sim.h0);
		prog_h_pert.loadPlaneDataPhysical(prog_h_pert_phys);
		prog_u.spectral_set_zero();
		prog_v.spectral_set_zero();

		// Setup prog vars
		//prog_u.physical_set_all(0);
		//prog_v.physical_set_all(0);

		// Check if input parameters are adequate for this simulation
		if (simVars.disc.space_grid_use_c_staggering && simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.space_grid_use_c_staggering ||  !simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
#endif

		if (simVars.disc.space_grid_use_c_staggering)
			gridMapping.setup(simVars, planeDataConfig);

		swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, simVars, op);

		prog_h_pert = t0_prog_h_pert;
		prog_u = t0_prog_u;
		prog_v = t0_prog_v;
		
		// Load data, if requested
		if (simVars.iodata.initial_condition_data_filenames.size() > 0)
		{
			prog_h_pert_phys.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[0].c_str(), simVars.iodata.initial_condition_input_data_binary);
			prog_h_pert.loadPlaneDataPhysical(prog_h_pert_phys);
		}

		if (simVars.iodata.initial_condition_data_filenames.size() > 1)
		{
			prog_u_phys.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[1].c_str(), simVars.iodata.initial_condition_input_data_binary);
			prog_u.loadPlaneDataPhysical(prog_u_phys);
		}

		if (simVars.iodata.initial_condition_data_filenames.size() > 2)
		{
			prog_v_phys.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[2].c_str(), simVars.iodata.initial_condition_input_data_binary);
			prog_v.loadPlaneDataPhysical(prog_v_phys);
		}

		timeSteppers.setup(
				simVars.disc.timestepping_method,
				simVars.disc.timestepping_order,
				simVars.disc.timestepping_order2,
				op,
				simVars
			);

		SimulationBenchmarkTimings::getInstance().main_setup.stop();
	}


	/**
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_h_pert, prog_u, prog_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

		// Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral

		if (simVars.sim.viscosity != 0 && simVars.misc.use_nonlinear_only_visc == 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE //TODO: this needs checking

			double dt = simVars.timecontrol.current_timestep_size;

			prog_u = prog_u + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_v = prog_v + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_h_pert = prog_h_pert + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity;
#else
			prog_u = op.implicit_diffusion(prog_u, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_v = op.implicit_diffusion(prog_v, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_h_pert = op.implicit_diffusion(prog_h_pert, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}

		// advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			SWEETError("Max simulation time exceeded!");

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
		return -1;

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");

	int initial_spectral_modes = simVars.disc.space_res_spectral[0];
	if (initial_spectral_modes <= 0)
	{
		SWEETError("Please specify the number of MODES");
	}

	if (simVars.timecontrol.current_timestep_size < 0)
		SWEETError("Timestep size not set");


	double dt = simVars.timecontrol.current_timestep_size;
	double Tmax = simVars.timecontrol.max_simulation_time;

	simVars.misc.verbosity = 6;

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);
	std::vector<SimulationInstance*> simulations;

        for (int i = 0; i < 4; i++)
	{

		simulations.push_back(new SimulationInstance);
		planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);


		if (i == 0)
		{
			simVars.rexi.exp_direct_precompute_phin = 0;
		}
		else if (i == 1)
		{
			simVars.rexi.exp_direct_precompute_phin = 1;
		}
		else if (i == 2)
		{
			simVars.rexi.exp_direct_precompute_phin = 0;
			simVars.timecontrol.max_simulation_time = Tmax + simVars.timecontrol.current_timestep_size / 2.;
		}
		else if (i == 3)
		{
			simVars.rexi.exp_direct_precompute_phin = 1;
			simVars.timecontrol.max_simulation_time = Tmax + simVars.timecontrol.current_timestep_size / 2.;
		}

		std::cout << std::endl;
		std::cout << " --> Running simulation with Tmax = " << simVars.timecontrol.max_simulation_time << "; precompute = " << simVars.rexi.exp_direct_precompute_phin << std::endl;
		while (!simulations[i]->should_quit())
			simulations[i]->run_timestep();
		std::cout << std::endl;

	}

	// Check if solutions are the same with and without precomputing
	double eps = 1e-13;
	for (int i = 0; i < 2; ++i)
	{
		double err_h = (simulations[2 * i]->prog_h_pert - simulations[2 * i + 1]->prog_h_pert).toPhys().physical_reduce_max_abs();
		double err_u = (simulations[2 * i]->prog_u - simulations[2 * i + 1]->prog_u).toPhys().physical_reduce_max_abs();
		double err_v = (simulations[2 * i]->prog_v - simulations[2 * i + 1]->prog_v).toPhys().physical_reduce_max_abs();

		std::cout << " --> Comparing solutions with Tmax = " << Tmax + i * simVars.timecontrol.current_timestep_size / 2. << ", dt = " << simVars.timecontrol.current_timestep_size << std::endl;
		std::cout << "  ** Error with / without precomputing phin:" << std::endl;
		std::cout << "    * Error h: " << err_h << std::endl;
		std::cout << "    * Error u: " << err_u << std::endl;
		std::cout << "    * Error v: " << err_v << std::endl;

		assert(err_h < eps);
		assert(err_u < eps);
		assert(err_v < eps);

		std::cout << "TEST OK" << std::endl << std::endl;

	}

	for (int i = 0; i < 4; ++i)
	{
		delete simulations[i];
		simulations[i] = nullptr;
	}

	return 0;
}
