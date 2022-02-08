/*
 * test_sphere_parareal_sl.cpp
 *
 *  Created on: 08 Feb 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks
 */


#if SWEET_GUI
#	error	"GUI not supported"
#endif

#if !SWEET_PARAREAL
#	error	"Parareal must be enabled in this unit test"
#endif

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include "../programs/swe_sphere_timeintegrators/SWE_Sphere_TimeSteppers.hpp"
#include "../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"

SphereData_Config sphereDataConfigInstance;
SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;

SimulationVariables simVars;

class SimulationInstance
{
public:
	SphereData_Spectral prog_phi, prog_vrt, prog_div;

	SWE_Sphere_TimeSteppers timeSteppers;

	BenchmarksSphereSWE sphereBenchmarks;

	SphereOperators_SphereData op;

	int set_previous_solution;

public:
	SimulationInstance(int i_set_previous_solution)	:
		prog_phi(sphereDataConfig),
		prog_vrt(sphereDataConfig),
		prog_div(sphereDataConfig),

		op(sphereDataConfig, &(simVars.sim))

	{
		set_previous_solution = i_set_previous_solution;
		reset();
	}

	~SimulationInstance()
	{
	}

	void reset()
	{
		simVars.reset();

		sphereBenchmarks.setup(simVars, op);
		sphereBenchmarks.master->get_initial_state(prog_phi, prog_vrt, prog_div);

		// setup planeDataconfig instance again
		sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans, simVars.misc.verbosity);

		timeSteppers.setup(simVars.disc.timestepping_method,
				   op,
				   simVars);

		if (simVars.misc.verbosity > 2)
			simVars.outputConfig();
	}

	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_phi, prog_vrt, prog_div,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

		double dt = simVars.timecontrol.current_timestep_size;

		// advance in time
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

	}



	void run_simulation()
	{
		reset();

		timeSteppers.master->set_previous_solution(prog_phi, prog_vrt, prog_div);
		while (!should_quit()) {

			SphereData_Spectral phi_prev = prog_phi;
			SphereData_Spectral vrt_prev = prog_vrt;
			SphereData_Spectral div_prev = prog_div;
			run_timestep();

			if (set_previous_solution == 0)
			{
				// do nothing
			} else
			if (set_previous_solution == 1)
			{
				timeSteppers.master->set_previous_solution(phi_prev, vrt_prev, div_prev);
				// set correct previous solution
			} else if (set_previous_solution == 2)
			{
				// set wrong_previous_solution
				SphereData_Spectral phi2 = 0.5 * phi_prev;
				SphereData_Spectral vrt2 = 0.5 * vrt_prev;
				SphereData_Spectral div2 = 0.5 * div_prev;
				timeSteppers.master->set_previous_solution(phi2, vrt2, div2);
			}

		}
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


	double compute_error(
				SphereData_Spectral phi2,
				SphereData_Spectral vrt2,
				SphereData_Spectral div2)
	{
		return std::max(
				(prog_phi - phi2).spectral_reduce_max_abs(),
				std::max(
					(prog_vrt - vrt2).spectral_reduce_max_abs(),
					(prog_div - div2).spectral_reduce_max_abs()
					)
				);
	}

};

int main(int i_argc, char *i_argv[])
{

	// Testing function set_previous_solution(), which provides t_{n-1} for SL in parareal 

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if ( simVars.disc.timestepping_method.compare("lg_exp_na_sl_lc_nr_etd_uv") &&
	     simVars.disc.timestepping_method.compare("l_irk_na_sl_nr_settls_uv_only") &&
	     simVars.disc.timestepping_method.compare("l_irk_na_sl_nr_settls_vd_only") &&
	     simVars.disc.timestepping_method.compare("l_irk_na_sl_settls_uv_only") &&
	     simVars.disc.timestepping_method.compare("l_irk_na_sl_settls_vd_only") &&
	     simVars.disc.timestepping_method.compare("ln_sl_exp_settls_uv") &&
	     simVars.disc.timestepping_method.compare("ln_sl_exp_settls_vd") )
		SWEETError("This unit test should be executed with a SL timestepping method.");


	sphereDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);


	SimulationInstance sim_ref(0);
	SimulationInstance sim_correct(1);
	SimulationInstance sim_wrong(2);

	std::cout << "Running ref simulation" << std::endl;
	sim_ref.run_simulation();

	std::cout << "Running simulation with correct previous solution" << std::endl;
	sim_correct.run_simulation();

	std::cout << "Running simulation with wrong previous solution" << std::endl;
	sim_wrong.run_simulation();

	double eps = 1e-14;
	double error_correct = sim_ref.compute_error(sim_correct.prog_phi, sim_correct.prog_vrt, sim_correct.prog_div);
	double error_wrong = sim_ref.compute_error(sim_wrong.prog_phi, sim_wrong.prog_vrt, sim_wrong.prog_div);

	std::cout << "Error sim_correct: " << error_correct << std::endl;
	if ( error_correct > eps )
		SWEETError("Simulations should provide the same solution");
	else
		std::cout << " -- Passed" << std::endl;

	std::cout << "Error sim_wrong: " << error_wrong << std::endl;
	if ( error_wrong < eps )
		SWEETError("Simulations should not provide the same solution");
	else
		std::cout << " -- Passed" << std::endl;


	return 0;

}
