/*
 * test_plane_parareal_sl.cpp
 *
 *  Created on: 08 Feb 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_plane_timeintegrators
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_plane_benchmarks
 */


#if SWEET_GUI
#	error	"GUI not supported"
#endif

#if !SWEET_PARAREAL
#	error	"Parareal must be enabled in this unit test"
#endif

#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "../programs/swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include "../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"

PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

SimulationVariables simVars;


class SimulationInstance
{
public:
	PlaneData_Spectral prog_h, prog_u, prog_v;

	SWE_Plane_TimeSteppers timeSteppers;

	SWEPlaneBenchmarksCombined planeBenchmarkCombined;

	PlaneOperators op;

	int set_previous_solution;

public:
	SimulationInstance(int i_set_previous_solution)	:
		prog_h(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)

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

		planeBenchmarkCombined.setupInitialConditions(
				prog_h, prog_u, prog_v,
				simVars, op
			);

		// setup planeDataconfig instance again
		planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

		timeSteppers.setup(simVars.disc.timestepping_method,
				   simVars.disc.timestepping_order,
				   simVars.disc.timestepping_order2,
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
				prog_h, prog_u, prog_v,
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

		timeSteppers.master->set_previous_solution(prog_h, prog_u, prog_v);
		while (!should_quit()) {

			PlaneData_Spectral h_prev = prog_h;
			PlaneData_Spectral u_prev = prog_u;
			PlaneData_Spectral v_prev = prog_v;
			run_timestep();

			if (set_previous_solution == 0)
			{
				// do nothing
			} else
			if (set_previous_solution == 1)
			{
				timeSteppers.master->set_previous_solution(h_prev, u_prev, v_prev);
				// set correct previous solution
			} else if (set_previous_solution == 2)
			{
				// set wrong_previous_solution
				PlaneData_Spectral h2 = 0.5 * h_prev;
				PlaneData_Spectral u2 = 0.5 * u_prev;
				PlaneData_Spectral v2 = 0.5 * v_prev;
				timeSteppers.master->set_previous_solution(h2, u2, v2);
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
				PlaneData_Spectral h2,
				PlaneData_Spectral u2,
				PlaneData_Spectral v2)
	{
		return std::max(
				(prog_h - h2).reduce_maxAbs(),
				std::max(
					(prog_u - u2).reduce_maxAbs(),
					(prog_v - v2).reduce_maxAbs()
					)
				);
	}


};

int main(int i_argc, char *i_argv[])
{

	// Testing function set_previous_solution(), which provides t_{n-1} for SL in parareal 

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if ( simVars.disc.timestepping_method.compare("l_cn_na_sl_nd_settls") &&
	     simVars.disc.timestepping_method.compare("l_rexi_na_sl_nd_etdrk") &&
	     simVars.disc.timestepping_method.compare("l_rexi_na_sl_nd_settls")  )
		SWEETError("This unit test should be executed with a SL timestepping method.");


	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);


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
	double error_correct = sim_ref.compute_error(sim_correct.prog_h, sim_correct.prog_u, sim_correct.prog_v);
	double error_wrong = sim_ref.compute_error(sim_wrong.prog_h, sim_wrong.prog_u, sim_wrong.prog_v);

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
