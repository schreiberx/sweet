/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionPlane/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include "../programs/pde_advectionPlane/ProgramPDEAdvectionPlane.hpp"

int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionPlane sim(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(sim);

	sim.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(sim);

	// Simply test whether the clear and setup works properly
	sim.clear();
	sim.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(sim);


	int max_modes = 256;

	double dt = sim.shackTimestepControl->current_timestep_size;
	int initial_spectral_modes = sim.shackPlaneDataOps->space_res_spectral[0];

	int loop_counter = 0;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		// only clear and setup the spatial parts
		sim.clear();

		// setup stage 1 + 2
		sim.setup_1_registration();
		sim.setup_2_processArguments();

		/*
		 * Overwrite parameters
		 */
		sim.shackTimestepControl->current_timestep_size = dt/std::pow(2.0, loop_counter);

		/*
		 * We need higher resolution for
		 * - SL methods or if using
		 * - finite differences in space
		 */
		int expected_order;

		if (
				sim.shackTimeDisc->timestepping_method == "na_sl" ||
				sim.shackPlaneDataOps->space_use_spectral_basis_diffs == false
		)
		{
			sim.shackPlaneDataOps->space_res_spectral[0] = i;
			sim.shackPlaneDataOps->space_res_spectral[1] = i;

			sim.shackPlaneDataOps->space_res_physical[0] = 0;
			sim.shackPlaneDataOps->space_res_physical[1] = 0;

			// order is limited by the spatial interpolation order
			expected_order = std::min(sim.shackTimeDisc->timestepping_order, 2);
		}
		else
		{
			sim.shackPlaneDataOps->space_res_spectral[0] = initial_spectral_modes;
			sim.shackPlaneDataOps->space_res_spectral[1] = initial_spectral_modes;

			sim.shackPlaneDataOps->space_res_physical[0] = 0;
			sim.shackPlaneDataOps->space_res_physical[1] = 0;

			// order is directly the time discretization order
			expected_order = sim.shackTimeDisc->timestepping_order;
		}

		// setup 3rd part
		sim.setup_3_data();

		std::cout << "Testing with " << sim.data.planeDataConfig.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << sim.shackTimestepControl->current_timestep_size << std::endl;

		{
			while (!sim.should_quit())
				sim.run_timestep();

			double max_error = sim.getLMaxError();

			std::cout << "Lmax error compared to initial condition: " << max_error << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / max_error;
				std::cout << "prev_max_error: " << prev_max_error << std::endl;
				std::cout << "max_error: " << max_error << std::endl;
				std::cout << "Convergence: " << conv << std::endl;
				std::cout << "expected_order: " << expected_order << std::endl;

				if (conv < std::pow(2.0, (double)expected_order)*0.9)
				{
					std::cerr << "Convergence too low!" << std::endl;
					exit(1);
				}

				if (conv > std::pow(2.0, (double)(expected_order))*1.1)
				{
					std::cerr << "Convergence too high, stopping here!" << std::endl;
					exit(1);
				}
			}
			prev_max_error = max_error;
			std::cout << std::endl;
		}

		loop_counter++;
	}

	return 0;
}
