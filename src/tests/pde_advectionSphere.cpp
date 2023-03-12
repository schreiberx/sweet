/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */

#include "../programs/pde_advectionSphere/ProgramPDEAdvectionSphere.hpp"

int main(int i_argc, char *i_argv[])
{
	ProgramPDEAdvectionSphere progPDEAdvSphere(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(progPDEAdvSphere);

	progPDEAdvSphere.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(progPDEAdvSphere);

	// Simply test whether the clear and setup works properly
	progPDEAdvSphere.clear();
	progPDEAdvSphere.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(progPDEAdvSphere);


	int max_modes = 256;

	double dt = progPDEAdvSphere.shackTimestepControl->current_timestep_size;
	int initial_spectral_modes = progPDEAdvSphere.shackSphereDataOps->space_res_spectral[0];

	int loop_counter = 0;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		// only clear and setup the spatial parts
		progPDEAdvSphere.clear();

		// setup stage 1 + 2
		progPDEAdvSphere.setup_1_shackRegistration();
		progPDEAdvSphere.setup_2_processArguments();

		/*
		 * Overwrite parameters
		 */
		progPDEAdvSphere.shackTimestepControl->current_timestep_size = dt/std::pow(2.0, loop_counter);

		/*
		 * We need higher resolution for
		 * - SL methods or if using
		 * - finite differences in space
		 */
		int expected_order;

		if (
				progPDEAdvSphere.shackTimeDisc->timestepping_method == "na_sl"
		)
		{
			progPDEAdvSphere.shackSphereDataOps->space_res_spectral[0] = i;
			progPDEAdvSphere.shackSphereDataOps->space_res_spectral[1] = i;

			progPDEAdvSphere.shackSphereDataOps->space_res_physical[0] = 0;
			progPDEAdvSphere.shackSphereDataOps->space_res_physical[1] = 0;

			// order is limited by the spatial interpolation order
			expected_order = std::min(progPDEAdvSphere.shackTimeDisc->timestepping_order, 2);
		}
		else
		{
			progPDEAdvSphere.shackSphereDataOps->space_res_spectral[0] = initial_spectral_modes;
			progPDEAdvSphere.shackSphereDataOps->space_res_spectral[1] = initial_spectral_modes;

			progPDEAdvSphere.shackSphereDataOps->space_res_physical[0] = 0;
			progPDEAdvSphere.shackSphereDataOps->space_res_physical[1] = 0;

			// order is directly the time discretization order
			expected_order = progPDEAdvSphere.shackTimeDisc->timestepping_order;
		}

		// setup 3rd part
		progPDEAdvSphere.setup_3_dataOpsEtc();

		std::cout << "Testing with " << progPDEAdvSphere.dataConfigOps.sphereDataConfig.getUniqueIDString() << std::endl;
		std::cout << "Testing with dt=" << progPDEAdvSphere.shackTimestepControl->current_timestep_size << std::endl;

		{
			while (!progPDEAdvSphere.should_quit())
				progPDEAdvSphere.runTimestep();

			double max_error = progPDEAdvSphere.getErrorLMaxOnH();

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
