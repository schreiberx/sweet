/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_advectionSphere/benchmarks/
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */


#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/time/TimesteppingSemiLagrangianSphereData.hpp>
#include <sweet/core/time/ShackTimesteppingSemiLagrangianSphereData.hpp>


#include "../programs/pde_advectionSphere/PDEAdvectionSphereBenchmarksCombined.hpp"
#include "../programs/pde_advectionSphere/PDEAdvectionSphereTimeSteppers.hpp"
#include "../programs/pde_advectionSphere/time/ShackPDEAdvectionSphereTimeDiscretization.hpp"


#include "../programs/pde_advectionSphere/ProgramPDEAdvectionSphere.hpp"


int main(int i_argc, char *i_argv[])
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackTimestepControl *shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackTimesteppingSemiLagrangianSphereData *shackTimesteppingSemiLagrangianSphereData = shackProgArgDict.getAutoRegistration<sweet::ShackTimesteppingSemiLagrangianSphereData>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	ShackPDEAdvectionSphereTimeDiscretization *shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphereTimeDiscretization>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	int initial_spectral_modes = shackSphereDataOps->space_res_spectral[0];

	if (shackTimestepControl->current_timestep_size < 0)
		SWEETError("Timestep size not set");

	int max_modes = 256;

	if (shackTimeDisc->timestepping_order == 1)
		max_modes = 512;
	else if (shackTimeDisc->timestepping_order == 2)
		max_modes = 256;

	double prev_max_error = -1;
	for (int i = initial_spectral_modes; i <= max_modes; i *= 2)
	{
		shackTimestepControl->current_timestep_size *= 0.5;
		//shackTimestepControl->setup_timestep_size = shackTimestepControl->current_timestep_size;

		if (shackTimeDisc->timestepping_method == "na_sl")
		{
			shackSphereDataOps->space_res_spectral[0] = i;
			shackSphereDataOps->space_res_spectral[1] = i;

			shackSphereDataOps->space_res_physical[0] = 2*i;
			shackSphereDataOps->space_res_physical[1] = i;
		}
		else
		{
			shackSphereDataOps->space_res_spectral[0] = initial_spectral_modes;
			shackSphereDataOps->space_res_spectral[1] = initial_spectral_modes;

			shackSphereDataOps->space_res_physical[0] = 0;
			shackSphereDataOps->space_res_physical[1] = 0;
		}


		if (1)
		{
			sweet::SphereData_Config sphereDataConfigInstance;
			sphereDataConfigInstance.setupAuto(
					shackSphereDataOps
			);
			ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(sphereDataConfigInstance);

			sweet::SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;

			std::cout << "Checking for right velocity" << std::endl;
			sweet::TimesteppingSemiLagrangianSphereData sl;
			sl.setup(sphereDataConfig, shackTimesteppingSemiLagrangianSphereData, shackTimeDisc->timestepping_order);


			/*
			 * Convert to Cartesian velocity space
			 */
			sweet::ScalarDataArray V_lon_D(sphereDataConfig->physical_array_data_number_of_elements);
			sweet::ScalarDataArray V_lat_D(sphereDataConfig->physical_array_data_number_of_elements);
			V_lon_D.set_all(1.0);
			V_lat_D.set_all(2.0);

			sweet::ScalarDataArray V_x_D, V_y_D, V_z_D;
			sweet::VectorMath::velocity_latlon_to_cartesian__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_lon_D,
					V_lat_D,
					V_x_D,
					V_y_D,
					V_z_D
				);

			sweet::ScalarDataArray V_lon_tmp, V_lat_tmp;
			sweet::VectorMath::velocity_cartesian_to_latlon__array(
					sl.pos_lon_A,
					sl.pos_lat_A,
					V_x_D,
					V_y_D,
					V_z_D,
					V_lon_tmp, V_lat_tmp
			);

			double err_lon = (V_lon_D - V_lon_tmp).reduce_maxAbs();
			double err_lat = (V_lat_D - V_lat_tmp).reduce_maxAbs();

			if (err_lon > 1e-10)
			{
				std::cerr << "Error: " << err_lon << std::endl;
				SWEETError("Error lon too high!");
			}

			if (err_lat > 1e-10)
			{
				std::cerr << "Error: " << err_lat << std::endl;
				SWEETError("Error lat too high!");
			}
		}

		std::cout << "Testing with SphereDataOps:" << std::endl;
		shackSphereDataOps->printShack();
		std::cout << std::endl;

		std::cout << "Testing with dt=" << shackTimestepControl->current_timestep_size << std::endl;

		ProgramPDEAdvectionSphere simulation(i_argc, i_argv);

		simulation.setup_1_shackRegistration();
		simulation.setup_2_processArguments();

		// physical resolution
		simulation.shackSphereDataOps->space_res_physical[0] = shackSphereDataOps->space_res_physical[0];
		simulation.shackSphereDataOps->space_res_physical[1] = shackSphereDataOps->space_res_physical[1];

		// spectral resolution
		simulation.shackSphereDataOps->space_res_spectral[0] = shackSphereDataOps->space_res_spectral[0];
		simulation.shackSphereDataOps->space_res_spectral[1] = shackSphereDataOps->space_res_spectral[1];

		// time step size
		simulation.shackTimestepControl->current_timestep_size = shackTimestepControl->current_timestep_size;

		// Setup the data and operators
		simulation.setup_3_dataOpsEtc();

		{
			while (!simulation.should_quit())
				simulation.runTimestep();

			double error_lmax = simulation.getErrorLMaxOnH();
			double error_rms = simulation.getErrorRMSOnH();

			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << error_lmax << std::endl;
			std::cout << "RMS error: " << error_lmax << std::endl;

			if (prev_max_error >= 0)
			{
				//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
				double conv = prev_max_error / error_lmax;
				std::cout << "Convergence: " << conv << std::endl;

				if (conv*1.1 < std::pow(2.0, (double)shackTimeDisc->timestepping_order))
					SWEETError("Convergence not given!");
			}

			if (error_lmax  > 1e10)
				SWEETError("Lmax error exceeded threshold!");

			prev_max_error = error_lmax;

			std::cout << "*********************************************" << std::endl;
		}
	}

	return 0;
}
