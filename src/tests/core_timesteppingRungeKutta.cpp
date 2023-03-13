/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif

#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/ProgramArguments.hpp>

#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>




class SimulationTestRK
{
public:
	sweet::ErrorBase error;

	sweet::TimesteppingExplicitRKPlaneData timestepping;

	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneData_Config planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h;
		sweet::PlaneData_Spectral prog_u;
		sweet::PlaneData_Spectral prog_v;

		sweet::PlaneData_Physical prog_h_phys;
		sweet::PlaneData_Physical prog_u_phys;
		sweet::PlaneData_Physical prog_v_phys;


		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			if (!planeDataConfig.setupAuto(*i_shackPlaneDataOps))
				return error.forwardWithPositiveReturn(planeDataConfig.error);

			if (!ops.setup(
					planeDataConfig,
					*i_shackPlaneDataOps
				))
				return error.forwardWithPositiveReturn(ops.error);

			prog_h.setup(planeDataConfig);
			prog_u.setup(planeDataConfig);
			prog_v.setup(planeDataConfig);

			prog_h_phys.setup(planeDataConfig);
			prog_u_phys.setup(planeDataConfig);
			prog_v_phys.setup(planeDataConfig);

			return true;
		}

		void clear()
		{
			prog_h_phys.clear();
			prog_u_phys.clear();
			prog_v_phys.clear();

			prog_h.clear();
			prog_u.clear();
			prog_v.clear();

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Simulation data
	Data data;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackIOData *shackIOData;

	int function_order = -1;
	int timestepping_order = -1;

public:
	SimulationTestRK(
			int i_argc,
			char *const * const i_argv,
			int i_function_order,
			int i_timestepping_order
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr),
		shackTimestepControl(nullptr),
		function_order(i_function_order),
		timestepping_order(i_timestepping_order)
	{
		ERROR_CHECK_COND_RETURN(shackProgArgDict);
	}


	bool setup()
	{
		shackProgArgDict.setup();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackTimestepControl->validateMaxSimulationTimeOrTimestepNr();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		data.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(data);

		data.prog_h_phys.physical_set_all_value(test_function(function_order, 0));
		data.prog_u_phys.physical_set_all_value(0);
		data.prog_v_phys.physical_set_all_value(0);

		data.prog_h.loadPlaneDataPhysical(data.prog_h_phys);
		data.prog_u.loadPlaneDataPhysical(data.prog_u_phys);
		data.prog_v.loadPlaneDataPhysical(data.prog_v_phys);

		return true;
	}


	void clear()
	{
		data.clear();

		shackPlaneDataOps = nullptr;
		shackProgArgDict.clear();

		shackPlaneDataOps = nullptr;
		shackTimestepControl = nullptr;
	}


	/**
	 * Function of 4th order.
	 */
public:
	double test_function(
			int i_order,
			double z
	)
	{
		switch(i_order)
		{
		case 0:
			return	2.0;

		case 1:
			return	z
					+2.0;

		case 2:
			return	+(5./12.)*z*z
					+z
					+2.0;

		case 3:
			return	-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;

		case 4:
			return	(1./12.)*z*z*z*z
					-(1./2.)*z*z*z
					+(5./12.)*z*z
					+z
					+2.0;
		}

		SWEETError("Not supported (i_order)");
		return 0;
	}



public:
	double test_function_diff(int i_order, double z)
	{
		switch(i_order)
		{
		case 0:
			return	0.0;

		case 1:
			return	1.0;

		case 2:
			return	+(10./12.)*z
					+1.0;

		case 3:
			return	-(3./2.)*z*z
					+(10./12.)*z
					+1.0;

		case 4:
			return	(4./12.)*z*z*z
					-(3./2.)*z*z
					+(10./12.)*z
					+1.0;
		}
		SWEETError("Not supported test_function_diff (i_order)");
		return 0;
	}

	void p_run_euler_timestep_update(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_current_timestamp = -1
	)
	{
		sweet::PlaneData_Physical o_h_t_phys(data.planeDataConfig);
		sweet::PlaneData_Physical o_u_t_phys(data.planeDataConfig);
		sweet::PlaneData_Physical o_v_t_phys(data.planeDataConfig);

		o_h_t_phys.physical_set_all_value(test_function_diff(function_order, i_current_timestamp));
		o_u_t_phys.physical_set_all_value(0);
		o_v_t_phys.physical_set_all_value(0);

		o_h_t.loadPlaneDataPhysical(o_h_t_phys);
		o_u_t.loadPlaneDataPhysical(o_u_t_phys);
		o_v_t.loadPlaneDataPhysical(o_v_t_phys);

		shackTimestepControl->current_timestep_nr++;
	}


	void runTimestep()
	{
		// either set time step size to 0 for autodetection or to
		// a positive value to use a fixed time step size
		assert(shackTimestepControl->current_timestep_size > 0);

		shackTimestepControl->timestepHelperStart();

		timestepping.runTimestep(
				this,
				&SimulationTestRK::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				data.prog_h, data.prog_u, data.prog_v,
				shackTimestepControl->current_timestep_size,
				timestepping_order,
				shackTimestepControl->current_simulation_time
			);

		shackTimestepControl->timestepHelperEnd();
	}


	bool should_quit()
	{
		return false;
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{
	for (int fun_order = 0; fun_order <= 4; fun_order++)
	{
		for (int timestepping_order = 1; timestepping_order <= 4; timestepping_order++)
		{
			SimulationTestRK simulation(i_argc, i_argv, fun_order, timestepping_order);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

			simulation.setup();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(simulation);

			while(true)
			{
				if (simulation.shackTimestepControl->isFinalTimestepReached())
					break;

				simulation.runTimestep();


				double value_numerical = simulation.data.prog_h.toPhys().physical_get(0,0);
				double value_exact = simulation.test_function(simulation.function_order, simulation.shackTimestepControl->current_simulation_time);

				std::cout << "t=" << simulation.shackTimestepControl->current_simulation_time;
				std::cout << ", ";
				std::cout << "num=" << value_numerical;
				std::cout << ", ";
				std::cout << "exact=" << value_exact;
				std::cout << std::endl;
			}

			double value_numerical = simulation.data.prog_h.toPhys().physical_get(0,0);

			double value_exact = simulation.test_function(simulation.function_order, simulation.shackTimestepControl->current_simulation_time);

			double error = std::abs(value_numerical - value_exact);

			if (fun_order <= timestepping_order)
			{
				if (error > 1e-8)
				{
					std::cout << "ERROR threshold exceeded!" << std::endl;
					return 1;
				}

				std::cout << "OK" << std::endl;
			}
			else
			{
				std::cout << "OK (errors expected)" << std::endl;

				if (error < 1e-8)
				{
					std::cout << "WARNING: Error expected to be larger, however relatively small error found" << std::endl;
					return 1;
				}

			}
		}
	}


	return 0;
}
