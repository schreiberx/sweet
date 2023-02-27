/*
 *  Created on: 2nd April 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/ProgramArguments.hpp>


#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>
#include <sweet/core/plane/Convert_PlaneDataPhysical_to_PlaneDataPhysicalComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataPhysicalComplex_to_PlaneDataPhysical.hpp>


class Core_planeSamplerInterpolation
{
public:
	sweet::ErrorBase error;


	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneDataConfig planeDataConfig;
		sweet::PlaneDataConfig planeDataConfigOversampling;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h;

		sweet::ScalarDataArray posx_a, posy_a;

		double *gaussianCenter;
		double gaussianExpFac = 50.0;

		sweet::PlaneDataSampler planeDataSampler;


		bool setup(
				sweet::ShackPlaneDataOps *i_shackPlaneDataOps,
				sweet::ShackPlaneDataOps *i_shackPlaneDataOpsOversampling,
				double i_gaussianCenter[]
		)
		{
			gaussianCenter = i_gaussianCenter;

			/*
			 * Setup Plane Data Config & Operators
			 */
			planeDataConfig.setupAuto(*i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(ops);

			planeDataConfigOversampling.setupAuto(*i_shackPlaneDataOpsOversampling);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(planeDataConfigOversampling);

			prog_h.setup(planeDataConfig);

			sweet::PlaneData_Physical prog_h_phys(planeDataConfig);

			prog_h_phys.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
					double x = (double)i*(i_shackPlaneDataOps->plane_domain_size[0]/(double)i_shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(i_shackPlaneDataOps->plane_domain_size[1]/(double)i_shackPlaneDataOps->space_res_physical[1]);

					io_data = gaussianValue(
							i_shackPlaneDataOps,
							gaussianCenter,
							x,
							y,
							gaussianExpFac
						);
				}
			);

			prog_h.loadPlaneDataPhysical(prog_h_phys);

			posx_a.setup(planeDataConfigOversampling.physical_array_data_number_of_elements);
			posy_a.setup(planeDataConfigOversampling.physical_array_data_number_of_elements);

			// setup some test sampling points
			// we use 2 arrays - one for each sampling position

			posx_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
				{
					int i = idx % planeDataConfigOversampling.physical_res[0];
					//int j = idx / planeDataConfig->physical_data_size[0];

					io_data = i_shackPlaneDataOps->plane_domain_size[0]*(double)i/(double)planeDataConfigOversampling.physical_res[0];
					assert(io_data >= 0);
					assert(io_data < i_shackPlaneDataOps->plane_domain_size[0]);
				}
			);
			posy_a.update_lambda_array_indices(
					[&](int idx, double &io_data)
				{
					//int i = idx % planeDataConfig->physical_data_size[0];
					int j = idx / planeDataConfigOversampling.physical_res[0];

					io_data = i_shackPlaneDataOps->plane_domain_size[1]*(double)j/(double)planeDataConfigOversampling.physical_res[1];

					assert(io_data >= -M_PI*0.5);
					assert(io_data < i_shackPlaneDataOps->plane_domain_size[1]);
				}
			);

			planeDataSampler.setup(i_shackPlaneDataOps->plane_domain_size, &planeDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();

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

	int interpolation_order = 3;

	double max_error;

public:
	Core_planeSamplerInterpolation(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr)
	{
		ERROR_CHECK_WITH_RETURN(shackProgArgDict);
	}

	bool setup(
			double i_gaussianCenter[],
			int i_interpolationOrder,
			int i_specModes[2],
			int i_specModesOversampled[2]
	)
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.setup();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

#if 0
		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);
#endif

		shackPlaneDataOps->space_res_spectral[0] = i_specModes[0];
		shackPlaneDataOps->space_res_spectral[1] = i_specModes[1];

		sweet::ShackPlaneDataOps shackPlaneDataOpsOversampling = *shackPlaneDataOps;

		shackPlaneDataOpsOversampling.space_res_spectral[0] = i_specModesOversampled[0];
		shackPlaneDataOpsOversampling.space_res_spectral[1] = i_specModesOversampled[1];

#if 0
		std::cout << "Shack Data regular:" << std::endl;
		shackPlaneDataOps->printShack("  ");

		std::cout << "Shack Data oversampled:" << std::endl;
		shackPlaneDataOpsOversampling.printShack("  ");
#endif
		data.setup(
				shackPlaneDataOps,
				&shackPlaneDataOpsOversampling,
				i_gaussianCenter
		);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(data);

		return true;
	}


	void clear()
	{
		data.clear();

		shackPlaneDataOps = nullptr;
		shackProgArgDict.clear();
	}


	static
	double gaussianValue(
			sweet::ShackPlaneDataOps *i_shackPlaneDataOps,
			double *i_center,
			double i_x, double i_y,
			double i_exp_fac
	)
	{
		double sx = i_shackPlaneDataOps->plane_domain_size[0];
		double sy = i_shackPlaneDataOps->plane_domain_size[1];

		// Gaussian
		double dx = i_x-i_center[0]*sx;
		double dy = i_y-i_center[1]*sy;

		if (dx > 0.5*i_shackPlaneDataOps->plane_domain_size[0])
			dx -= i_shackPlaneDataOps->plane_domain_size[0];
		else if (dx < -0.5*i_shackPlaneDataOps->plane_domain_size[0])
			dx += i_shackPlaneDataOps->plane_domain_size[0];

		if (dy > 0.5*i_shackPlaneDataOps->plane_domain_size[1])
			dy -= i_shackPlaneDataOps->plane_domain_size[1];
		else if (dy < -0.5*i_shackPlaneDataOps->plane_domain_size[1])
			dy += i_shackPlaneDataOps->plane_domain_size[1];

		dx /= sx;
		dy /= sy;

		return std::exp(-i_exp_fac*(dx*dx + dy*dy));
	}



	void run_tests()
	{
		sweet::ScalarDataArray out_data;
		out_data.setup(data.posx_a.number_of_elements);

		if (interpolation_order == 2)
		{
			data.planeDataSampler.bilinear_scalar(
					data.prog_h.toPhys(),
					data.posx_a,
					data.posy_a,
					out_data
			);
		}
		else if (interpolation_order == 3)
		{
			data.planeDataSampler.bicubic_scalar(
					data.prog_h.toPhys(),
					data.posx_a,
					data.posy_a,
					out_data
			);
		}
		else
		{
			SWEETError("Interpolation order not available");
		}

		/*
		 * Compute errors
		 */
		max_error = 0;
		for (std::size_t i = 0; i < data.posx_a.number_of_elements; i++)
		{
			double value = gaussianValue(
					shackPlaneDataOps,
					data.gaussianCenter,
					data.posx_a.scalar_data[i],
					data.posy_a.scalar_data[i],
					data.gaussianExpFac
				);

			max_error = std::max(max_error, std::abs(value - out_data.scalar_data[i]));
		}
	}


	bool should_quit()
	{
		return false;
	}
};




int main(int i_argc, char *i_argv[])
{
	/*
	 * SHACK: Register classes which we require
	 */

	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackPlaneDataOps *shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(shackProgArgDict);


	int initialSpectralModes = shackPlaneDataOps->space_res_spectral[0];

	double gaussianCenterArray[6][2] = {
			{0.5, 0.5},
			{0.9, 0.4},
			{0.4, 0.9},
			{0.3, 0.8},
			{0.7, 0.3},
			{0.3, 0.7},
	};

	for (int i_gaussians = 2; i_gaussians < 6; i_gaussians++)
	{
		std::cout << "*********************************************************" << std::endl;
		std::cout << "* Running studies for Gaussian at " << gaussianCenterArray[i_gaussians][0] << ", " << gaussianCenterArray[i_gaussians][1] << std::endl;
		std::cout << "*********************************************************" << std::endl;

		for (int interpolation_order = 2; interpolation_order <= 3; interpolation_order++)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* Running studies for interpolation of order " << interpolation_order << std::endl;
			std::cout << "*********************************************************" << std::endl;

			int oversamplingFactor = 13;
			std::cout << "Using oversampling of " << oversamplingFactor << std::endl;

			double prev_max_error = -1;
			for (int specModes = initialSpectralModes; specModes <= 256; specModes *= 2)
			{

				Core_planeSamplerInterpolation simulation(i_argc, i_argv);
				ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

				int specModes_[2] = {specModes, specModes};
				int specModesOversampled_[2] = {specModes*oversamplingFactor, specModes*oversamplingFactor};

				simulation.setup(
						gaussianCenterArray[i_gaussians],
						interpolation_order,
						specModes_,
						specModesOversampled_
					);
				ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

				simulation.run_tests();
				ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

				std::cout << "Error: " << simulation.max_error << std::endl;
				if (prev_max_error >= 0)
				{
					if (std::isnan(simulation.max_error))
						SWEETError("NaN detected");

					//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
					double conv = prev_max_error / simulation.max_error;
					std::cout << "Convergence: " << conv << std::endl;

					if (conv*1.1 < std::pow(2.0, interpolation_order))
						SWEETError("Convergence not given!");
				}
				prev_max_error = simulation.max_error;

				simulation.clear();
				ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

			}
		}
	}


	std::cout << "FIN" << std::endl;
	return 0;
}
