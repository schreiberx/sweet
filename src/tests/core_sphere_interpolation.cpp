/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
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



class InterpolationTests
{
public:
	sweet::SphereData_Config *sphereDataConfig;
	sweet::SphereData_Config *sphereDataConfigOversampled;

	sweet::SphereData_Spectral prog_h;

	int interpolation_order = 3;

	bool use_limiter = false;
	bool use_poles_pseudo_points = false;

	sweet::ScalarDataArray posx_a, posy_a;

	#define MAX_GAUSSIANS 9
	double gaussian_center_array[MAX_GAUSSIANS][2] = {
			{0.0, M_PI*0.5},	// top
			{0.0, 0.0},		// equator
			{0.0, -M_PI*0.5},	// bottom

			{1.1, M_PI*0.5*0.8},	// a-kind top
			{3.0, 0.1},		// a-kind equator
			{2.0, -M_PI*0.5*0.75},	// a-kind bottom

			{0.1, M_PI*0.3},	// misc
			{1.0, M_PI*0.4},	// misc
			{2.3, -M_PI*0.34},	// misc
	};

	int gaussian_id = 0;
	//double center_lon = 0;
	//double center_lat = 0;
	double exp_fac = 20.0;

	sweet::SphereOperators_Sampler_SphereDataPhysical sphereDataSampler;


public:
	InterpolationTests()	:
		sphereDataConfig(nullptr),
		sphereDataConfigOversampled(nullptr)
	{
	}

	void setup(
			sweet::SphereData_Config &i_sphereDataConfig,
			sweet::SphereData_Config &i_sphereDataConfigOversampled
	)
	{
		sphereDataConfig = &i_sphereDataConfig;
		sphereDataConfigOversampled = &i_sphereDataConfigOversampled;

		prog_h.setup(sphereDataConfig);

		sweet::SphereData_Spectral tmp_vort(sphereDataConfig);
		sweet::SphereData_Spectral tmp_div(sphereDataConfig);

		sweet::SphereData_Physical prog_h_phys(sphereDataConfig);

		prog_h_phys.physical_update_lambda(
			[&](double i_lon, double i_lat, double &o_data)
			{
				o_data = gaussianValue(i_lon, i_lat, exp_fac);
			}
		);
		prog_h.loadSphereDataPhysical(prog_h_phys);

		posx_a.setup(sphereDataConfigOversampled->physical_array_data_number_of_elements);
		posy_a.setup(sphereDataConfigOversampled->physical_array_data_number_of_elements);

		// setup some test sampling points
		// we use 2 arrays - one for each sampling position
		posx_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % sphereDataConfigOversampled->physical_num_lon;

				io_data = 2.0*M_PI*(double)i/(double)sphereDataConfigOversampled->physical_num_lon;
				assert(io_data >= 0);
				assert(io_data < 2.0*M_PI);
			}
		);
		posy_a.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				//int i = idx % sphereDataConfig->physical_data_size[0];
				int j = idx / sphereDataConfigOversampled->physical_num_lon;

				io_data = sphereDataConfigOversampled->lat[j];

				assert(io_data >= -M_PI*0.5);
				assert(io_data <= M_PI*0.5);
			}
		);

		sphereDataSampler.setup(sphereDataConfig);
	}


	double gaussianValue(
			double i_lon, double i_lat,
			double i_exp_fac
	)
	{
		if (gaussian_id >= 0)
		{
			return gaussianValue_(gaussian_center_array[gaussian_id][0], gaussian_center_array[gaussian_id][1], i_lon, i_lat, i_exp_fac);
		}
		else
		{
			double o_data = 0;

			for (int i_gaussians = 0; i_gaussians < MAX_GAUSSIANS; i_gaussians++)
			{
				double scalar = std::cos(i_lat*2)*std::cos(i_lon*4);
				scalar = 1.0;
				o_data += scalar*gaussianValue_(gaussian_center_array[i_gaussians][0], gaussian_center_array[i_gaussians][1], i_lon, i_lat, i_exp_fac);
				//o_data = scalar;
			}

			o_data /= MAX_GAUSSIANS;
			return o_data;
		}
	}

	double gaussianValue_(
			double i_center_lon, double i_center_lat,
			double i_lon, double i_lat,
			double i_exp_fac
	)
	{
#if 1

		double x0[3];
		sweet::VectorMath::point_latlon_to_cartesian__scalar(i_center_lon, i_center_lat, x0[0], x0[1], x0[2]);

		double x[3];
		sweet::VectorMath::point_latlon_to_cartesian__scalar(i_lon, i_lat, x[0], x[1], x[2]);

#if 0
		double d =	(x[0] - x0[0])*(x[0] - x0[0])*(2.0+i_lon*0.1) +
					(x[1] - x0[1])*(x[1] - x0[1])*(2.0+i_lat*0.1) +
					(x[2] - x0[2])*(x[2] - x0[2])*(2.0+i_lon*i_lat*0.1);
#else

		double d =	(x[0] - x0[0])*(x[0] - x0[0])*(2.0) +
					(x[1] - x0[1])*(x[1] - x0[1])*(3.0) +
					(x[2] - x0[2])*(x[2] - x0[2])*(4.0);
#endif

		return std::exp(-20*d);

#else
		double center_lon = i_center_lon;
		double center_lat = i_center_lat;

		double mu = std::sin(i_lat);
		double phi1 = asin(mu);
		double phi2 = center_lat;
		double lambda1 = i_lon;
		double lambda2 = center_lon;

		double d = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		double d1 = acos(sin(phi1)*sin(phi2));
		double d2 = acos(cos(phi1)*cos(phi2)*cos(lambda1-lambda2));

		//return std::exp(-d*d*i_exp_fac)*(1.0-(d1*d1)/(M_PI*M_PI));//*0.1*simVars.sim.h0;// + simVars.sim.h0;
		return std::exp(-d*d*i_exp_fac);//*0.1*simVars.sim.h0;// + simVars.sim.h0;
#endif
	}


	double runTests()
	{
		sweet::ScalarDataArray out_data(posx_a.number_of_elements);

		if (interpolation_order == 2)
		{
			sphereDataSampler.bilinear_scalar(
					prog_h.toPhys(),
					posx_a,
					posy_a,
					out_data.scalar_data,
					false,
					use_poles_pseudo_points
			);
		}
		else if (interpolation_order == 3)
		{
			sphereDataSampler.bicubic_scalar(
					prog_h.toPhys(),
					posx_a,
					posy_a,
					out_data.scalar_data,
					false,
					use_poles_pseudo_points,
					use_limiter
			);
		}
		else
		{
			SWEETError("Interpolation order not available");
		}

		double max_error = 0;
		assert(posx_a.number_of_elements != 0);
		for (std::size_t i = 0; i < posx_a.number_of_elements; i++)
		{
			double value = gaussianValue(posx_a.scalar_data[i], posy_a.scalar_data[i], exp_fac);
			double err = std::abs(value - out_data.scalar_data[i]);
			max_error = std::max(max_error, err);
		}
		return max_error;
	}


	bool should_quit()
	{
		return false;
	}
};


int main(int i_argc, char *i_argv[])
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	int initial_spectral_modes = shackSphereDataOps->space_res_spectral[0];

	for (
			int use_poles_pseudo_points = 0;
			use_poles_pseudo_points < 2;
			use_poles_pseudo_points++
	)
	{
		std::cout << std::endl;
		std::cout << "*********************************************************" << std::endl;
		std::cout << "* Running studies with or without pseudo pole points " << use_poles_pseudo_points << std::endl;
		std::cout << "*********************************************************" << std::endl;

		//int i_gaussians_start = -1;

		for (int gaussian_id = -1; gaussian_id < 9; gaussian_id++)
		//for (int i_gaussians = 8; i_gaussians >= -1; i_gaussians--)
		{
			std::cout << std::endl;
			std::cout << "*********************************************************" << std::endl;
			std::cout << "* Running studies for Gaussian type " << gaussian_id << std::endl;
			std::cout << "*********************************************************" << std::endl;

			for (int interpolation_order = 2; interpolation_order <= 3; interpolation_order++)
			{
				std::cout << std::endl;
				std::cout << "*********************************************************" << std::endl;
				std::cout << "* Running studies for interpolation of order " << interpolation_order << std::endl;
				std::cout << "*********************************************************" << std::endl;

				int oversampling = 5;
				std::cout << "Using oversampling of " << oversampling << std::endl;

				double prev_max_error = -1;
				for (int i = initial_spectral_modes; i <= 256*2; i *= 2)
				{

					shackSphereDataOps->space_res_physical[0] = 2*i;
					shackSphereDataOps->space_res_physical[1] = i;

					shackSphereDataOps->space_res_spectral[0] = i;
					shackSphereDataOps->space_res_spectral[1] = i;

					sweet::SphereData_Config sphereDataConfig;
					sphereDataConfig.setupAuto(shackSphereDataOps);

					/*
					 * Generate higher resolution
					 */
					sweet::ShackSphereDataOps shackSphereDataOpsOversampled;
					shackSphereDataOpsOversampled = *shackSphereDataOps;

					shackSphereDataOpsOversampled.space_res_physical[0] *= oversampling;
					shackSphereDataOpsOversampled.space_res_physical[1] *= oversampling;

					shackSphereDataOpsOversampled.space_res_spectral[0] *= oversampling;
					shackSphereDataOpsOversampled.space_res_spectral[1] *= oversampling;

					sweet::SphereData_Config sphereDataConfigOversampled;
					sphereDataConfigOversampled.setupAuto(shackSphereDataOpsOversampled);


					{
						InterpolationTests interpolationTests;

						// Update interpolation order
						interpolationTests.interpolation_order = interpolation_order;

						// center of Gaussian bump
						interpolationTests.gaussian_id = gaussian_id;

						interpolationTests.use_poles_pseudo_points = use_poles_pseudo_points;

						interpolationTests.setup(sphereDataConfig, sphereDataConfigOversampled);

						double max_error = interpolationTests.runTests();
						{
							std::cout << "Lmax error: " << max_error << std::endl;

							if (prev_max_error >= 0)
							{
								//double conv = (prev_max_error - simulation.max_error) / simulation.max_error;
								double conv = prev_max_error / max_error;
								std::cout << "Convergence: " << conv << std::endl;

								if (use_poles_pseudo_points == 1 && gaussian_id == -1)
								{
									if (conv*1.2 < std::pow(2.0, interpolation_order))
										SWEETError("Convergence not given!");
								}
								else
								{
									if (conv*1.1 < std::pow(2.0, interpolation_order))
										SWEETError("Convergence not given!");
								}
							}
							prev_max_error = max_error;
						}
					}
				}
			}
		}
	}

	return 0;
}
