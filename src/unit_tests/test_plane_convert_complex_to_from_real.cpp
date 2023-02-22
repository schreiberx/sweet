/*
 * test_planedata_convert_complex_to_from_real.cpp
 *
 *  Created on: 18 July 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>
#include <sweet/plane/Convert_PlaneDataPhysical_to_PlaneDataPhysicalComplex.hpp>
#include <sweet/plane/Convert_PlaneDataPhysicalComplex_to_PlaneDataPhysical.hpp>

SimulationVariables simVars;

PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

class PlaneDataModes
{
public:
	PlaneOperators op;

	void test_planedata_planedatacomplex_physicalphysical_convert()
	{
		PlaneData_Physical test(planeDataConfig);
		PlaneData_PhysicalComplex testcplx(planeDataConfig);
		PlaneData_Spectral test_spc(planeDataConfig);
		PlaneData_SpectralComplex test_spc_cplx(planeDataConfig);

		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.physical_set_zero();
				test.physical_set_value(y, x, 1.0);

				testcplx = Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);

				test_spc_cplx.loadPlaneDataPhysical(testcplx);

				test_spc_cplx.test_realphysical();
				PlaneData_Physical tmp = Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

				double error = (test-tmp).physical_reduce_max_abs();

				if (error > 1e-8)
				{
					test_spc.loadPlaneDataPhysical(test);
					PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
					tmp_spc.loadPlaneDataPhysical(tmp);

					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test_spc.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					test_spc_cplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp_spc.print_spectralData_zeroNumZero();

					SWEETError("Inconsistency detected a");
				}
			}
		}


		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_set_zero();
					test_spc.spectral_set(y, x, 1.0);
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toPhys();

					testcplx = Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					PlaneData_Physical tmp = Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
						tmp_spc.loadPlaneDataPhysical(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected b");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_set_zero();
					test_spc.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toPhys();

					testcplx = Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					PlaneData_Physical tmp = Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
						tmp_spc.loadPlaneDataPhysical(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected c");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_set_zero();
					test_spc.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toPhys();

					testcplx = Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					PlaneData_Physical tmp = Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
						tmp_spc.loadPlaneDataPhysical(tmp);

						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test_spc.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						test_spc_cplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp_spc.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected d");
					}
				}
			}
		}
	}


	void test_planedata_planedatacomplex_physicalspectral_convert()
	{
		PlaneData_Spectral test(planeDataConfig);
		PlaneData_SpectralComplex testcplx(planeDataConfig);
		PlaneData_Physical test_phys(planeDataConfig);
		PlaneData_PhysicalComplex test_phys_cplx(planeDataConfig);

		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.physical_set_zero();
				test_phys.physical_set_value(y, x, 1.0);

				test.loadPlaneDataPhysical(test_phys);

				//testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx.test_realphysical();

				//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
				PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

				double error = (test-tmp).spectral_reduce_max_abs();

				if (error > 1e-8)
				{
					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					testcplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp.print_spectralData_zeroNumZero();

					SWEETError("Inconsistency detected e");
				}
			}
		}


		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					test.loadPlaneDataPhysical(test.toPhys()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected f");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test.spectral_zeroAliasingModes();

					test.loadPlaneDataPhysical(test.toPhys()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected g");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test.spectral_zeroAliasingModes();

					test.loadPlaneDataPhysical(test.toPhys()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected h");
					}
				}
			}
		}
	}

	void test_planedata_planedatacomplex_spectralphysical_convert()
	{
		PlaneData_Spectral test(planeDataConfig);
		PlaneData_SpectralComplex testcplx(planeDataConfig);
		PlaneData_Physical test_phys(planeDataConfig);
		PlaneData_PhysicalComplex test_phys_cplx(planeDataConfig);


		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.physical_set_zero();
				test_phys.physical_set_value(y, x, 1.0);
				test.loadPlaneDataPhysical(test_phys);

				//testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
				testcplx.test_realphysical();

				//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
				PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

				double error = (test-tmp).spectral_reduce_max_abs();

				if (error > 1e-8)
				{
					std::cout << std::endl;
					std::cout << "Original" << std::endl;
					test.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original complex" << std::endl;
					testcplx.print_spectralData_zeroNumZero();

					std::cout << std::endl;
					std::cout << "Original converted" << std::endl;
					tmp.print_spectralData_zeroNumZero();

					SWEETError("Inconsistency detected i");
				}
			}
		}


		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected j");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, std::complex<double>(0.0, 1.0));
					test.spectral_zeroAliasingModes();

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);
					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_imag(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected k");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, std::complex<double>(1.0, 1.0));
					test.spectral_zeroAliasingModes();

					testcplx = Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);
					//PlaneData_Spectral tmp = Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);

					double error = (test-tmp).spectral_reduce_max_abs();

					if (error > 1e-8)
					{
						std::cout << std::endl;
						std::cout << "Original" << std::endl;
						test.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original complex" << std::endl;
						testcplx.print_spectralData_zeroNumZero();

						std::cout << std::endl;
						std::cout << "Original converted" << std::endl;
						tmp.print_spectralData_zeroNumZero();

						SWEETError("Inconsistency detected l");
					}
				}
			}
		}
	}

	PlaneDataModes()	:
		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		test_planedata_planedatacomplex_physicalphysical_convert();
		test_planedata_planedatacomplex_physicalspectral_convert();
		test_planedata_planedatacomplex_spectralphysical_convert();
	}
};



int main(
		int i_argc,
		char *i_argv[]
)
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	simVars.outputConfig();

	PlaneDataModes planeDataModes;

	return 0;
}
