/*
 * output_plane_data_modes.cpp
 *
 *  Created on: 18 Juli 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

SimulationVariables simVars;

PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

class PlaneDataModes
{
public:
	PlaneOperators op;

	void test_planedata_planedatacomplex_physicalphysical_convert()
	{
		PlaneData test(planeDataConfig);
		PlaneDataComplex testcplx(planeDataConfig);

		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
				testcplx.test_realphysical();

				PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

				double error = (test-tmp).reduce_maxAbs();

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

					FatalError("Inconsistency detected");
				}
			}
		}


		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, 1.0);

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(0.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}
			}
		}
	}


	void test_planedata_planedatacomplex_physicalspectral_convert()
	{
		PlaneData test(planeDataConfig);
		PlaneDataComplex testcplx(planeDataConfig);

		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
				testcplx.test_realphysical();

				PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(testcplx);

				double error = (test-tmp).reduce_maxAbs();

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

					FatalError("Inconsistency detected");
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
					test.p_spectral_set(y, x, 1.0);

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(0.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::physical_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::spectral_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}
			}
		}
	}

	void test_planedata_planedatacomplex_spectralphysical_convert()
	{
		PlaneData test(planeDataConfig);
		PlaneDataComplex testcplx(planeDataConfig);

		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				testcplx = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(test);
				testcplx.test_realphysical();

				PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

				double error = (test-tmp).reduce_maxAbs();

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

					FatalError("Inconsistency detected");
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
					test.p_spectral_set(y, x, 1.0);

					testcplx = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 0.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(0.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}

				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 1.0+1.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

					testcplx = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(test);
					testcplx.test_realphysical();

					PlaneData tmp = Convert_PlaneDataComplex_To_PlaneData::physical_convert(testcplx);

					double error = (test-tmp).reduce_maxAbs();

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

						FatalError("Inconsistency detected");
					}
				}
			}
		}
	}


	void output_planedata()
	{
		PlaneData test(planeDataConfig);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Phys coord : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalArrayData();
				std::cout << "***" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}
	}



	void output_planedatacomplex()
	{
		PlaneDataComplex test(planeDataConfig);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->physical_data_size[0]; x++)
			{
				test.physical_set_zero();
				test.p_physical_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Phys coord : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, 1.0);

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		for (std::size_t y = 0; y < planeDataConfig->spectral_complex_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < planeDataConfig->spectral_complex_data_size[0]; x++)
			{
				test.spectral_set_zero();
				test.p_spectral_set(y, x, std::complex<double>(1.0, 1.0));

				std::cout << "***************************************" << std::endl;
				std::cout << "* Modes : " << y << ", " << x << std::endl;
				std::cout << "***************************************" << std::endl;
				std::cout << "Spectral space:" << std::endl;
				test.print_spectralData_zeroNumZero();
				std::cout << "***" << std::endl;
				std::cout << "Physical space:" << std::endl;
				test.print_physicalData_zeroNumZero();
			}
		}
	}

	PlaneDataModes()	:
		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		if (simVars.pde.id == 0)
		{
			test_planedata_planedatacomplex_physicalphysical_convert();
			test_planedata_planedatacomplex_physicalspectral_convert();
			test_planedata_planedatacomplex_spectralphysical_convert();
		}
		else
		{
			output_planedata();
			output_planedatacomplex();
		}
	}
};



int main(
		int i_argc,
		char *i_argv[]
)
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	simVars.outputConfig();

	PlaneDataModes planeDataModes;

	return 0;
}
