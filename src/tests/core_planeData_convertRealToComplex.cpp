/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *      
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <sweet/core/defaultPrecompilerValues.hpp>

#if SWEET_GUI
#	error	"GUI not supported"
#endif


#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/ProgramArguments.hpp>


#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneData_PhysicalComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>
#include <sweet/core/plane/Convert_PlaneDataPhysical_to_PlaneDataPhysicalComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataPhysicalComplex_to_PlaneDataPhysical.hpp>


class TestPlaneDataModes
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

		sweet::PlaneData_Config planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h;


		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			planeDataConfig.setupAuto(*i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(ops);

			prog_h.setup(planeDataConfig);

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

public:
	TestPlaneDataModes(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr)
	{
		ERROR_CHECK_WITH_RETURN(shackProgArgDict);
	}


	bool setup()
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

		shackProgArgDict.printShackData();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		data.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(data);

		return true;
	}


	void clear()
	{
		data.clear();

		shackPlaneDataOps = nullptr;
		shackProgArgDict.clear();
	}

	void test_planedata_planedatacomplex_physicalphysical_convert()
	{
		sweet::PlaneData_Physical test(data.planeDataConfig);
		sweet::PlaneData_PhysicalComplex testcplx(data.planeDataConfig);
		sweet::PlaneData_Spectral test_spc(data.planeDataConfig);
		sweet::PlaneData_SpectralComplex test_spc_cplx(data.planeDataConfig);

		for (std::size_t y = 0; y < data.planeDataConfig.physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test.physical_set_zero();
				test.physical_set_value(y, x, 1.0);

				testcplx = sweet::Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);

				test_spc_cplx.loadPlaneDataPhysical(testcplx);

				test_spc_cplx.test_realphysical();
				sweet::PlaneData_Physical tmp = sweet::Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

				double error = (test-tmp).physical_reduce_max_abs();

				if (error > 1e-8)
				{
					test_spc.loadPlaneDataPhysical(test);
					sweet::PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
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


		for (std::size_t y = 0; y < data.planeDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_physicalphysical_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test_spc.spectral_set_zero();
					test_spc.spectral_set(y, x, 1.0);
					test_spc.spectral_zeroAliasingModes();

					test = test_spc.toPhys();

					testcplx = sweet::Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::PlaneData_Physical tmp = sweet::Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						sweet::PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
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

					testcplx = sweet::Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::PlaneData_Physical tmp = sweet::Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						sweet::PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
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

					testcplx = sweet::Convert_PlaneDataPhysical_To_PlaneDataPhysicalComplex::physical_convert(test);
					test_spc_cplx.loadPlaneDataPhysical(testcplx);
					test_spc_cplx.test_realphysical();

					sweet::PlaneData_Physical tmp = sweet::Convert_PlaneDataPhysicalComplex_To_PlaneDataPhysical::physical_convert_real(testcplx);

					double error = (test-tmp).physical_reduce_max_abs();

					if (error > 1e-8)
					{
						test_spc.loadPlaneDataPhysical(test);
						sweet::PlaneData_Spectral tmp_spc(tmp.planeDataConfig);
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
		sweet::PlaneData_Spectral test(data.planeDataConfig);
		sweet::PlaneData_SpectralComplex testcplx(data.planeDataConfig);
		sweet::PlaneData_Physical test_phys(data.planeDataConfig);
		sweet::PlaneData_PhysicalComplex test_phys_cplx(data.planeDataConfig);

		for (std::size_t y = 0; y < data.planeDataConfig.physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.physical_set_zero();
				test_phys.physical_set_value(y, x, 1.0);

				test.loadPlaneDataPhysical(test_phys);

				//testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx.test_realphysical();

				//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
				sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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


		for (std::size_t y = 0; y < data.planeDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_physicalspectral_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					test.loadPlaneDataPhysical(test.toPhys()); // EXTRA DEALIASING REQUIRED FOR CORRECT RESULT??

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
					testcplx.test_realphysical();

					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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
		sweet::PlaneData_Spectral test(data.planeDataConfig);
		sweet::PlaneData_SpectralComplex testcplx(data.planeDataConfig);
		sweet::PlaneData_Physical test_phys(data.planeDataConfig);
		sweet::PlaneData_PhysicalComplex test_phys_cplx(data.planeDataConfig);


		for (std::size_t y = 0; y < data.planeDataConfig.physical_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.physical_data_size[0]; x++)
			{
				std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing physical 1.0 value at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
				test_phys.physical_set_zero();
				test_phys.physical_set_value(y, x, 1.0);
				test.loadPlaneDataPhysical(test_phys);

				//testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::physical_convert(test);
				testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
				testcplx.test_realphysical();

				//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
				sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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


		for (std::size_t y = 0; y < data.planeDataConfig.spectral_data_size[1]; y++)
		{
			for (std::size_t x = 0; x < data.planeDataConfig.spectral_data_size[0]; x++)
			{
				{
					std::cout << "test_planedata_planedatacomplex_spectralphysical_convert: Testing spectral value 1.0+0.0i at (" << y << ", " << x << ") with physical/physical conversion" << std::endl;
					test.spectral_set_zero();
					test.spectral_set(y, x, 1.0);
					test.spectral_zeroAliasingModes();

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);
					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);

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

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);
					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_imag(testcplx);

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

					testcplx = sweet::Convert_PlaneDataSpectral_To_PlaneDataSpectralComplex::spectral_convert(test);
					testcplx.test_realphysical();

					sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::spectral_convert_physical_real_only(testcplx);
					//sweet::PlaneData_Spectral tmp = sweet::Convert_PlaneDataSpectralComplex_To_PlaneDataSpectral::physical_convert_real(testcplx);

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

	void run_tests()
	{
		test_planedata_planedatacomplex_physicalphysical_convert();
		test_planedata_planedatacomplex_physicalspectral_convert();
		test_planedata_planedatacomplex_spectralphysical_convert();
	}
};



int main(int i_argc, char *i_argv[])
{
	TestPlaneDataModes simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.run_tests();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.clear();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);



	std::cout << "FIN" << std::endl;
	return 0;
}
