#include <iostream>
#include <stdlib.h>
#include <sweet/SWEETFileDict.hpp>
#include <sweet/SimulationVariables.hpp>


int main(
		int i_argc,
		char *const i_argv[]
)
{
	std::string filename;

	if (i_argc >= 2)
	{
		filename = i_argv[1];
	}
	else
	{
		filename = "output_shared_dict.sweet";

		std::cout << "Assuming that '" << filename << "' should be processed." << std::endl;
		std::cout << "Override with 1st program parameter" << std::endl;
	}

	SimulationVariables simVars;

	const char *user_defined_prog_params[] = {{"sweet-file-dict"}, {"unit-test-sweet-file-dict"}, nullptr};
	if (!simVars.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params, false))
		return EXIT_FAILURE;

	if (simVars.user_defined.var[0] == "")
	{
		std::cerr << "Please specify sweet file dictionary with --sweet-file-dict=[filename]" << std::endl;
		SWEETError("EXIT");
	}

	SWEETFileDict fd(simVars.user_defined.var[0], true);

	std::cout << fd << std::endl;

	if (simVars.user_defined.var[1] != "")
	{
		/*
		 * Run unit tests
		 */

		std::cout << "Validating:" << std::endl;

		{
			std::string value;
			fd.getValue("str", value);

			if (value != "hello world")
				SWEETError("string mismatch");

			std::cout << " + str: OK - " << value << std::endl;
		}

		{
			long long value;
			fd.getValue("int64", value);

			if (value != 123)
				SWEETError("int64 mismatch");

			std::cout << " + int64: OK - " << value << std::endl;
		}

		{
			double value;
			fd.getValue("float64", value);

			if (value != 123.456)
				SWEETError("float64 mismatch");

			std::cout << " + float64: OK - " << value << std::endl;
		}

		{
			SWEETArray<1,SWEETFileDict::float64> value;
			fd.getValue("array_float64_1d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				SWEETFileDict::float64 test_value = 1 + i0;
				if (value.get(i0) != test_value)
				{
					std::ostringstream ss;
					ss << "array_float64_1d: value (" << i0 << ") mismatch: " << value.get(i0) << " vs. " << test_value;
					SWEETError(ss.str());
				}
			}

			std::cout << " + array_1d_float64: OK - " << value << std::endl;
		}


		{
			SWEETArray<2,SWEETFileDict::float64> value;
			fd.getValue("array_float64_2d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				for (int i1 = 0; i1 < value.shape()[1]; i1++)
				{
					SWEETFileDict::float64 test_value = 1 + i0*value.shape()[1] + i1;
					if (value.get(i0, i1) != test_value)
					{
						std::ostringstream ss;
						ss << "array_float64_2d: value (" << i0 << ", " << i1 << ") mismatch: " << value.get(i0, i1) << " vs. " << test_value;
						SWEETError(ss.str());
					}
				}
			}

			std::cout << " + array_2d_float64: OK - " << value << std::endl;
		}


		{
			SWEETArray<3,SWEETFileDict::float64> value;
			fd.getValue("array_float64_3d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				for (int i1 = 0; i1 < value.shape()[1]; i1++)
				{
					for (int i2 = 0; i2 < value.shape()[2]; i2++)
					{
						SWEETFileDict::float64 test_value = 1 + i0*value.shape()[1]*value.shape()[2] + i1*value.shape()[2] + i2;
						if (value.get(i0, i1, i2) != test_value)
						{
							std::ostringstream ss;
							ss << "array_float64_3d: value (" << i0 << ", " << i1 <<  ", " << i2 << ") mismatch: " << value.get(i0, i1, i2) << " vs. " << test_value;
							SWEETError(ss.str());
						}
					}
				}
			}

			std::cout << " + array_3d_float64: OK - " << value << std::endl;
		}



		{
			SWEETArray<1,SWEETFileDict::complex128> value;
			fd.getValue("array_complex128_1d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				SWEETFileDict::complex128 test_value = 1 + i0;
				test_value += SWEETFileDict::complex128(0, 101+test_value.real());
				if (value.get(i0) != test_value)
				{
					std::ostringstream ss;
					ss << "array_complex128_1d: value (" << i0 << ") mismatch: " << value.get(i0) << " vs. " << test_value;
					SWEETError(ss.str());
				}
			}

			std::cout << " + array_1d_complex128: OK - " << value << std::endl;
		}


		{
			SWEETArray<2,SWEETFileDict::complex128> value;
			fd.getValue("array_complex128_2d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				for (int i1 = 0; i1 < value.shape()[1]; i1++)
				{
					SWEETFileDict::complex128 test_value = 1 + i0*value.shape()[1] + i1;
					test_value += SWEETFileDict::complex128(0, 101+test_value.real());
					if (value.get(i0, i1) != test_value)
					{
						std::ostringstream ss;
						ss << "array_complex128_2d: value (" << i0 << ", " << i1 << ") mismatch: " << value.get(i0, i1) << " vs. " << test_value;
						SWEETError(ss.str());
					}
				}
			}

			std::cout << " + array_2d_complex128: OK - " << value << std::endl;
		}


		{
			SWEETArray<3,SWEETFileDict::complex128> value;
			fd.getValue("array_complex128_3d", value);

			for (int i0 = 0; i0 < value.shape()[0]; i0++)
			{
				for (int i1 = 0; i1 < value.shape()[1]; i1++)
				{
					for (int i2 = 0; i2 < value.shape()[2]; i2++)
					{
						SWEETFileDict::complex128 test_value = 1 + i0*value.shape()[1]*value.shape()[2] + i1*value.shape()[2] + i2;
						test_value += SWEETFileDict::complex128(0, 101+test_value.real());
						if (value.get(i0, i1, i2) != test_value)
						{
							std::ostringstream ss;
							ss << "array_complex128_3d: value (" << i0 << ", " << i1 <<  ", " << i2 << ") mismatch: " << value.get(i0, i1, i2) << " vs. " << test_value;
							SWEETError(ss.str());
						}
					}
				}
			}

			std::cout << " + array_3d_complex128: OK - " << value << std::endl;
		}


	}

	return EXIT_SUCCESS;
}
