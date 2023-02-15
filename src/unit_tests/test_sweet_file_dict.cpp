#include <string>
#include <iostream>
#include <stdlib.h>
#include <sweet/SWEETFileDict.hpp>
#include <sweet/SimulationVariables.hpp>


/*
 * Read the file given by the parameter and compare it to predefined values
 */
void run_compare_file_dict(
		const SWEETFileDict &fd
)
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


void run_assignment_test()
{
	SWEETArray<2,double> test_array;

#if 1
	test_array.setup({2, 3});
#else
	const double shape[] = {2, 3};
	test_array.setup(shape);
#endif


	const double tmp[] = {
			1, 2, 3,
			4, 5, 6
	};

	test_array = tmp;

	/*
	 * Validate
	 */
	for (int i0 = 0; i0 < test_array.shape()[0]; i0++)
	{
		for (int i1 = 0; i1 < test_array.shape()[1]; i1++)
		{
			SWEETFileDict::float64 test_value = 1 + i0*test_array.shape()[1] + i1;
			SWEETFileDict::float64 real_value = test_array.get(i0, i1);
			if (real_value != test_value)
			{
				std::ostringstream ss;
				ss << "flat array assignment error: value (" << i0 << ", " << i1 << ") mismatch: " << real_value << " vs. " << test_value;
				SWEETError(ss.str());
			}
		}
	}

	std::cout << " + flat array assignment: OK - " << test_array << std::endl;

	std::cout << test_array << std::endl;
}


int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;

	// Do not use this, since it's not supported by, e.g., gcc-8 compiler
	// const char *user_defined_prog_params[] = {{"sweet-file-dict"}, {"run-unit-test-sweet-file-dict"}, nullptr};

	const char *user_defined_prog_params[3];
	user_defined_prog_params[0] = "sweet-file-dict";
	user_defined_prog_params[1] = "run-unit-test-sweet-file-dict";
	user_defined_prog_params[2] = nullptr;

	if (!simVars.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params, false))
		return EXIT_FAILURE;

	/*
	 * Convert program parameters
	 */
	const std::string &sweet_file_dict_filename = simVars.user_defined.var[0];
	bool run_unit_test = atoi(simVars.user_defined.var[1].c_str());


	if (sweet_file_dict_filename != "")
	{
		/*
		 * Simply read the file dictionary
		 */
		SWEETFileDict fd(sweet_file_dict_filename, true);
		std::cout << fd << std::endl;

		/*
		 * Should we also run tests on it?
		 */
		if (run_unit_test)
		{
			run_compare_file_dict(sweet_file_dict_filename);
		}
	}

	run_assignment_test();

	return EXIT_SUCCESS;
}
