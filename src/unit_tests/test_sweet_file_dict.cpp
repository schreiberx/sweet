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
		fd.get("str", value);

		if (value != "hello world")
			SWEETError("string mismatch");

		std::cout << " + str: OK - " << value << std::endl;
	}

	{
		long long value;
		fd.get("int64", value);

		if (value != 123)
			SWEETError("int64 mismatch");

		std::cout << " + int64: OK - " << value << std::endl;
	}

	{
		double value;
		fd.get("float64", value);

		if (value != 123.456)
			SWEETError("float64 mismatch");

		std::cout << " + float64: OK - " << value << std::endl;
	}

	{
		SWEETArray<1,SWEETFileDict::float64> value;
		fd.get("array_float64_1d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			SWEETFileDict::float64 test_value = 1 + i0;
			if (value(i0) != test_value)
			{
				std::ostringstream ss;
				ss << "array_float64_1d: value (" << i0 << ") mismatch: " << value(i0) << " vs. " << test_value;
				SWEETError(ss.str());
			}
		}

		std::cout << " + array_1d_float64: OK - " << value << std::endl;
	}


	{
		SWEETArray<2,SWEETFileDict::float64> value;
		fd.get("array_float64_2d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			for (int i1 = 0; i1 < value.shape()[1]; i1++)
			{
				SWEETFileDict::float64 test_value = 1 + i0*value.shape()[1] + i1;
				if (value(i0, i1) != test_value)
				{
					std::ostringstream ss;
					ss << "array_float64_2d: value (" << i0 << ", " << i1 << ") mismatch: " << value(i0, i1) << " vs. " << test_value;
					SWEETError(ss.str());
				}
			}
		}

		std::cout << " + array_2d_float64: OK - " << value << std::endl;
	}


	{
		SWEETArray<3,SWEETFileDict::float64> value;
		fd.get("array_float64_3d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			for (int i1 = 0; i1 < value.shape()[1]; i1++)
			{
				for (int i2 = 0; i2 < value.shape()[2]; i2++)
				{
					SWEETFileDict::float64 test_value = 1 + i0*value.shape()[1]*value.shape()[2] + i1*value.shape()[2] + i2;
					if (value(i0, i1, i2) != test_value)
					{
						std::ostringstream ss;
						ss << "array_float64_3d: value (" << i0 << ", " << i1 <<  ", " << i2 << ") mismatch: " << value(i0, i1, i2) << " vs. " << test_value;
						SWEETError(ss.str());
					}
				}
			}
		}

		std::cout << " + array_3d_float64: OK - " << value << std::endl;
	}



	{
		SWEETArray<1,SWEETFileDict::complex128> value;
		fd.get("array_complex128_1d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			SWEETFileDict::complex128 test_value = 1 + i0;
			test_value += SWEETFileDict::complex128(0, 101+test_value.real());
			if (value(i0) != test_value)
			{
				std::ostringstream ss;
				ss << "array_complex128_1d: value (" << i0 << ") mismatch: " << value(i0) << " vs. " << test_value;
				SWEETError(ss.str());
			}
		}

		std::cout << " + array_1d_complex128: OK - " << value << std::endl;
	}


	{
		SWEETArray<2,SWEETFileDict::complex128> value;
		fd.get("array_complex128_2d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			for (int i1 = 0; i1 < value.shape()[1]; i1++)
			{
				SWEETFileDict::complex128 test_value = 1 + i0*value.shape()[1] + i1;
				test_value += SWEETFileDict::complex128(0, 101+test_value.real());
				if (value(i0, i1) != test_value)
				{
					std::ostringstream ss;
					ss << "array_complex128_2d: value (" << i0 << ", " << i1 << ") mismatch: " << value(i0, i1) << " vs. " << test_value;
					SWEETError(ss.str());
				}
			}
		}

		std::cout << " + array_2d_complex128: OK - " << value << std::endl;
	}


	{
		SWEETArray<3,SWEETFileDict::complex128> value;
		fd.get("array_complex128_3d", value);

		for (int i0 = 0; i0 < value.shape()[0]; i0++)
		{
			for (int i1 = 0; i1 < value.shape()[1]; i1++)
			{
				for (int i2 = 0; i2 < value.shape()[2]; i2++)
				{
					SWEETFileDict::complex128 test_value = 1 + i0*value.shape()[1]*value.shape()[2] + i1*value.shape()[2] + i2;
					test_value += SWEETFileDict::complex128(0, 101+test_value.real());
					if (value(i0, i1, i2) != test_value)
					{
						std::ostringstream ss;
						ss << "array_complex128_3d: value (" << i0 << ", " << i1 <<  ", " << i2 << ") mismatch: " << value(i0, i1, i2) << " vs. " << test_value;
						SWEETError(ss.str());
					}
				}
			}
		}

		std::cout << " + array_3d_complex128: OK - " << value << std::endl;
	}
}


void run_sweet_array_assignment_test()
{
	std::cout << " + run_sweet_array_assignment_test()" << std::endl;
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
			SWEETFileDict::float64 real_value = test_array(i0, i1);
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


void _setup_dict_test_data(
		int size_i0, int size_i1, int size_i2,
		int version_id,
		SWEETFileDict &io_fileDict
)
{
	io_fileDict.set("str", "hello world");

	io_fileDict.set("int64", 123);
	io_fileDict.set("float64", 123.456);

	std::array<int,1> shape1d({size_i2});
	std::array<int,2> shape2d({size_i1, size_i2});
	std::array<int,3> shape3d({size_i0, size_i1, size_i2});

	int size1d = size_i2;
	int size2d = size_i1*size_i2;
	int size3d = size_i0*size_i1*size_i2;

	{
		SWEETFileDict::float64 data[3];
		for (int i = 0; i < size1d; i++)
			data[i] = i+1;

		if (0)
		{
			SWEETArray<1,SWEETFileDict::float64> array(shape1d, data);

			io_fileDict.set("array_float64_1d", array);
		}

		if (1)
		{
			SWEETArray<1,SWEETFileDict::float64> array;
			array.setup(shape1d);
			array = data;

			io_fileDict.set("array_float64_1d", array);
		}

		if (2)
		{
			SWEETArray<1,SWEETFileDict::float64> array;
			array.setup(shape1d);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
			{
				array.set(i0, data[c]);
				c++;
			}

			io_fileDict.set("array_float64_1d", array);
		}
	}

	{
		// io_fileDict.set("array_float64_2d", np.array([[1., 2., 3.], [4., 5., 6.]]));

		std::array<int,2> shape(shape2d);
		SWEETFileDict::float64 data[size2d];

		for (int i = 0; i < size2d; i++)
			data[i] = i+1;

		if (0)
		{
			SWEETArray<2,SWEETFileDict::float64> array(shape, data);

			io_fileDict.set("array_float64_2d", array);
		}

		if (1)
		{
			SWEETArray<2,SWEETFileDict::float64> array;
			array.setup(shape);
			array = data;

			io_fileDict.set("array_float64_2d", array);
		}

		if (2)
		{
			SWEETArray<2,SWEETFileDict::float64> array;
			array.setup(shape);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
				for (int i1 = 0; i1 < array.shape()[1]; i1++)
				{
					array.set(i0, i1, data[c]);
					c++;
				}

			io_fileDict.set("array_float64_2d", array);
		}
	}

	{
		SWEETFileDict::float64 data[size3d];
		for (int i = 0; i < sizeof(data)/sizeof(SWEETFileDict::float64); i++)
			data[i] = i+1;

		if (0)
		{
			SWEETArray<3,SWEETFileDict::float64> array(shape3d, data);

			io_fileDict.set("array_float64_3d", array);
		}

		if (1)
		{
			SWEETArray<3,SWEETFileDict::float64> array;
			array.setup(shape3d);
			array = data;

			io_fileDict.set("array_float64_3d", array);
		}

		if (2)
		{
			SWEETArray<3,SWEETFileDict::float64> array;
			array.setup(shape3d);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
				for (int i1 = 0; i1 < array.shape()[1]; i1++)
					for (int i2 = 0; i2 < array.shape()[2]; i2++)
					{
						array.set3(i0, i1, i2, data[c]);
						c++;
					}

			io_fileDict.set("array_float64_3d", array);
		}
	}


	{
		SWEETFileDict::complex128 data[size1d];
		for (int i = 0; i < size1d; i++)
			data[i] = SWEETFileDict::complex128(i+1, 100+i);

		if (0)
		{
			SWEETArray<1,SWEETFileDict::complex128> array(shape1d, data);

			io_fileDict.set("array_complex128_1d", array);
		}

		if (1)
		{
			SWEETArray<1,SWEETFileDict::complex128> array;
			array.setup(shape1d);
			array = data;

			io_fileDict.set("array_complex128_1d", array);
		}

		if (2)
		{
			SWEETArray<1,SWEETFileDict::complex128> array;
			array.setup(shape1d);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
			{
				array.set(i0, data[c]);
				c++;
			}

			io_fileDict.set("array_complex128_1d", array);
		}
	}

	{
		// io_fileDict.set("array_complex128_2d", np.array([[1.+102j, 2.+103j, 3.+104j], [4.+105j, 5.+106j, 6.+107j]]))

		SWEETFileDict::complex128 data[size2d];

		for (int i = 0; i < size2d; i++)
			data[i] = SWEETFileDict::complex128(i+1, 100+i);

		if (0)
		{
			SWEETArray<2,SWEETFileDict::complex128> array(shape2d, data);

			io_fileDict.set("array_complex128_2d", array);
		}

		if (1)
		{
			SWEETArray<2,SWEETFileDict::complex128> array;
			array.setup(shape2d);
			array = data;

			io_fileDict.set("array_complex128_2d", array);
		}

		if (2)
		{
			SWEETArray<2,SWEETFileDict::complex128> array;
			array.setup(shape2d);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
				for (int i1 = 0; i1 < array.shape()[1]; i1++)
				{
					array.set(i0, i1, data[c]);
					c++;
				}

			io_fileDict.set("array_complex128_2d", array);
		}
	}

	{
		// io_fileDict.set("array_complex128_3d", np.array([[[1.+102j, 2.+103j, 3.+104j], [4.+105j, 5.+106j, 6.+107j]], [[7.+108j, 8.+109j, 9.+110j], [10.+111j, 11.+112j, 12.+113j]]]))

		SWEETFileDict::complex128 data[size3d];
		for (int i = 0; i < size3d; i++)
			data[i] = SWEETFileDict::complex128(i+1, 100+i);

		if (0)
		{
			SWEETArray<3,SWEETFileDict::complex128> array(shape3d, data);

			io_fileDict.set("array_complex128_3d", array);
		}

		if (1)
		{
			SWEETArray<3,SWEETFileDict::complex128> array;
			array.setup(shape3d);
			array = data;

			io_fileDict.set("array_complex128_3d", array);
		}

		if (2)
		{
			SWEETArray<3,SWEETFileDict::complex128> array;
			array.setup(shape3d);

			int c = 0;
			for (int i0 = 0; i0 < array.shape()[0]; i0++)
				for (int i1 = 0; i1 < array.shape()[1]; i1++)
					for (int i2 = 0; i2 < array.shape()[2]; i2++)
					{
						array.set3(i0, i1, i2, data[c]);
						c++;
					}

			io_fileDict.set("array_complex128_3d", array);
		}
	}
}

void run_sweet_file_dict_creation_test()
{
	std::cout << " + run_sweet_file_dict_creation_test()" << std::endl;

	SWEETFileDict fileDict;


	for (int version_id = 0; version_id < 3; version_id++)
	{

		_setup_dict_test_data(2, 2, 3, version_id, fileDict);


		std::cout << fileDict << std::endl;
	}
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

	run_sweet_array_assignment_test();

	run_sweet_file_dict_creation_test();

	std::cout << "FIN" << std::endl;

	return EXIT_SUCCESS;
}
