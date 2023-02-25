#include <string>
#include <iostream>
#include <stdlib.h>
#include <sweet/core/dict/Dict.hpp>
#include <sweet/core/SimulationVariables.hpp>


template <int D>
void create_reference_array_float(
		sweet::DictArrayND<D,sweet::Dict::float64> &io_array,
		sweet::Dict::float64 i_offset = 0
)
{
	std::array<int,3> shape_iter;

	if (D == 1)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = 1;
		shape_iter[2] = 1;
	}
	else if (D == 2)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = io_array.shape()[1];
		shape_iter[2] = 1;
	}
	else if (D == 3)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = io_array.shape()[1];
		shape_iter[2] = io_array.shape()[2];
	}

	int c = 0;
	for (int i0 = 0; i0 < shape_iter[0]; i0++)
	{
		for (int i1 = 0; i1 < shape_iter[1]; i1++)
		{
			for (int i2 = 0; i2 < shape_iter[2]; i2++)
			{
				sweet::Dict::float64 test_value = i_offset + i0*shape_iter[1]*shape_iter[2] + i1*shape_iter[2] + i2;
				io_array.data()[c] = test_value;

				c++;
			}
		}
	}
}


template <int D>
void create_reference_array_complex(
		sweet::DictArrayND<D,sweet::Dict::complex128> &io_array,
		sweet::Dict::float64 i_offset = 0
)
{
	std::array<int,3> shape_iter;

	if (D == 1)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = 1;
		shape_iter[2] = 1;
	}
	else if (D == 2)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = io_array.shape()[1];
		shape_iter[2] = 1;
	}
	else if (D == 3)
	{
		shape_iter[0] = io_array.shape()[0];
		shape_iter[1] = io_array.shape()[1];
		shape_iter[2] = io_array.shape()[2];
	}

	int c = 0;
	for (int i0 = 0; i0 < shape_iter[0]; i0++)
	{
		for (int i1 = 0; i1 < shape_iter[1]; i1++)
		{
			for (int i2 = 0; i2 < shape_iter[2]; i2++)
			{
				sweet::Dict::complex128 test_value = sweet::Dict::complex128(i_offset + c, 101+ i_offset + c);
				io_array.data()[c] = test_value;
				c++;
			}
		}
	}
}

void create_reference_dict(
	sweet::Dict &io_dict,
	int i_offset = 0,
	int i0 = 2,
	int i1 = 2,
	int i2 = 3
)
{
	io_dict.set("str", "hello world");
	io_dict.set<sweet::Dict::int64>("int64", 123);
	io_dict.set<sweet::Dict::float64>("float64", 123.456);
	io_dict.set<sweet::Dict::complex128>("complex128", std::complex<double>(123.456, 123.456+100));

	sweet::DictArrayND<1,sweet::Dict::float64> array1dFloat64(i2);
	create_reference_array_float(array1dFloat64, 1);
	io_dict.set("array_float64_1d", array1dFloat64);

	sweet::DictArrayND<2,sweet::Dict::float64> array2dFloat64(i1, i2);
	create_reference_array_float(array2dFloat64, 1);
	io_dict.set("array_float64_2d", array2dFloat64);

	sweet::DictArrayND<3,sweet::Dict::float64> array3dFloat64(i0, i1, i2);
	create_reference_array_float(array3dFloat64, 1);
	io_dict.set("array_float64_3d", array3dFloat64);


	sweet::DictArrayND<1,sweet::Dict::complex128> array1dComplex64(i2);
	create_reference_array_complex(array1dComplex64, 1);
	io_dict.set("array_complex128_1d", array1dComplex64);

	sweet::DictArrayND<2,sweet::Dict::complex128> array2dComplex64(i1, i2);
	create_reference_array_complex(array2dComplex64, 1);
	io_dict.set("array_complex128_2d", array2dComplex64);

	sweet::DictArrayND<3,sweet::Dict::complex128> array3dComplex64(i0, i1, i2);
	create_reference_array_complex(array3dComplex64, 1);
	io_dict.set("array_complex128_3d", array3dComplex64);
}



void compareDict1EntriesEqualDict2(
		const sweet::Dict &i_dict1,
		const sweet::Dict &i_dict2,
		const std::string &i_prefix
)
{
	for (int i = 0; i < i_dict1.size(); i++)
	{
		const sweet::_DictElement &e1 = i_dict1[i];

		const sweet::_DictElementTypeEnum::Enum typeId = e1.getTypeID();

		std::cout <<i_prefix << "Comparing '" << e1.getKey() << "'" << std::endl;
		switch(typeId)
		{
			default:
				SWEETError("Unknown type");
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_STRING:
			{
				typedef std::string T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_INT64:
			{
				typedef sweet::Dict::int64 T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_FLOAT64:
			{
				typedef sweet::Dict::float64 T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_COMPLEX128:
			{
				typedef sweet::Dict::complex128 T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
			{
				typedef sweet::DictArrayND<1,sweet::Dict::float64> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
			{
				typedef sweet::DictArrayND<2,sweet::Dict::float64> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
			{
				typedef sweet::DictArrayND<3,sweet::Dict::float64> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
			{
				typedef sweet::DictArrayND<1,sweet::Dict::complex128> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
			{
				typedef sweet::DictArrayND<2,sweet::Dict::complex128> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;

			case sweet::_DictElementTypeEnum::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
			{
				typedef sweet::DictArrayND<3,sweet::Dict::complex128> T;

				T value1;	e1.get(value1);
				T value2;	i_dict2.get(e1.getKey(), value2);
				if (value1 != value2)
					SWEETError("Values don't match");
			}
				break;
		}
	}
}

/*
 * Read the file given by the parameter and compare it to predefined values
 */
void run_compare_file_dict_with_reference(
		const sweet::Dict &i_dict
)
{
	sweet::Dict ref_dict;
	create_reference_dict(ref_dict);

	std::cout << " + Comparing dictionaries: dict with ref" << std::endl;
	compareDict1EntriesEqualDict2(i_dict, ref_dict, "   + ");
	std::cout << + "    + OK" << std::endl;

	std::cout << " + Comparing dictionaries: ref with dict" << std::endl;
	compareDict1EntriesEqualDict2(ref_dict, i_dict, "   + ");
	std::cout << + "    + OK" << std::endl;
}


void run_sweet_array_assignment_test()
{
	std::cout << " + run_sweet_array_assignment_test()" << std::endl;
	sweet::DictArrayND<2,double> test_array;

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
			sweet::Dict::float64 test_value = 1 + i0*test_array.shape()[1] + i1;
			sweet::Dict::float64 real_value = test_array.get(i0, i1);
			if (real_value != test_value)
			{
				std::ostringstream ss;
				ss << "flat array assignment error: value (" << i0 << ", " << i1 << ") mismatch: " << real_value << " vs. " << test_value;
				SWEETError(ss.str());
			}
		}
	}

	std::cout << " + flat array assignment: OK - " << test_array << std::endl;
}



/*
 * Test different variants to set the data in a dictionary
 */
void _setup_dict_test_data_different_variants(
		int size_i0, int size_i1, int size_i2,
		int version_id,
		sweet::Dict &io_fileDict
)
{
	io_fileDict.set("str", "hello world");

	io_fileDict.set("int64", 123);
	io_fileDict.set("float64", 123.456);
	io_fileDict.set("complex128", std::complex<double>(123.456, 123.456+100));

	std::array<int,1> shape1d({size_i2});
	std::array<int,2> shape2d({size_i1, size_i2});
	std::array<int,3> shape3d({size_i0, size_i1, size_i2});

	int size1d = size_i2;
	int size2d = size_i1*size_i2;
	int size3d = size_i0*size_i1*size_i2;

	{
		sweet::Dict::float64 data[3];
		for (int i = 0; i < size1d; i++)
			data[i] = i+1;

		if (version_id == 0)
		{
			sweet::DictArrayND<1,sweet::Dict::float64> array(shape1d, data);

			io_fileDict.set("array_float64_1d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<1,sweet::Dict::float64> array;
			array.setup(shape1d);
			array = data;

			io_fileDict.set("array_float64_1d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<1,sweet::Dict::float64> array;
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
		sweet::Dict::float64 data[size2d];

		for (int i = 0; i < size2d; i++)
			data[i] = i+1;

		if (version_id == 0)
		{
			sweet::DictArrayND<2,sweet::Dict::float64> array(shape, data);

			io_fileDict.set("array_float64_2d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<2,sweet::Dict::float64> array;
			array.setup(shape);
			array = data;

			io_fileDict.set("array_float64_2d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<2,sweet::Dict::float64> array;
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
		sweet::Dict::float64 data[size3d];
		for (std::size_t i = 0; i < sizeof(data)/sizeof(sweet::Dict::float64); i++)
			data[i] = i+1;

		if (version_id == 0)
		{
			sweet::DictArrayND<3,sweet::Dict::float64> array(shape3d, data);

			io_fileDict.set("array_float64_3d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<3,sweet::Dict::float64> array;
			array.setup(shape3d);
			array = data;

			io_fileDict.set("array_float64_3d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<3,sweet::Dict::float64> array;
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
		sweet::Dict::complex128 data[size1d];
		for (int i = 0; i < size1d; i++)
			data[i] = sweet::Dict::complex128(i+1, 100+i);

		if (version_id == 0)
		{
			sweet::DictArrayND<1,sweet::Dict::complex128> array(shape1d, data);

			io_fileDict.set("array_complex128_1d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<1,sweet::Dict::complex128> array;
			array.setup(shape1d);
			array = data;

			io_fileDict.set("array_complex128_1d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<1,sweet::Dict::complex128> array;
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

		sweet::Dict::complex128 data[size2d];

		for (int i = 0; i < size2d; i++)
			data[i] = sweet::Dict::complex128(i+1, 100+i);

		if (version_id == 0)
		{
			sweet::DictArrayND<2,sweet::Dict::complex128> array(shape2d, data);

			io_fileDict.set("array_complex128_2d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<2,sweet::Dict::complex128> array;
			array.setup(shape2d);
			array = data;

			io_fileDict.set("array_complex128_2d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<2,sweet::Dict::complex128> array;
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

		sweet::Dict::complex128 data[size3d];
		for (int i = 0; i < size3d; i++)
			data[i] = sweet::Dict::complex128(i+1, 100+i);

		if (version_id == 0)
		{
			sweet::DictArrayND<3,sweet::Dict::complex128> array(shape3d, data);

			io_fileDict.set("array_complex128_3d", array);
		}

		if (version_id == 1)
		{
			sweet::DictArrayND<3,sweet::Dict::complex128> array;
			array.setup(shape3d);
			array = data;

			io_fileDict.set("array_complex128_3d", array);
		}

		if (version_id == 2)
		{
			sweet::DictArrayND<3,sweet::Dict::complex128> array;
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

	sweet::Dict fileDict0;
	sweet::Dict fileDict1;
	sweet::Dict fileDict2;


	_setup_dict_test_data_different_variants(2, 2, 3, 0, fileDict0);
	_setup_dict_test_data_different_variants(2, 2, 3, 1, fileDict1);
	_setup_dict_test_data_different_variants(2, 2, 3, 2, fileDict2);

	std::cout << "   + Comparing dictionaries: dict0 with dict1" << std::endl;
	compareDict1EntriesEqualDict2(fileDict0, fileDict1, "   + ");
	std::cout << "   + OK" << std::endl;

	std::cout << "   + Comparing dictionaries: dict1 with dict2" << std::endl;
	compareDict1EntriesEqualDict2(fileDict1, fileDict2, "   + ");
	std::cout << "   + OK" << std::endl;

}


void run_sweet_file_dict_write_read()
{
	std::cout << " + run_sweet_file_dict_write_read()" << std::endl;

	int verbosity = 10;
	sweet::Dict save_dict(verbosity);
	create_reference_dict(save_dict);

	std::cout << "   + Writing dictionary" << std::endl;
	save_dict.fileSave("output_test_write.sweet");

	sweet::Dict load_dict(verbosity);
	std::cout << "   + Loading dictionary" << std::endl;
	load_dict.fileLoad("output_test_write.sweet");

	std::cout << "   + Comparing dictionaries: dict with ref" << std::endl;
	compareDict1EntriesEqualDict2(save_dict, load_dict, "   + ");
	std::cout << "   + OK" << std::endl;

	std::cout << "   + Comparing dictionaries: ref with dict" << std::endl;
	compareDict1EntriesEqualDict2(load_dict, save_dict, "   + ");
	std::cout << "   + OK" << std::endl;
}

void hline_star()
{
	std::cout << "********************************************************************************" << std::endl;
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
		hline_star();
		std::cout << "+ Loading from file '" << sweet_file_dict_filename << "'" << std::endl;
		hline_star();

		sweet::Dict fd(sweet_file_dict_filename, false);
		//std::cout << fd << std::endl;

		/*
		 * Should we also run tests on it?
		 */
		if (run_unit_test)
		{
			hline_star();
			std::cout << " + Create reference dictionary and compare with data from file" << std::endl;
			hline_star();
			run_compare_file_dict_with_reference(sweet_file_dict_filename);
		}
	}

	hline_star();
	std::cout << " + Assignment tests" << std::endl;
	hline_star();
	run_sweet_array_assignment_test();

	hline_star();
	std::cout << " + Different dictionary creation tests" << std::endl;
	hline_star();
	run_sweet_file_dict_creation_test();

	hline_star();
	std::cout << " + File write / read tests" << std::endl;
	hline_star();
	run_sweet_file_dict_write_read();

	std::cout << "FIN" << std::endl;

	return EXIT_SUCCESS;
}
