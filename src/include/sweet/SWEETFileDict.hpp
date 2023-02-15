/*
 * SWEETFileDict.hpp
 *
 *  Created on: Feb 13, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_FILE_DICT_HPP__
#define SRC_INCLUDE_SWEET_FILE_DICT_HPP__

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/SWEETError.hpp>
#include <sweet/SWEETArray.hpp>

class SWEETFileDict
{
public:
	typedef long long int64;
	typedef double float64;
	typedef std::complex<double> complex128;

private:
	enum SWEETFileDict_ElementTypes
	{
		SWEET_FILE_DICT_STRING = 100,
		SWEET_FILE_DICT_INT64 = 230,
		SWEET_FILE_DICT_FLOAT64 = 240,
		SWEET_FILE_DICT_ARRAY_1D_FLOAT64 = 401,
		SWEET_FILE_DICT_ARRAY_2D_FLOAT64 = 402,
		SWEET_FILE_DICT_ARRAY_3D_FLOAT64 = 403,
		SWEET_FILE_DICT_ARRAY_1D_COMPLEX128 = 501,
		SWEET_FILE_DICT_ARRAY_2D_COMPLEX128 = 502,
		SWEET_FILE_DICT_ARRAY_3D_COMPLEX128 = 503,
	};


	class SWEETFileDict_Element
	{
	public:
		std::string key;

		SWEETFileDict_ElementTypes type_id;

// Union doesn't work with C++ classes :-(
//		union
//		{
			std::string value_str;
			int64 value_scalar_int64;
			double value_scalar_float64;

			SWEETArray<1,double> value_array_1d_float64;
			SWEETArray<2,double> value_array_2d_float64;
			SWEETArray<3,double> value_array_3d_float64;

			SWEETArray<1,complex128> value_array_1d_complex128;
			SWEETArray<2,complex128> value_array_2d_complex128;
			SWEETArray<3,complex128> value_array_3d_complex128;
//		};

		template <typename T>
		void getValue(T &o_value)	const;

		void getValue(std::string &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_STRING)
				SWEETError("Type mismatch!");

			o_value = value_str;
		}

		void getValue(int64 &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = value_scalar_int64;
		}

		void getValue(float64 &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = value_scalar_float64;
		}

		void getValue(SWEETArray<1,float64> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_1D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = value_array_1d_float64;
		}

		void getValue(SWEETArray<2,float64> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_2D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = value_array_2d_float64;
		}

		void getValue(SWEETArray<3,float64> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_3D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = value_array_3d_float64;
		}

		void getValue(SWEETArray<1,complex128> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_1D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = value_array_1d_complex128;
		}

		void getValue(SWEETArray<2,complex128> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_2D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = value_array_2d_complex128;
		}

		void getValue(SWEETArray<3,complex128> &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_ARRAY_3D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = value_array_3d_complex128;
		}

		/*
		 * Special handlers which automatically convert values
		 */
		void getValue(bool &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = (bool)value_scalar_int64;
		}
		void getValue(int &o_value)	const
		{
			if (type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = (int)value_scalar_int64;
		}
	};


	std::vector<SWEETFileDict_Element> _dict;

	bool _debug;

public:
	SWEETFileDict()	:
		_debug(false)
	{
	}

	SWEETFileDict(const std::string &i_filename, bool i_debug = false)
	{
		_debug = i_debug;

		loadFromFile(i_filename);
	}


	/*
	 * Load from a SWEETFileDict
	 *
	 * See corresponding Python file for description of file format
	 */
	bool loadFromFile(const std::string &i_filename)
	{
		if (sizeof(int64) != 8)
			SWEETError("Something is weird!");

		std::ifstream f(i_filename, std::ios::in | std::ios::binary);

		if (!f.is_open())
			SWEETError(std::string("Unable to open file ")+i_filename);

		std::string magic_start = _read_str0(f);

		if (magic_start != "SWEETFileDict")
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_start << "' at beginning!";
			SWEETError(ss.str());
		}

		std::size_t num_entries = _read<int64>(f);

		if (_debug)
			std::cout << "SWEETFileDict: Found " << num_entries << " dictionary entries" << std::endl;

		_dict.resize(num_entries);
		for (std::size_t i = 0; i < num_entries; i++)
		{
			SWEETFileDict_Element &e = _dict[i];

			e.key = _read_str0(f);

			e.type_id = (SWEETFileDict_ElementTypes)_read<int64>(f);

			switch(e.type_id)
			{
				case SWEET_FILE_DICT_STRING:
					e.value_str = _read_str0(f);
					break;

				case SWEET_FILE_DICT_INT64:
					e.value_scalar_int64 = _read<int64>(f);
					break;

				case SWEET_FILE_DICT_FLOAT64:
					e.value_scalar_float64 = _read<double>(f);
					break;

				case SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
				{
					std::array<int,1> shape;
					shape[0] = _read<int64>(f);

					e.value_array_1d_float64.setup(shape);
					_read_array(
							f,
							e.value_array_1d_float64
						);
				}
					break;

				case SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
				{
					std::array<int,2> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);

					e.value_array_2d_float64.setup(shape);
					_read_array(
							f,
							e.value_array_2d_float64
						);
				}
					break;

				case SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
				{
					std::array<int,3> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					shape[2] = _read<int64>(f);

					e.value_array_3d_float64.setup(shape);
					_read_array(
							f,
							e.value_array_3d_float64
						);
				}
					break;

				case SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
				{
					std::array<int,1> shape;
					shape[0] = _read<int64>(f);

					e.value_array_1d_complex128.setup(shape);
					_read_array(
							f,
							e.value_array_1d_complex128
						);
				}
					break;

				case SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
				{
					std::array<int,2> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);

					e.value_array_2d_complex128.setup(shape);
					_read_array(
							f,
							e.value_array_2d_complex128
						);
				}
					break;

				case SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
				{
					std::array<int,3> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					shape[2] = _read<int64>(f);

					e.value_array_3d_complex128.setup(shape);
					_read_array(
							f,
							e.value_array_3d_complex128
						);
				}
					break;

				default:
				{
					std::ostringstream ss;
					ss << "Unknown type id '" << e.type_id << "'";
					SWEETError(ss.str());
				}
			}

			if (_debug)
				std::cout << " + Processing key '" << e.key << "'" << std::endl;
		}

		std::string magic_end = _read_str0(f);

		if (magic_start != "SWEETFileDict")
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_end << "' at end!";
			SWEETError(ss.str());
		}

		return false;
	}

	void print(std::ostream& os = std::cout)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			const SWEETFileDict_Element &e = _dict[i];

			os << " + key: '" << e.key << "'" << std::endl;
			os << "   value: ";

			switch(e.type_id)
			{
				default:
					SWEETError("Unknown type");
					break;

				case SWEET_FILE_DICT_STRING:
					os << "'" << e.value_str << "'" << std::endl;
					break;

				case SWEET_FILE_DICT_INT64:
					os << "'" << e.value_scalar_int64 << "'" << std::endl;
					break;

				case SWEET_FILE_DICT_FLOAT64:
					os << "'" << e.value_scalar_float64 << "'" << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
					os << std::endl;
					os << e.value_array_1d_float64 << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
					os << std::endl;
					os << e.value_array_2d_float64 << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
					os << std::endl;
					os << e.value_array_3d_float64 << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
					os << std::endl;
					os << e.value_array_1d_complex128 << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
					os << std::endl;
					os << e.value_array_2d_complex128 << std::endl;
					break;

				case SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
					os << std::endl;
					os << e.value_array_3d_complex128 << std::endl;
					break;
			}
		}

	}


	template <typename T>
	void getValue(const std::string &i_key, T &o_value)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].key == i_key)
			{
				_dict[i].getValue(o_value);
				return;
			}
		}

		std::ostringstream ss;
		ss << "Key '" << i_key << "' not found!";
		SWEETError(ss.str());
	}


public:
	friend
	std::ostream& operator<<(std::ostream& os, const SWEETFileDict &fd)
	{
		fd.print(os);
		return os;
	}


private:
	/**
	 * Read 0 terminated string from file
	 */
	std::string _read_str0(std::ifstream &f)
	{
		std::vector<char> buffer;
		buffer.resize(1024);

		bool found = false;
		for (std::size_t i = 0; i < buffer.size()-1; i++)
		{
			f.read((char*)&(buffer[i]), sizeof(char));
			if (buffer[i] == 0)
			{
				found = true;
				break;
			}
		}

		if (!found)
			SWEETError("0 termination character not found, stopping here");

		// convert to std::string
		return (std::string)(char*)buffer.data();
	}

	template <typename T>
	T _read(std::ifstream &f)
	{
		T retval;
		f.read((char*)&retval, sizeof(retval));
		return retval;
	}


	template <int D, typename T>
	void _read_array(
			std::ifstream &f,
			SWEETArray<D,T> &array
	)
	{
		for (std::size_t i = 0; i < array.size(); i++)
			array.data()[i] = _read<T>(f);
	}
};



#endif
