/*
 * Dict.hpp
 *
 *  Created on: Feb 13, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_FILE_DICTIONARY_HPP__
#define SRC_INCLUDE_SWEET_FILE_DICTIONARY_HPP__

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/SWEETError.hpp>
#include <sweet/DictArrayND.hpp>

namespace sweet
{

class Dict
{
public:
	typedef long long int64;
	typedef double float64;
	typedef std::complex<double> complex128;


	class Dict_Element
	{
	public:
		enum ElementTypes
		{
			SWEET_FILE_DICT_NONE = 0,

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
	private:

		std::string _key;

		ElementTypes _type_id;

		std::string _value_str;
		int64 _value_scalar_int64;
		double _value_scalar_float64;

		sweet::ArrayND<1,double> _value_array_1d_float64;
		sweet::ArrayND<2,double> _value_array_2d_float64;
		sweet::ArrayND<3,double> _value_array_3d_float64;

		sweet::ArrayND<1,complex128> _value_array_1d_complex128;
		sweet::ArrayND<2,complex128> _value_array_2d_complex128;
		sweet::ArrayND<3,complex128> _value_array_3d_complex128;

public:
		template <typename T>
		Dict_Element(
				const std::string &i_key,
				const T &i_value
		)
		{
			set(i_key, i_value);
		}


		Dict_Element()	:
			_type_id(ElementTypes::SWEET_FILE_DICT_NONE)
		{
		}

		template <typename T>
		void get(T &o_value)	const;

		template <typename T>
		void set(const std::string &i_key, const T &i_value);

		void get(std::string &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_STRING)
				SWEETError("Type mismatch!");

			o_value = _value_str;
		}

		void set(const std::string &i_key, const std::string &i_value)
		{
			_type_id = SWEET_FILE_DICT_STRING;
			_key = i_key;
			_value_str = i_value;
		}

		void get(int64 &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = _value_scalar_int64;
		}

		void set(const std::string &i_key, const int64 &i_value)
		{
			_type_id = SWEET_FILE_DICT_INT64;
			_key = i_key;
			_value_scalar_int64 = i_value;
		}

		void get(float64 &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = _value_scalar_float64;
		}

		void set(const std::string &i_key, const float64 &i_value)
		{
			_type_id = SWEET_FILE_DICT_FLOAT64;
			_key = i_key;
			_value_scalar_float64 = i_value;
		}

		void get(sweet::ArrayND<1,float64> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_1D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = _value_array_1d_float64;
		}

		void set(const std::string &i_key, const sweet::ArrayND<1,float64> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_1D_FLOAT64;
			_key = i_key;
			_value_array_1d_float64 = i_value;
		}

		void get(sweet::ArrayND<2,float64> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_2D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = _value_array_2d_float64;
		}

		void set(const std::string &i_key, const sweet::ArrayND<2,float64> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_2D_FLOAT64;
			_key = i_key;
			_value_array_2d_float64 = i_value;
		}

		void get(sweet::ArrayND<3,float64> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_3D_FLOAT64)
				SWEETError("Type mismatch!");

			o_value = _value_array_3d_float64;
		}

		void set(const std::string &i_key, const sweet::ArrayND<3,float64> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_3D_FLOAT64;
			_key = i_key;
			_value_array_3d_float64 = i_value;
		}

		void get(sweet::ArrayND<1,complex128> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_1D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = _value_array_1d_complex128;
		}

		void set(const std::string &i_key, const sweet::ArrayND<1,complex128> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_1D_COMPLEX128;
			_key = i_key;
			_value_array_1d_complex128 = i_value;
		}

		void get(sweet::ArrayND<2,complex128> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_2D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = _value_array_2d_complex128;
		}

		void set(const std::string &i_key, const sweet::ArrayND<2,complex128> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_2D_COMPLEX128;
			_key = i_key;
			_value_array_2d_complex128 = i_value;
		}

		void get(sweet::ArrayND<3,complex128> &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_ARRAY_3D_COMPLEX128)
				SWEETError("Type mismatch!");

			o_value = _value_array_3d_complex128;
		}

		void set(const std::string &i_key, const sweet::ArrayND<3,complex128> &i_value)
		{
			_type_id = SWEET_FILE_DICT_ARRAY_3D_COMPLEX128;
			_key = i_key;
			_value_array_3d_complex128 = i_value;
		}

		/*
		 * Special handlers which automatically convert values
		 */

		void set(const std::string &i_key, const char *i_value)
		{
			_type_id = SWEET_FILE_DICT_STRING;
			_key = i_key;
			_value_str = i_value;
		}

		void get(bool &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = (bool)_value_scalar_int64;
		}

#if 0
		void set(const std::string &i_key, const bool &i_value)
		{
			_type_id = SWEET_FILE_DICT_INT64;
			_key = i_key;
			_value_scalar_int64 = (int64)i_value;
		}
#endif

		void get(int &o_value)	const
		{
			if (_type_id != SWEET_FILE_DICT_INT64)
				SWEETError("Type mismatch!");

			o_value = (int)_value_scalar_int64;
		}

		void set(const std::string &i_key, const int &i_value)
		{
			_type_id = SWEET_FILE_DICT_INT64;
			_key = i_key;
			_value_scalar_int64 = (int64)i_value;
		}


		const std::string &getKey() const
		{
			return _key;
		}


		friend
		std::ostream& operator<<(std::ostream& os, const Dict_Element &e)
		{
			os << "'" << e._key << "' => ";

			switch(e._type_id)
			{
				default:
					SWEETError("Unknown type");
					break;

				case SWEET_FILE_DICT_STRING:
					os << "'" << e._value_str << "'";
					break;

				case SWEET_FILE_DICT_INT64:
					os << "'" << e._value_scalar_int64 << "'";
					break;

				case SWEET_FILE_DICT_FLOAT64:
					os << "'" << e._value_scalar_float64 << "'";
					break;

				case SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
					os << std::endl;
					os << e._value_array_1d_float64;
					break;

				case SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
					os << std::endl;
					os << e._value_array_2d_float64;
					break;

				case SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
					os << std::endl;
					os << e._value_array_3d_float64;
					break;

				case SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
					os << std::endl;
					os << e._value_array_1d_complex128;
					break;

				case SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
					os << std::endl;
					os << e._value_array_2d_complex128;
					break;

				case SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
					os << std::endl;
					os << e._value_array_3d_complex128;
					break;
			}

			return os;
		}
	};

	std::vector<Dict_Element> _dict;

	bool _debug;

public:
	Dict()	:
		_debug(false)
	{
	}

	Dict(const std::string &i_filename, bool i_debug = false)
	{
		_debug = i_debug;

		loadFromFile(i_filename);
	}


	/*
	 * Load from a Dict
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

		if (magic_start != "Dict")
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_start << "' at beginning!";
			SWEETError(ss.str());
		}

		std::size_t num_entries = _read<int64>(f);

		if (_debug)
			std::cout << "Dict: Found " << num_entries << " dictionary entries" << std::endl;

		_dict.resize(num_entries);
		for (std::size_t i = 0; i < num_entries; i++)
		{
			Dict_Element &e = _dict[i];

			std::string key = _read_str0(f);

			Dict_Element::ElementTypes type_id = (Dict_Element::ElementTypes)_read<int64>(f);

			switch(type_id)
			{
				case Dict_Element::SWEET_FILE_DICT_STRING:
					e.set(key, _read_str0(f));
					break;

				case Dict_Element::SWEET_FILE_DICT_INT64:
					e.set(key, _read<int64>(f));
					break;

				case Dict_Element::SWEET_FILE_DICT_FLOAT64:
					e.set(key, _read<double>(f));
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_1D_FLOAT64:
				{
					std::array<int,1> shape;
					shape[0] = _read<int64>(f);
					sweet::ArrayND<1, float64> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_2D_FLOAT64:
				{
					std::array<int,2> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					sweet::ArrayND<2, float64> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_3D_FLOAT64:
				{
					std::array<int,3> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					shape[2] = _read<int64>(f);
					sweet::ArrayND<3, float64> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_1D_COMPLEX128:
				{
					std::array<int,1> shape;
					shape[0] = _read<int64>(f);
					sweet::ArrayND<1, complex128> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_2D_COMPLEX128:
				{
					std::array<int,2> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					sweet::ArrayND<2, complex128> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				case Dict_Element::SWEET_FILE_DICT_ARRAY_3D_COMPLEX128:
				{
					std::array<int,3> shape;
					shape[0] = _read<int64>(f);
					shape[1] = _read<int64>(f);
					shape[2] = _read<int64>(f);
					sweet::ArrayND<3, complex128> array(shape);

					_read_array(f, array);
					e.set(key, array);
				}
					break;

				default:
				{
					std::ostringstream ss;
					ss << "Unknown type id '" << type_id << "'";
					SWEETError(ss.str());
				}
			}

			if (_debug)
				std::cout << " + Processing key '" << key << "'" << std::endl;
		}

		std::string magic_end = _read_str0(f);

		if (magic_start != "Dict")
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
			const Dict_Element &e = _dict[i];

			std::cout << " + " << e << std::endl;
		}

	}


	template <typename T>
	void get(const std::string &i_key, T &o_value)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
			{
				_dict[i].get(o_value);
				return;
			}
		}

		std::ostringstream ss;
		ss << "Key '" << i_key << "' not found!";
		SWEETError(ss.str());
	}



	bool keyExists(const std::string &i_key)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return true;
		}

		return false;
	}

	int keyIndex(const std::string &i_key)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return i;
		}

		return -1;
	}


	template <typename T>
	void set(const std::string &i_key, const T &i_value, bool i_ignore_existing = true)
	{
		int idx = keyIndex(i_key);

		if (idx < 0)
		{
			_dict.push_back(Dict_Element(i_key, i_value));
			return;
		}

		if (!i_ignore_existing)
			SWEETError("Key already exists");

		_dict[idx].set(i_key, i_value);
	}


public:
	friend
	std::ostream& operator<<(std::ostream& os, const Dict &fd)
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
			sweet::ArrayND<D,T> &array
	)
	{
		for (std::size_t i = 0; i < array.size(); i++)
			array.data()[i] = _read<T>(f);
	}
};

}

#endif
