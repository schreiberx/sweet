/*
 *
 * ArrayND.hpp
 *  Created on: Feb 13, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_ARRAYND_HPP_
#define SRC_INCLUDE_SWEET_ARRAYND_HPP_

#include <array>
#include <vector>
#include <ostream>


/*
 * SWEET's 1D/2D/3D array data class for dictionaries.
 *
 * This is not intended for HPC, but just to have
 * some D dimensional container for arbitrary types.
 */

namespace sweet
{

template <int D, typename T>
class ArrayND
{
	std::array<int,D> _shape;
	std::size_t _size;
	std::vector<T> _data;

public:
	ArrayND()	:
		_size(0)
	{
		for (int i = 0; i < D; i++)
			_shape[i] = 0;
	}

	ArrayND(const std::array<int,D> &i_shape)
	{
		setup(i_shape);
	}


	ArrayND(const std::array<int,D> &i_shape, const T i_data[])
	{
		setup(i_shape);
		operator=(i_data);
	}

	void setup(const std::array<int,D> &i_shape)
	{
		if (D < 1 || D > 3)
			SWEETError("Only 1D, 2D or 3D are supported!");

		_shape = i_shape;

		_size = 1;
		for (int i = 0; i < D; i++)
			_size *= _shape[i];

		_data.resize(_size);
	}

	void setup(const std::array<int,D> &i_shape, T i_data[])
	{
		setup(i_shape);
		operator=(i_data);
	}

	std::size_t size()
	{
		return _size;
	}

	T *data()
	{
		return _data.data();
	}

	const std::array<int,D> &shape()
	{
		return _shape;
	}


	inline
	void set(int i0, const T &i_value)
	{
		if (D != 1)
			SWEETError("Only for 1D");

		_data[i0] = i_value;
	}

	inline
	void set(int i0, int i1, const T &i_value)
	{
		if (D != 2)
			SWEETError("Only for 2D");

		_data[i0*_shape[1] + i1] = i_value;
	}

	inline
	void set3(int i0, int i1, int i2, const T &i_value)
	{
		if (D != 3)
			SWEETError("Only for 3D");

		_data[i0*_shape[1]*_shape[2] + i1*_shape[2] + i2] = i_value;
	}

	/*
	 * Special getter which is just constant
	 *
	 * This is required if called by other constant functions from this class.
	 */
	const T& getConst(int i0, int i1=-1, int i2=-1) const
	{
		if (D == 1)
		{
			return _data[i0];
		}
		else if (D == 2)
		{
			return _data[i0*_shape[1] + i1];
		}
		else if (D == 3)
		{
			return _data[i0*_shape[1]*_shape[2] + i1*_shape[2] + i2];
		}
		else
		{
			SWEETError("Not supported!");
		}
	}

	T& get(int i0, int i1=-1, int i2=-1)
	{
		if (D == 1)
		{
			return _data[i0];
		}
		else if (D == 2)
		{
			return _data[i0*_shape[1] + i1];
		}
		else if (D == 3)
		{
			return _data[i0*_shape[1]*_shape[2] + i1*_shape[2] + i2];
		}
		else
		{
			SWEETError("Not supported!");
		}
	}


	/*
	 * This is just a convenience handler to use allow a very compact
	 * access to this array rather than using .get(...)
	 *
	 * The array operator[] is no option for us since this only supports one
	 * argument in C++11
	 */
	inline
	const T& operator()(int i0, int i1=-1, int i2=-1)	const
	{
		return getConst(i0, i1, i2);
	}

#if 0
	/*
	 * 1D Array access
	 */
	inline
	T operator[](int i0) const {
		if (D != 1)
			SWEETError("Only 1D supported");

		return _data[i0];
	}
#endif

	inline
	ArrayND<D,T>& operator=(const T *i_values_flat)
	{
		for (int i = 0; i < D; i++)
			if (_shape[i] == 0)
				SWEETError("Shape is 0, you need to resize array before assigning raw data!");

		// We simply hope that the data is properly allocated
		for (std::size_t i = 0; i < _size; i++)
			_data[i] = i_values_flat[i];

		return *this;
	}


	inline
	ArrayND<D,T>& operator=(const ArrayND<D,T> &a)
	{
		this->setup(a._shape);

		for (std::size_t i = 0; i < _size; i++)
			_data[i] = a._data[i];

		return *this;
	}

public:
	friend
	std::ostream& operator<<(std::ostream& os, const ArrayND<D,T> &a)
	{
		if (a._size == 0)
		{
			os << "(empty)";
			return os;
		}

		if (D == 1)
		{
			std::cout << "[";
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				os << a.getConst(i0);

				if (i0 != a._shape[0]-1)
					os << ",\t";
			}
			std::cout << "]";
			os << std::endl;
		}
		else if (D == 2)
		{
			std::cout << "[" << std::endl;
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				std::cout << "\t[";
				for (int i1 = 0; i1 < a._shape[1]; i1++)
				{
					os << a.getConst(i0, i1);

					if (i1 != a._shape[1]-1)
						os << ",\t";
				}
				std::cout << "]";
				os << std::endl;
			}
			std::cout << "]";
		}
		else if (D == 3)
		{
			std::cout << "[" << std::endl;
			for (int i0 = 0; i0 < a._shape[0]; i0++)
			{
				std::cout << "\t[" << std::endl;
				for (int i1 = 0; i1 < a._shape[1]; i1++)
				{
					std::cout << "\t\t[";
					for (int i2 = 0; i2 < a._shape[2]; i2++)
					{
						os << a.getConst(i0, i1, i2);

						if (i2 != a._shape[2]-1)
							os << ",\t";
					}
					std::cout << "]";
					os << std::endl;
				}
				std::cout << "\t]";
				os << std::endl;
			}
			std::cout << "]";
		}
		else
		{
			SWEETError("Not supported!");
		}
		return os;
	}
};

}

#endif /* SRC_INCLUDE_SWEET_ArrayND_HPP_ */
