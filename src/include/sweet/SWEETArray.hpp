/*
 * SWEETArray.hpp
 *
 *  Created on: Feb 13, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SWEETARRAY_HPP_
#define SRC_INCLUDE_SWEET_SWEETARRAY_HPP_

#include <array>
#include <vector>
#include <ostream>


/*
 * SWEET's 1D/2D/3D array data class.
 *
 * This is not intended for HPC, but just to have
 * some 2D container for arbitrary types.
 */

template <int D, typename T>
class SWEETArray
{
	std::array<int,D> _shape;
	std::size_t _size;
	std::vector<T> _data;

public:
	SWEETArray()	:
		_size(0)
	{
		for (int i = 0; i < D; i++)
			_shape[i] = 0;
	}

	SWEETArray(const std::array<int,D> &i_shape)
	{
		setup(i_shape);
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
	T get(int i0, int i1 = -1, int i2 = -1)	const
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

	inline
	T operator()(int i0, int i1=-1, int i2=-1) const {
		return get(i0, i1, i2);
	}


	inline
	T operator[](int i0) const {
		if (D != 1)
			SWEETError("Only 1D supported");

		return _data[i0];
	}

	inline
	SWEETArray<D,T>& operator=(const T *i_values_flat)
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
	SWEETArray<D,T>& operator=(const SWEETArray<D,T> &a)
	{
		this->setup(a._shape);

		for (std::size_t i = 0; i < _size; i++)
			_data[i] = a._data[i];

		return *this;
	}

public:
	friend
	std::ostream& operator<<(std::ostream& os, const SWEETArray<D,T> &a)
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
				os << a.get(i0);

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
					os << a.get(i0, i1);

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
						os << a.get(i0, i1, i2);

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

#endif /* SRC_INCLUDE_SWEET_SWEETARRAY_HPP_ */
