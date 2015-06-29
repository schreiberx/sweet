/*
 * Copyright 2010 Martin Schreiber
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *      Author: Martin Schreiber (schreiberx@gmail.com)
 */

/**
 * CHANGELOG:
 *
 * 2008-10-28: martin schreiber
 * 	created
 *
 */

#ifndef __CMATRIX_HH
	#error "dont include CMatrix2.hpp directly - use CMatrix.hpp instead!"
#endif

#ifndef CMATRIX2_HPP
#define CMATRIX2_HPP

#include <limits>
#include <iostream>
#include <cmath>
#include "CVector.hpp"


/**
 * \brief	2x2 Matrix Handler
 */
template <typename T>
class CMatrix2
{
public:
	T matrix[2][2];		///< storage for matrix data

/******************************************************
 ******************* CONSTRUCTORS *********************
 ******************************************************/

	/**
	 * constructor: load identity matrix
	 */
	inline CMatrix2()
	{
		loadIdentity();
	}


	/**
	 * constructor: load matrix entries
	 */
	inline CMatrix2(	T v00, T v01,
						T v10, T v11
					)
	{
		matrix[0][0] = v00;
		matrix[0][1] = v01;

		matrix[1][0] = v10;
		matrix[1][1] = v11;
	}


/******************************************************
 ******************* GENERAL FUNCTIONS ****************
 ******************************************************/


	/**
	 * load identity matrix
	 */
	inline void loadIdentity()
	{
		matrix[0][0] = (T)1;
		matrix[0][1] = (T)0;

		matrix[1][0] = (T)0;
		matrix[1][1] = (T)1;
	}



	/**
	 * return inverse of matrix
	 * taken from crystalspace engine
	 * http://www.crystalspace3d.org/docs/online/api/matrix3_8h-source.html
	 *
	 * the inversion makes use of cramers rule
	 */
	CMatrix2<T> getInverse() const
	{
	    return CMatrix2<T>(	matrix[1][1], -matrix[0][1],
							-matrix[1][0], matrix[0][0])
						/	(matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]);
	}

	/**
	 * return transpose of matrix
	 */
	CMatrix2<T> getTranspose() const
	{
		return CMatrix2<T> (
			matrix[0][0], matrix[1][0],
			matrix[0][1], matrix[1][1]
		);
	}

	/**
	 * return the inverse transpose of the matrix
	 */
	CMatrix2<T> getInverseTranspose() const
	{
		return getInverse().getTranspose();
	}


/******************************************************
 ******************* OPERATORS ************************
 ******************************************************/

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix2<T>& operator=(const CMatrix2<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];

		return *this;
	};


	/**
	 * multiply matrix with vector
	 *
	 * \param v	vector
	 */
	inline CMatrix2<T>& operator*=(CVector<2,T> &v)
	{
		matrix[0][0] *= v.data[0];
		matrix[0][1] *= v.data[1];

		matrix[1][0] *= v.data[0];
		matrix[1][1] *= v.data[1];

		return *this;
	};


	/**
	 * multiply matrix with scalar v
	 *
	 * \param v	divisor
	 */
	inline CMatrix2<T> operator*(const T v)
	{
		CMatrix2<T> nm;

		nm[0][0] = matrix[0][0] * v;
		nm[0][1] = matrix[0][1] * v;

		nm[1][0] = matrix[1][0] * v;
		nm[1][1] = matrix[1][1] * v;

		return nm;
	};


	/**
	 * divide each matrix entry by v
	 *
	 * \param v	divisor
	 */
	inline CMatrix2<T> operator/(const T v)
	{
		CMatrix2<T> nm;

		nm[0][0] = matrix[0][0] / v;
		nm[0][1] = matrix[0][1] / v;

		nm[1][0] = matrix[1][0] / v;
		nm[1][1] = matrix[1][1] / v;

		return nm;
	};



	/**
	 * matrix vector product
	 *
	 * \return matrix
	 */
	inline CMatrix2<T>	operator*=(const CVector<2,T> &v)
	{
		matrix[0][0] *= v.data[0];
		matrix[0][1] *= v.data[1];

		matrix[1][0] *= v.data[0];
		matrix[1][1] *= v.data[1];

		return *this;
	}

	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix2<T>	operator*(const CMatrix2<T> &m2)
	{
		CMatrix2<T> m;

		m[0][0] = matrix[0][0]*m2.matrix[0][0] + matrix[0][1]*m2.matrix[1][0];
		m[0][1] = matrix[0][0]*m2.matrix[0][1] + matrix[0][1]*m2.matrix[1][1];

		m[1][0] = matrix[1][0]*m2.matrix[0][0] + matrix[1][1]*m2.matrix[1][0];
		m[1][1] = matrix[1][0]*m2.matrix[0][1] + matrix[1][1]*m2.matrix[1][1];

		return m;
	}


	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix2<T>	operator*=(CMatrix2<T> &m2)
	{
		CMatrix2 m;

		m[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0];
		m[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1];

		m[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0];
		m[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1];

		*this = m;
		return *this;
	}


	/**
	 * linear array access to matrix components
	 */
	inline T* operator[](const int i)
	{
#if DEBUG
		if (i < 0 || i >= 2)
		{
			std::cerr << "OUT OF ARRAY ACCESS!!! creating null exception..." << std::endl;
			*((int*)(0)) = 0;
		}
#endif
		return matrix[i];
	}
};



template <class T>
inline
::std::ostream&
operator<<(::std::ostream &co, const CMatrix2<T> &m)
{
	return co	<< "[" << m.matrix[0][0] << ", " << m.matrix[0][1] << "]" << ::std::endl
				<< "[" << m.matrix[1][0] << ", " << m.matrix[1][1] << "]" << ::std::endl	;
}



#endif
