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
 */

#ifndef __CMATRIX_HH
	#error "dont include CMatrix3.hpp directly - use CMatrix.hpp instead!"
#endif

#ifndef CMATRIX3_HPP
#define CMATRIX3_HPP

#include <limits>
#include <iostream>
#include <cmath>
#include "CVector.hpp"

template <typename T>
class CMatrix4;

/**
 * \brief	3x3 Matrix Handler
 *
 * note: this is everything else than an optimized implementation
 *
 * use expression templates if this operations are frequently used in
 * your program!.
 */
template <typename T>
class CMatrix3
{
public:
	T matrix[3][3];	///< storage for data



/******************************************************
 ******************* CONSTRUCTORS *********************
 ******************************************************/

	/**
	 * constructor: load identity matrix
	 */
	inline CMatrix3()
	{
		loadIdentity();
	}


	/**
	 * constructor: load matrix entries
	 */
	inline CMatrix3(	T v00, T v01, T v02,
						T v10, T v11, T v12,
						T v20, T v21, T v22
					)
	{
		matrix[0][0] = v00;
		matrix[0][1] = v01;
		matrix[0][2] = v02;

		matrix[1][0] = v10;
		matrix[1][1] = v11;
		matrix[1][2] = v12;

		matrix[2][0] = v20;
		matrix[2][1] = v21;
		matrix[2][2] = v22;

	}


	/**
	 * constructor: initialize all matrix values to p_value
	 */
	inline CMatrix3(	T p_value	)
	{
		matrix[0][0] = p_value;
		matrix[0][1] = p_value;
		matrix[0][2] = p_value;

		matrix[1][0] = p_value;
		matrix[1][1] = p_value;
		matrix[1][2] = p_value;

		matrix[2][0] = p_value;
		matrix[2][1] = p_value;
		matrix[2][2] = p_value;

	}



/******************************************************
 ******************* GENERAL FUNCTIONS ****************
 ******************************************************/

	/**
	 * load identity matrix
	 */
	inline CMatrix3<T>& loadIdentity()
	{
		matrix[0][0] = (T)1;
		matrix[0][1] = (T)0;
		matrix[0][2] = (T)0;

		matrix[1][0] = (T)0;
		matrix[1][1] = (T)1;
		matrix[1][2] = (T)0;

		matrix[2][0] = (T)0;
		matrix[2][1] = (T)0;
		matrix[2][2] = (T)1;

		return *this;
	}


	/**
	 * create a rotation matrix
	 *
	 * taken from opengl manpage:
	 *
	 *	x^2(1-c)+c	xy(1-c)-zs	xz(1-c)+ys
	 *	yx(1-c)+zs	y^2(1-c)+c	yz(1-c)-xs
	 *	xz(1-c)-ys	yz(1-c)+xs	z^2(1-c)+c
	 */
	void genRotation(	T angle, 		// rotation angle in radians
						const CVector<3,T> axis	// axis of rotation
			)
	{
		T c = std::cos(angle);
		T s = std::sin(angle);
		T cm = T(1)-c;

		CVector<3,T> naxis = axis;
		naxis.normalize();

		matrix[0][0] = naxis[0]*naxis[0]*cm + c;
		matrix[0][1] = naxis[0]*naxis[1]*cm - naxis[2]*s;
		matrix[0][2] = naxis[0]*naxis[2]*cm + naxis[1]*s;

		matrix[1][0] = naxis[1]*naxis[0]*cm + naxis[2]*s;
		matrix[1][1] = naxis[1]*naxis[1]*cm + c;
		matrix[1][2] = naxis[1]*naxis[2]*cm - naxis[0]*s;

		matrix[2][0] = naxis[2]*naxis[0]*cm - naxis[1]*s;
		matrix[2][1] = naxis[2]*naxis[1]*cm + naxis[0]*s;
		matrix[2][2] = naxis[2]*naxis[2]*cm + c;
	}


	/**
	 * return inverse of matrix
	 * taken from crystalspace engine
	 * http://www.crystalspace3d.org/docs/online/api/matrix3_8h-source.html
	 *
	 * the inversion makes use of cramers rule
	 */
	inline CMatrix3<T> getInverse() const
	{
	    CMatrix3<T> nm(
	             (matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]), -(matrix[0][1]*matrix[2][2] - matrix[0][2]*matrix[2][1]),  (matrix[0][1]*matrix[1][2] - matrix[0][2]*matrix[1][1]),
	            -(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]),  (matrix[0][0]*matrix[2][2] - matrix[0][2]*matrix[2][0]), -(matrix[0][0]*matrix[1][2] - matrix[0][2]*matrix[1][0]),
	             (matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]), -(matrix[0][0]*matrix[2][1] - matrix[0][1]*matrix[2][0]),  (matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0])
	            );
	    T s = (T)1/(matrix[0][0]*nm.matrix[0][0] + matrix[0][1]*nm.matrix[1][0] + matrix[0][2]*nm.matrix[2][0]);

	    return nm*s;
	}

	/**
	 * return transpose of matrix
	 */
	CMatrix3<T> getTranspose() const
	{
		return CMatrix3<T> (
			matrix[0][0], matrix[1][0], matrix[2][0],
			matrix[0][1], matrix[1][1], matrix[2][1],
			matrix[0][2], matrix[1][2], matrix[2][2]
		);
	}

	/**
	 * return the inverse transpose of the matrix
	 */
	CMatrix3<T> getInverseTranspose() const
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
	inline CMatrix3<T>& operator=(CMatrix3<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];

		return *this;
	};


	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3<T>& operator=(CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];

		return *this;
	};



	/**
	 * multiply matrix with vector
	 *
	 * \param v	vector
	 */
	inline CMatrix3<T>& operator*=(const CVector<3,T> &v)
	{
		matrix[0][0] *= v[0];
		matrix[0][1] *= v[1];
		matrix[0][2] *= v[2];

		matrix[1][0] *= v[0];
		matrix[1][1] *= v[1];
		matrix[1][2] *= v[2];

		matrix[2][0] *= v[0];
		matrix[2][1] *= v[1];
		matrix[2][2] *= v[2];

		return *this;
	};



	/**
	 * add a matrix component wise
	 *
	 * \param v	vector
	 */
	inline CMatrix3<T>& operator+=(const CMatrix3<T> &m)
	{
		matrix[0][0] += m.matrix[0][0];
		matrix[0][1] += m.matrix[0][1];
		matrix[0][2] += m.matrix[0][2];

		matrix[1][0] += m.matrix[1][0];
		matrix[1][1] += m.matrix[1][1];
		matrix[1][2] += m.matrix[1][2];

		matrix[2][0] += m.matrix[2][0];
		matrix[2][1] += m.matrix[2][1];
		matrix[2][2] += m.matrix[2][2];

		return *this;
	};



	/**
	 * multiply matrix with vector
	 *
	 * \param v	vector
	 */
	inline CMatrix3<T> operator*(const T v)
	{
		CMatrix3<T> nm;

		nm[0][0] = matrix[0][0] * v;
		nm[0][1] = matrix[0][1] * v;
		nm[0][2] = matrix[0][2] * v;

		nm[1][0] = matrix[1][0] * v;
		nm[1][1] = matrix[1][1] * v;
		nm[1][2] = matrix[1][2] * v;

		nm[2][0] = matrix[2][0] * v;
		nm[2][1] = matrix[2][1] * v;
		nm[2][2] = matrix[2][2] * v;

		return nm;
	};

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3<T>& operator=(const CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];
		return *this;
	};

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3(const CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];
	};

	/**
	 * matrix vector product
	 *
	 * use 1 as 3rd component of input vector v
	 *
	 * avoid using parameters by references because this cannot be handled by
	 * subsequenced orders like a+(o*i)
	 */
	inline CVector<2,T>	operator*(	CVector<2,T> v	)		const
	{
		return CVector<2,T>(	matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2],
								matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]
				);
	}

	/**
	 * matrix vector product
	 */
	inline CVector<3,T>	operator*(	const CVector<3,T> &v	)	const
	{
		return CVector<3,T>(
					matrix[0][0]*v.data[0] + matrix[0][1]*v.data[1] + matrix[0][2]*v.data[2],
					matrix[1][0]*v.data[0] + matrix[1][1]*v.data[1] + matrix[1][2]*v.data[2],
					matrix[2][0]*v.data[0] + matrix[2][1]*v.data[1] + matrix[2][2]*v.data[2]
				);
	}


	/**
	 * matrix vector product
	 *
	 * \return matrix
	 */
	inline CMatrix3<T>	operator*=(CVector<3,T> &v)
	{
		matrix[0][0] *= v[0];
		matrix[0][1] *= v[1];
		matrix[0][2] *= v[2];

		matrix[1][0] *= v[0];
		matrix[1][1] *= v[1];
		matrix[1][2] *= v[2];

		matrix[2][0] *= v[0];
		matrix[2][1] *= v[1];
		matrix[2][2] *= v[2];

		return *this;
	}

	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix3<T>	operator*(const CMatrix3<T> &m2)	const
	{
		CMatrix3 m;

		m[0][0] = matrix[0][0]*m2.matrix[0][0] + matrix[0][1]*m2.matrix[1][0] + matrix[0][2]*m2.matrix[2][0];
		m[0][1] = matrix[0][0]*m2.matrix[0][1] + matrix[0][1]*m2.matrix[1][1] + matrix[0][2]*m2.matrix[2][1];
		m[0][2] = matrix[0][0]*m2.matrix[0][2] + matrix[0][1]*m2.matrix[1][2] + matrix[0][2]*m2.matrix[2][2];

		m[1][0] = matrix[1][0]*m2.matrix[0][0] + matrix[1][1]*m2.matrix[1][0] + matrix[1][2]*m2.matrix[2][0];
		m[1][1] = matrix[1][0]*m2.matrix[0][1] + matrix[1][1]*m2.matrix[1][1] + matrix[1][2]*m2.matrix[2][1];
		m[1][2] = matrix[1][0]*m2.matrix[0][2] + matrix[1][1]*m2.matrix[1][2] + matrix[1][2]*m2.matrix[2][2];

		m[2][0] = matrix[2][0]*m2.matrix[0][0] + matrix[2][1]*m2.matrix[1][0] + matrix[2][2]*m2.matrix[2][0];
		m[2][1] = matrix[2][0]*m2.matrix[0][1] + matrix[2][1]*m2.matrix[1][1] + matrix[2][2]*m2.matrix[2][1];
		m[2][2] = matrix[2][0]*m2.matrix[0][2] + matrix[2][1]*m2.matrix[1][2] + matrix[2][2]*m2.matrix[2][2];

		return m;
	}


	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix3<T>	operator*=(CMatrix3<T> &m2)
	{
		CMatrix3 m;

		m[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0] + matrix[0][2]*m2[2][0];
		m[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1] + matrix[0][2]*m2[2][1];
		m[0][2] = matrix[0][0]*m2[0][2] + matrix[0][1]*m2[1][2] + matrix[0][2]*m2[2][2];

		m[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0] + matrix[1][2]*m2[2][0];
		m[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1] + matrix[1][2]*m2[2][1];
		m[1][2] = matrix[1][0]*m2[0][2] + matrix[1][1]*m2[1][2] + matrix[1][2]*m2[2][2];

		m[2][0] = matrix[2][0]*m2[0][0] + matrix[2][1]*m2[1][0] + matrix[2][2]*m2[2][0];
		m[2][1] = matrix[2][0]*m2[0][1] + matrix[2][1]*m2[1][1] + matrix[2][2]*m2[2][1];
		m[2][2] = matrix[2][0]*m2[0][2] + matrix[2][1]*m2[1][2] + matrix[2][2]*m2[2][2];

		*this = m;
		return *this;
	}



	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix3<T>	operator*=(CMatrix4<T> &m2)
	{
		CMatrix3 m;

		m[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0] + matrix[0][2]*m2[2][0];
		m[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1] + matrix[0][2]*m2[2][1];
		m[0][2] = matrix[0][0]*m2[0][2] + matrix[0][1]*m2[1][2] + matrix[0][2]*m2[2][2];

		m[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0] + matrix[1][2]*m2[2][0];
		m[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1] + matrix[1][2]*m2[2][1];
		m[1][2] = matrix[1][0]*m2[0][2] + matrix[1][1]*m2[1][2] + matrix[1][2]*m2[2][2];

		m[2][0] = matrix[2][0]*m2[0][0] + matrix[2][1]*m2[1][0] + matrix[2][2]*m2[2][0];
		m[2][1] = matrix[2][0]*m2[0][1] + matrix[2][1]*m2[1][1] + matrix[2][2]*m2[2][1];
		m[2][2] = matrix[2][0]*m2[0][2] + matrix[2][1]*m2[1][2] + matrix[2][2]*m2[2][2];

		*this = m;
		return *this;
	}

	/**
	 * set every entry to zero
	 */
	inline CMatrix3<T>& setZero()
	{
		matrix[0][0] = 0;
		matrix[0][1] = 0;
		matrix[0][2] = 0;

		matrix[1][0] = 0;
		matrix[1][1] = 0;
		matrix[1][2] = 0;

		matrix[2][0] = 0;
		matrix[2][1] = 0;
		matrix[2][2] = 0;

		return *this;
	}

	/**
	 * setup matrix to behave like the left side of a cross product
	 */
	inline CMatrix3<float>& setupCrossProduct(CVector<3,float> &v)
	{
		matrix[0][0] = 0;
		matrix[0][1] = -v[2];
		matrix[0][2] = v[1];

		matrix[1][0] = v[2];
		matrix[1][1] = 0;
		matrix[1][2] = -v[0];

		matrix[2][0] = -v[1];
		matrix[2][1] = v[0];
		matrix[2][2] = 0;

		return *this;
	}

	/**
	 * linear array access to matrix components
	 */
	inline T* operator[](const int i)
	{
#if DEBUG
		if (i < 0 || i >= 3)
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
operator<<(::std::ostream &co, const CMatrix3<T> &m)
{
	return co	<< "[" << m.matrix[0][0] << ", " << m.matrix[0][1] << ", " << m.matrix[0][2] << "]" << ::std::endl
				<< "[" << m.matrix[1][0] << ", " << m.matrix[1][1] << ", " << m.matrix[1][2] << "]" << ::std::endl
				<< "[" << m.matrix[2][0] << ", " << m.matrix[2][1] << ", " << m.matrix[2][2] << "]" << ::std::endl	;
}



#endif
