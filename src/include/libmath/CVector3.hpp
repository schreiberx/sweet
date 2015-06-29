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


/*
 * CHANGELOG:
 *
 * 2008-02-27: martin schreiber
 * 	dist2()
 *
 * 2009-06-01: martin schreiber
 * 	min(), max()
 * 
 * 2010-01-18: martin schreiber
 *  updates for doxygen
 */

#ifndef __CVECTOR_HH
	#error "dont include CVector3.hpp directly!"
#endif

#ifndef __CVECTOR3_HH
#define __CVECTOR3_HH

#include <iostream>
#include <cassert>
#include <cmath>


template <typename T> class CMatrix3;

/**
 * \brief	3D Vector handler
 */
template <typename T>
class CVector<3,T>
{
public:
	T data[3];	///< vector data


/******************************************************
 ******************* CONSTRUCTORS *********************
 ******************************************************/

	/**
	 * default constructor
	 */
	inline CVector()
	{
		setZero();
	}

	/**
	 * initialize vector with (x0, x1, x2)
	 */
	inline CVector(const T x0, const T x1, const T x2)
	{
		data[0] = x0;
		data[1] = x1;
		data[2] = x2;
	}

	/**
	 * initialize all vector components with the scalar value 'x'
	 */
	inline CVector(const T x)
	{
		data[0] = x;
		data[1] = x;
		data[2] = x;
	}

	/**
	 * initialize vector components with the array 'v'
	 */
	inline CVector(const CVector<3,T> &v)
	{
		data[0] = v.data[0];
		data[1] = v.data[1];
		data[2] = v.data[2];
	}

	/**
	 * initialize vector components with the array 'v'
	 */
	inline CVector(const CVector<4,T> &v)
	{
		data[0] = v.data[0];
		data[1] = v.data[1];
		data[2] = v.data[2];
	}

	/**
	 * initialize integer vector components with the array 'v' converting to the type (T)
	 */
/*
	inline CVector(const CVector<3,int> &v)
	{
		data[0] = (T)v.data[0];
		data[1] = (T)v.data[1];
		data[2] = (T)v.data[2];
	}
*/


	/**
	 * initialize vector components with the array 'v'
	 */
	inline CVector(const T v[3])
	{
		data[0] = v[0];
		data[1] = v[1];
		data[2] = v[2];
	}




/******************************************************
 ******************* GENERAL FUNCTIONS ****************
 ******************************************************/

	/**
	 * set all components of vector to 0
	 */
	inline void setZero()
	{
		data[0] = T(0);
		data[1] = T(0);
		data[2] = T(0);	
	}


	/**
	 * return normalized (length=1) vector
	 */
	inline CVector<3,T> normal()
	{
		CVector<3,T> v = *this;
		T inv_length = 1.0f/sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		v[0] *= inv_length;
		v[1] *= inv_length;
		v[2] *= inv_length;
		return v;
	}

	/**
	 * compute and return the dot product of this vector and 'v'
	 * \return dot product
	 */
	inline T dotProd(const CVector<3,T> &v)	const
	{
		return v.data[0]*data[0] + v.data[1]*data[1] + v.data[2]*data[2];
	}

	/**
	 * compute and return the cross product of this vector and 'a'
	 * \return cross product
	 */
	inline CVector<3,T> operator%(const CVector<3,T> &a)	const
	{
		return CVector<3,T>(	data[1]*a.data[2] - data[2]*a.data[1],
								data[2]*a.data[0] - data[0]*a.data[2],
								data[0]*a.data[1] - data[1]*a.data[0]
				);
	}

	/**
	 * return square distance to other point
	 */
	inline T dist2(const CVector<3,T> &v)	const
	{
		CVector<3,T> d = CVector<3,T>(v.data[0] - data[0], v.data[1] - data[1], v.data[2] - data[2]);
		return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
	}

	/**
	 * return distance to other point
	 */
	inline T dist(const CVector<3,T> &v)	const
	{
		CVector<3,T> d = CVector<3,T>(v.data[0] - data[0], v.data[1] - data[1], v.data[2] - data[2]);
		return std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
	}

	/**
	 * return elements in cube.
	 * e.g. if the vector components give the resolution of a cube, this function returns the number of cells
	 */
	inline T elements()	const
	{
		return (data[0]*data[1]*data[2]);
	}

	/**
	 * return length of vector
	 */
	inline T getLength()	const
	{
		return std::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
	}

	/**
	 * return (length*length) of vector
	 */
	inline T getLength2()	const
	{
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
	}

	/**
	 * return the minimum vector component
	 */
	inline T min()
	{
		T min = (data[0] < data[1] ? data[0] : data[1]);
		return (data[2] < min ? data[2] : min);
	}

	/**
	 * return the maximum vector component
	 */
	inline T max()
	{
		T max = (data[0] > data[1] ? data[0] : data[1]);
		return (data[2] > max ? data[2] : max);
	}

	/**
	 * normalize the vector
	 */
	inline void normalize()	
	{
		T il = 1/getLength();
		data[0] *= il;
		data[1] *= il;
		data[2] *= il;
	}

	/**
	 * return the normalized vector
	 */
	inline CVector<3,T> getNormalized()	const
	{
		T il = 1/getLength();
		return CVector<3,T>(data[0]*il, data[1]*il, data[2]*il);
	}

	/**
	 * clamp to values -1 and 1
	 */
	inline void clamp1_1()
	{
		data[0] = (data[0] < -1 ? -1 : (data[0] > 1 ? 1 : data[0]));
		data[1] = (data[1] < -1 ? -1 : (data[1] > 1 ? 1 : data[1]));
		data[2] = (data[2] < -1 ? -1 : (data[2] > 1 ? 1 : data[2]));
	}


/******************************************************
 ******************* OPERATORS ************************
 ******************************************************/

	/// assign values of a[2] to this vector and return reference to this vector
	inline CVector<3,T>&	operator=(const T a[3])	{	data[0] = a[0]; data[1] = a[1]; data[2] = a[2];	return *this;	};
	/// assign values of the vector a this vector and return reference to this vector
	inline CVector<3,T>&	operator=(CVector<3,T> const & a)	{	data[0] = a.data[0]; data[1] = a.data[1]; data[2] = a.data[2];	return *this;	};


// T
	/// return new vector (this+a)
	inline CVector<3,T>	operator+(const T a)	const {	return CVector<3,T>(data[0]+a, data[1]+a, data[2]+a);	}
	/// return new vector (this-a)
	inline CVector<3,T>	operator-(const T a)	const {	return CVector<3,T>(data[0]-a, data[1]-a, data[2]-a);	}
	/// return new vector with component wise (this*a)
	inline CVector<3,T>	operator*(const T a)	const {	return CVector<3,T>(data[0]*a, data[1]*a, data[2]*a);	}
	/// return new vector with component wise (this/a)
	inline CVector<3,T>	operator/(const T a)	const {	return CVector<3,T>(data[0]/a, data[1]/a, data[2]/a);	}
	/// add a to this vector and return reference to this vector
	inline CVector<3,T>&	operator+=(const T a)	{	data[0] += a; data[1] += a; data[2] += a;	return *this;	}
	/// subtract a from this vector and return reference to this vector
	inline CVector<3,T>&	operator-=(const T a)	{	data[0] -= a; data[1] -= a; data[2] -= a;	return *this;	}
	/// multiply each component of this vector with scalar a and return reference to this vector
	inline CVector<3,T>&	operator*=(const T a)	{	data[0] *= a; data[1] *= a; data[2] *= a;	return *this;	}
	/// divide each component of this vector by scalar a and return reference to this vector
	inline CVector<3,T>&	operator/=(const T a)	{	data[0] /= a; data[1] /= a; data[2] /= a;	return *this;	}


	/// return new vector with values of this vector multiplied component wise with vector v
	inline CVector<3,T>	operator*(const CVector<3,T> &v)	const
	{
		return CVector<3,T>(data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2]);
	}

	/// return new vector with values of this vector divided component wise by components of vector v
	inline CVector<3,T>	operator/(const CVector<3,T> &v)	const
	{
		return CVector<3,T>(data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2]);
	}

	/// return this vector after adding v
	inline CVector<3,T>&	operator+=(const CVector<3,T> &v)
	{
		data[0] += v.data[0]; data[1] += v.data[1]; data[2] += v.data[2];
		return *this;
	}

	/// return this vector after subtracting v
	inline CVector<3,T>&	operator-=(const CVector<3,T> &v)
	{
		data[0] -= v.data[0];
		data[1] -= v.data[1];
		data[2] -= v.data[2];
		return *this;
	}

	/// return true, if each component of the vector is equal to the corresponding component of vector v
	inline bool	operator==(const CVector<3,T> &v)	{	return bool(data[0] == v.data[0] && data[1] == v.data[1] && data[2] == v.data[2]);	}
	/// return true, if at lease component of the vector is not equal to the corresponding component of vector v
	inline bool	operator!=(const CVector<3,T> &v)	{	return bool(data[0] != v.data[0] || data[1] != v.data[1] || data[2] != v.data[2]);	}

	/**
	 * access element i
	 */
	inline T& operator[](const int i)
	{
		assert(i >= 0 && i < 3);
		return data[i];
	}

	/**
	 * multiply a vector with a matrix.
	 *
	 * this is equivalent to the operation "m^T * this" but without the necessity
	 * to transpose the matrix m.
	 */
	inline CVector<3,T>	operator*(const class CMatrix3<T> &m)	const
	{
		return CVector<3,T>(
				data[0]*m.matrix[0][0] + data[1]*m.matrix[1][0] + data[2]*m.matrix[2][0],
				data[0]*m.matrix[0][1] + data[1]*m.matrix[1][1] + data[2]*m.matrix[2][1],
				data[0]*m.matrix[0][2] + data[1]*m.matrix[1][2] + data[2]*m.matrix[2][2]
			);
	}

/******************************************************
 ******************* SPECIAL STUFF ********************
 ******************************************************/

	/**
	 * \brief	compare set for sort operation
	 */
	struct compareSet
	{
		/**
		 * compare set operator
		 */
		inline bool operator()(CVector<3,T> *v1, CVector<3,T> *v2)
		{
			if ((*v1)[0] != (*v2)[0])
				return (*v1)[0] < (*v2)[0];
			if ((*v1)[1] != (*v2)[1])
				return (*v1)[1] < (*v2)[1];
			if ((*v1)[2] != (*v2)[2])
				return (*v1)[2] < (*v2)[2];
			return false;
		}

		/**
		 * compare set operator
		 */
		inline bool operator()(const CVector<3,T> &v1, const CVector<3,T> &v2)
		{
			if (v1.data[0] != v2.data[0])
				return v1.data[0] < v2.data[0];
			if (v1.data[1] != v2.data[1])
				return v1.data[1] < v2.data[1];
			if (v1.data[2] != v2.data[2])
				return v1.data[2] < v2.data[2];
			return false;
		}
	};
};

/// return new vector subtracting vector v from this vector
template <class T>
inline CVector<3,T>	operator-(const CVector<3,T> &a, const CVector<3,T> &v)
{
	return CVector<3,T>(	a.data[0] - v.data[0],
							a.data[1] - v.data[1],
							a.data[2] - v.data[2]
						);
}

/// return new vector adding v to this vector
template <class T>
inline CVector<3,T>	operator+(const CVector<3,T> &a, const CVector<3,T> &v)
{
	return CVector<3,T>(	a.data[0] + v.data[0],
							a.data[1] + v.data[1],
							a.data[2] + v.data[2]
						);
}


/// return new vector -this
template <class T>
inline CVector<3,T>	operator-(const CVector<3,T> &a)
{
	return CVector<3,T>(-a.data[0], -a.data[1], -a.data[2]);
}

template <class T>
inline
::std::ostream&
operator<<(::std::ostream &co, const CVector<3,T> &v)
{
	return co << "[" << v.data[0] << ", " << v.data[1] << ", " << v.data[2] << "]";
}


#endif
