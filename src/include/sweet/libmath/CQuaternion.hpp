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
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> (schreiberx@gmail.com)
 */

#ifndef __CQUATERNION_HPP__
#define __CQUATERNION_HPP__

#include <array>

/**
 * implementation of a quaternion for rotations without gimbal lock
 */
template<typename T>
class CQuaternion
{
public:
	T i, j, k, w;

	CQuaternion()
	{
		reset();
	}


	/**
	 * constructor
	 *
	 * \param axix	the axis of the initial rotation
	 * \param angle	the angle of rotation around axis
	 */
	CQuaternion(const CVector<3,T> &axis, T angle)
	{
		reset();
		rotatePost(axis, angle);
	}


	/**
	 * constructor
	 *
	 * \param angular_rotation	this parameter specifies the axis of rotation
	 * and also the amount of rotation by its length
	 */
	CQuaternion(const CVector<3,T> &angular_rotation)
	{
		T length = angular_rotation.getLength();
		if (length <= 0)
		{
			reset();
			return;
		}

		setRotation(angular_rotation * (1.0/length), length);
	}

	void applyAngularRotation(const CVector<3,T> &angular_rotation)
	{
		CQuaternion n;

		T length = angular_rotation.getLength();
		if (length <= 0)
		{
			n.reset();
			return;
		}

		this->rotatePost(angular_rotation * (1.0/length), length);
	}

	/**
	 * reset quaternion to 0 degree rotation
	 */
	void reset()
	{
		i = j = k = (T)0;
		w = (T)1.0;
	}

	/**
	 * apply an existing quaternion rotation to this one
	 */
	CQuaternion<T> &operator*=(const CQuaternion<T> &q)
	{
		T tw = w*q.w - i*q.i - j*q.j - k*q.k;
		T ti = w*q.i + i*q.w + j*q.k - k*q.j;
		T tj = w*q.j - i*q.k + j*q.w + k*q.i;
		T tk = w*q.k + i*q.j - j*q.i + k*q.w;

		i = ti;
		j = tj;
		k = tk;
		w = tw;

		return *this;
	}

	/**
	 * apply an existing quaternion rotation to this one
	 */
	CQuaternion<T> operator*(const CQuaternion<T> &q)
	{
		CQuaternion<T> p;

		p.w = w*q.w - i*q.i - j*q.j - k*q.k;
		p.i = w*q.i + i*q.w + j*q.k - k*q.j;
		p.j = w*q.j - i*q.k + j*q.w + k*q.i;
		p.k = w*q.k + i*q.j - j*q.i + k*q.w;

		return q;
	}

	/**
	 * apply the rotation from the "right side".
	 *
	 * thus first the point is rotated with the rotation specified by the
	 * parameters, then the point is rotated with the quaternion stored
	 * in this class before this method was called.
	 */
	CQuaternion<T> &rotatePost(const CVector<3,T> &axis, T angle)
	{
		T sinhalf = std::sin(angle*(T)-0.5);
		T ni = sinhalf*axis.data[0];
		T nj = sinhalf*axis.data[1];
		T nk = sinhalf*axis.data[2];
		T nw = std::cos(angle*(T)-0.5);

		T tw = nw*w - ni*i - nj*j - nk*k;
		T ti = ni*w + nw*i + nk*j - nj*k;
		T tj = nj*w - nk*i + nw*j + ni*k;
		T tk = nk*w + nj*i - ni*j + nw*k;

		i = ti;
		j = tj;
		k = tk;
		w = tw;

		return *this;
	}


	/**
	 * update the quaternion to first rotate with the parameters given to this method,
	 * then around the rotation stored in this quaternion
	 */
	CQuaternion<T> &rotatePre(const CVector<3,T> &axis, T angle)
	{
		T sinhalf = std::sin(angle*(T)-0.5);
		T ni = sinhalf*axis.data[0];
		T nj = sinhalf*axis.data[1];
		T nk = sinhalf*axis.data[2];
		T nw = std::cos(angle*(T)-0.5);

		T tw = nw*w - ni*i - nj*j - nk*k;
		T ti = nw*i + ni*w + nj*k - nk*j;
		T tj = nw*j - ni*k + nj*w + nk*i;
		T tk = nw*k + ni*j - nj*i + nk*w;

		i = ti;
		j = tj;
		k = tk;
		w = tw;

		return *this;
	}

	/**
	 * set the quaternion to a rotation around axis
	 * \param angle	given in radians
	 * \param axis	normalized axis of rotation with length 1
	 */
	CQuaternion<T> &setRotation(const CVector<3,T> &axis, T angle)
	{
		T sinhalf = std::sin(angle*(T)-0.5);
		i = sinhalf*axis.data[0];
		j = sinhalf*axis.data[1];
		k = sinhalf*axis.data[2];
		w = std::cos(angle*(T)-0.5);

		return *this;
	}

	/**
	 * set the quaternion to a rotation around axis
	 * \param angle	given in radians
	 * \param axis	normalized axis of rotation with length 1
	 */
	CQuaternion<T> &setRotation(
			const std::array<T,3> &axis,
			T angle
	)
	{
		T sinhalf = std::sin(angle*(T)-0.5);
		i = sinhalf*axis[0];
		j = sinhalf*axis[1];
		k = sinhalf*axis[2];
		w = std::cos(angle*(T)-0.5);

		return *this;
	}

	/**
	 * return the rotation matrix which describes the rotation of this
	 * quaternion.
	 */
	CMatrix3<T> getRotationMatrix()
	{
		return CMatrix3<T>(
				(T)1.0-(T)2.0*(j*j+k*k),	(T)2.0*(i*j+k*w),		(T)2.0*(i*k-j*w),
				(T)2.0*(i*j-k*w),			(T)1.0-2.0*(i*i+k*k),	(T)2.0*(j*k+i*w),
				(T)2.0*(i*k+j*w),			(T)2.0*(j*k-i*w),		(T)1.0-2.0*(i*i+j*j)
				);
	}

	/**
	 * due to the finite accuracy of the T numbers, it's good to
	 * normalize the quaternion sometimes. to normalize the quaternion
	 * values, just call this function.
	 */
	void normalize()
	{
		T inv_length = (T)1.0/std::sqrt(i*i + j*j + k*k + w*w);

		if (inv_length <= 0)
		{
			// reset to 0 degree rotation
			i = j = k = 0.0;
			w = 1.0;
			return;
		}

		i *= inv_length;
		j *= inv_length;
		k *= inv_length;
		w *= inv_length;
	}

};



template <class T>
inline
std::ostream&
operator<<(std::ostream &co, const CQuaternion<T> &q)
{
	return co	<< "[" << q.i << ", " << q.j << ", " << q.k << ", " << q.w << "]";
}


#endif
