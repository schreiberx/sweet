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
 */

/*
 * CGlSlMath.hpp
 *
 *  Created on: Mar 5, 2010
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef CGLSLMATH_HPP_
#define CGLSLMATH_HPP_

#include <cmath>
#include <sweet/libmath/CMatrix.hpp>
#include <sweet/libmath/CVector.hpp>

/**
 * this class implements math operations similar to those available in the OpenGL shading language.
 *
 * the idea was taken from GLM (http://GLSL.g-truc.net/) which was not used to avoid the dependency.
 *
 * the only limitation for this class is the usage of float datatype instead of double datatype
 */
namespace GLSL
{
	typedef CMatrix2<float> mat2;	///< type definition for mat2
	typedef CMatrix3<float> mat3;	///< type definition for mat3
	typedef CMatrix4<float> mat4;	///< type definition for mat4

	typedef CVector<2,float> vec2;	///< type definition for vec2
	typedef CVector<3,float> vec3;	///< type definition for vec3
	typedef CVector<4,float> vec4;	///< type definition for vec4

	typedef CVector<2,int> ivec2;	///< type definition for vec2
	typedef CVector<3,int> ivec3;	///< type definition for vec3
	typedef CVector<4,int> ivec4;	///< type definition for vec4

	/**
	 * return the opengl translation matrix
	 */
	inline mat4 translate(
					vec3 t	///< translation
				)
	{
		mat4 m;
		m[0][3] = t[0];
		m[1][3] = t[1];
		m[2][3] = t[2];

		return m;
	}

	/**
	 * return the opengl translation matrix
	 */
	inline mat4 translate(	float tx,	///< x translation
							float ty,	///< y translation
							float tz	///< z translation
				)
	{
		mat4 m;
		m[0][3] = tx;
		m[1][3] = ty;
		m[2][3] = tz;

		return m;
	}


	/**
	 * return the opengl scaling matrix
	 */
	inline mat4 scale(	float sx,	///< x scale
						float sy,	///< y scale
						float sz	///< z scale
				)
	{
		mat4 m;
		m.loadIdentity();

		m[0][0] = sx;
		m[1][1] = sy;
		m[2][2] = sz;

		return m;
	}

	/**
	 * return the opengl rotation matrix
	 */
	inline mat4 rotate(	float angle,	///< rotation angle
						float x,		///< x component of rotation axis
						float y,		///< y component of rotation axis
						float z			///< z component of rotation axis
					)
	{
		// normalize rotation axis
		float div = 1.0/std::sqrt(x*x+y*y+z*z);

		x *= div;
		y *= div;
		z *= div;

		float mul = M_PI/180.0f;
		float c = std::cos(angle*mul);
		float s = std::sin(angle*mul);

		mat4 m;
		m[0][0] = x*x*(1-c)+c;
		m[0][1] = x*y*(1-c)-z*s;
		m[0][2] = x*z*(1-c)+y*s;
		m[0][3] = 0;

		m[1][0] = y*x*(1-c)+z*s;
		m[1][1] = y*y*(1-c)+c;
		m[1][2] = y*z*(1-c)-x*s;
		m[1][3] = 0;

		m[2][0] = x*z*(1-c)-y*s;
		m[2][1] = y*z*(1-c)+x*s;
		m[2][2] = z*z*(1-c)+c;
		m[2][3] = 0;

		m[3][0] = 0;
		m[3][1] = 0;
		m[3][2] = 0;
		m[3][3] = 1;

		return m;
	}


	/**
	 * return the inverse of 4x4 the matrix pm
	 */
	inline mat4 inverse(const mat4 &pm)
	{
		return pm.getInverse();
	}



	/**
	 * return the inverse transpose of 4x4 the matrix pm
	 */
	inline mat4 inverseTranspose(const mat4 &pm)
	{
		return pm.getInverseTranspose();
	}


	/**
	 * return the transpose of the 4x4 matrix pm
	 */
	inline mat4 transpose(const mat4 &pm)
	{
		return pm.getTranspose();
	}


	/**
	 * return the inverse of 3x3 the matrix pm
	 */
	inline mat3 inverse(const mat3 &pm)
	{
		return pm.getInverse();
	}


	/**
	 * return the inverse transpose of the 3x3 matrix pm
	 */
	inline mat3 inverseTranspose(const mat3 &pm)
	{
		return pm.getInverseTranspose();
	}


	/**
	 * return the transpose of the 3x3 matrix pm
	 */
	inline mat3 transpose(const mat3 &pm)
	{
		return pm.getTranspose();
	}

	/**
	 * return the length of the vector 3
	 */
	inline float length(const vec3 &v)
	{
		return v.getLength();
	}

	/**
	 * return the length^2 of the vector 3
	 */
	inline float length2(const vec3 &v)
	{
		return v.getLength2();
	}

	/**
	 * return the normalized vector v
	 */
	inline vec3 normalize(const vec3 &v)
	{
		return v.getNormalized();
	}

	/**
	 * return the absolute values of vector v
	 */
	inline vec3 abs(const vec3 &v)
	{
		return GLSL::vec3(std::abs(v.data[0]), std::abs(v.data[1]), std::abs(v.data[2]));
	}

	/**
	 * create orthogonal matrix
	 *
	 * http://www.opengl.org/sdk/docs/man/xhtml/glOrtho.xml
	 */
	inline mat4 ortho(	float left,
						float right,
						float bottom,
						float top,
						float nearval,
						float farval
					)
	{
		mat4 m;

		m[0][0] = 2.0/(right-left);
		m[0][1] = 0;
		m[0][2] = 0;
		m[0][3] = -(right+left)/(right-left);

		m[1][0] = 0;
		m[1][1] = 2.0/(top-bottom);
		m[1][2] = 0;
		m[1][3] = -(top+bottom)/(top-bottom);

		m[2][0] = 0;
		m[2][1] = 0;
		m[2][2] = -2.0/(farval - nearval);
		m[2][3] = -(farval+nearval)/(farval-nearval);

		m[3][0] = 0;
		m[3][1] = 0;
		m[3][2] = 0;
		m[3][3] = 1;

		return m;
	}

	/**
	 * return cross product
	 */
	inline vec3 crossProd(
						vec3 v1,
						vec3 v2
					)
	{
		return v1%v2;
	}

	/**
	 * create frustum matrix
	 *
	 * http://www.opengl.org/sdk/docs/man/xhtml/glFrustum.xml
	 */
	inline mat4 frustum(	float left,
							float right,
							float bottom,
							float top,
							float nearval,
							float farval
					)
	{
		mat4 m;

		m[0][0] = 2.0*nearval/(right-left);
		m[0][1] = 0;
		m[0][2] = (right+left)/(right-left);
		m[0][3] = 0;

		m[1][0] = 0;
		m[1][1] = 2.0*nearval/(top-bottom);
		m[1][2] = (top+bottom)/(top-bottom);
		m[1][3] = 0;

		m[2][0] = 0;
		m[2][1] = 0;
		m[2][2] = -(farval+nearval)/(farval-nearval);
		m[2][3] = -2.0*farval*nearval/(farval-nearval);

		m[3][0] = 0;
		m[3][1] = 0;
		m[3][2] = -1;
		m[3][3] = 0;

		return m;
	}

	/**
	 * create lookat matrix (same as gluLookAt)
	 */
	inline mat4 const lookAt(	const vec3 &eye,
								const vec3 &center,
								const vec3 &up
					)
	{
		vec3 f = (center-eye).getNormalized();
		vec3 u = up.getNormalized();

		// DAMN! the normalization is missing in the official glu manpage :/ :/ :/
		vec3 s = normalize(crossProd(f, u));
		u = crossProd(s, f);

		return mat4 (
					s[0],	s[1],	s[2],	0,
					u[0],	u[1],	u[2],	0,
					-f[0],	-f[1],	-f[2],	0,
					0,		0,		0,		1
				)*translate(-eye);
	}

	/**
	 * return pointer to raw data to use it with OpenGL. The ordering is ROW-MAJOR!!!
	 */
	inline float* value_ptr(const mat4 &m)
	{
		return (float*)m.matrix;
	}

	/**
	 * return pointer to raw data to use it with OpenGL. The ordering is ROW-MAJOR!!!
	 */
	inline float* value_ptr(const mat3 &m)
	{
		return (float*)m.matrix;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline float* value_ptr(const vec4 &v)
	{
		return (float*)v.data;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline float* value_ptr(const vec3 &v)
	{
		return (float*)v.data;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline float* value_ptr(const vec2 &v)
	{
		return (float*)v.data;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline const int* value_ptr(const ivec4 &v)
	{
		return (const int*)v.data;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline const int* value_ptr(const ivec3 &v)
	{
		return (const int*)v.data;
	}

	/**
	 * return pointer to raw data to use it with OpenGL
	 */
	inline const int* value_ptr(const ivec2 &v)
	{
		return (const int*)v.data;
	}
}

#endif /* CGLSLMATH_HPP_ */
