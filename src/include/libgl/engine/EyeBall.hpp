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



#ifndef CEYEBALL_HPP
#define CEYEBALL_HPP

#include "libmath/CVector.hpp"
#include "libmath/CMatrix.hpp"

/**
 * \brief	handle the mouse drag'n drop movement on a 3d ball
 *
 * the eye should rotate around the 3d ball and rotate around it's local x and y axis at the local position
 * the direction of the x and y axis in sphere coordinates is then given by the up / right coordinate
 *
 * the calling program just has to rotate everything around the cartesian x, y and z axis according to the
 * rotation degrees given by vec3f rotate.
 *
 * the rotation ordering is important: Z Y and then X
 *
 * TODO: this does not work as it should do but we can rotate somehow
 */
template <typename T>
class EyeBall
{
private:
public:
	// all vectors have length 1
	CVector<3,T> eye;		///< eye look direction - always to center
	CVector<3,T> up;		///< up vector
	CVector<3,T> right;		///< vector to the right

	CMatrix4<T> rotationMatrix;	///< rotation matrix for OpenGL matrix operations

	/**
	 * constructor: reset with default values
	 */
	EyeBall()
	{
		reset();
	}

	/**
	 * reset with default values, so that rotationMatrix is loaded with the identity matrix
	 */
	void reset()
	{
		eye = CVector<3,T>(0, 0, 1);
		up = CVector<3,T>(0, 1, 0);
		right = CVector<3,T>(1, 0, 0);

		reconstructRotationMatrix();
	}

	/**
	 * because we make slightly numerical changes, we always scale the base vectors
	 * to a length of 1 using this function
	 */
	void makeNumericStable()
	{
		right = up % eye;
		up = eye % right;
		eye = right % up;

		// TODO: use cross products to preserve 90 degree angles
		eye.normalize();
		eye.clamp1_1();

		up.normalize();
		up.clamp1_1();

		right.normalize();
		right.clamp1_1();
	}

	/**
	 * return the correct angle build by perpendicular lying vectors with length x and y
	 */
	T getATan(T x, T y, T extra_add = 0)
	{
		T v = y / x;
		T angle = atan(v);	// \in [-PI_2;PI_2]

		if (x < 0)
		{
			if (y < 0)	angle += -M_PI;
			else		angle += M_PI;
		}

		angle += extra_add;	// (0,0,1) is the origin eye position

		if (angle < -M_PI)	angle += 2*M_PI;
		else if (angle > M_PI)	angle -= 2*M_PI;

		return angle;
	}


	/**
	 * this function is called to recompute the rotation matrix which can then be used by OpenGL
	 */
	void reconstructRotationMatrix()
	{
		rotationMatrix[0][0] = right[0];
		rotationMatrix[1][0] = right[1];
		rotationMatrix[2][0] = right[2];
		rotationMatrix[3][0] = 0;

		rotationMatrix[0][1] = up[0];
		rotationMatrix[1][1] = up[1];
		rotationMatrix[2][1] = up[2];
		rotationMatrix[3][1] = 0;

		rotationMatrix[0][2] = eye[0];
		rotationMatrix[1][2] = eye[1];
		rotationMatrix[2][2] = eye[2];
		rotationMatrix[3][2] = 0;

		rotationMatrix[0][3] = 0;
		rotationMatrix[1][3] = 0;
		rotationMatrix[2][3] = 0; 
		rotationMatrix[3][3] = 1;
	}

	/**
	 * this function is called to recompute the inverse of the rotation matrix which can then be used by OpenGL
	 */
	void reconstructRotationMatrixInverse()
	{
		rotationMatrix[0][0] = right[0];
		rotationMatrix[1][0] = up[0];
		rotationMatrix[2][0] = eye[0];
		rotationMatrix[3][0] = 0;

		rotationMatrix[0][1] = right[1];
		rotationMatrix[1][1] = up[1];
		rotationMatrix[2][1] = eye[1];
		rotationMatrix[3][1] = 0;

		rotationMatrix[0][2] = right[2];
		rotationMatrix[1][2] = up[2];
		rotationMatrix[2][2] = eye[2];
		rotationMatrix[3][2] = 0;

		rotationMatrix[0][3] = 0;
		rotationMatrix[1][3] = 0;
		rotationMatrix[2][3] = 0;
		rotationMatrix[3][3] = 1;
	}

	/**
	 * rotate around given angles - the angles are given in degrees
	 * the degree parameters specify the rotation of a viewed object which lies in the center
	 *
	 * alpha:	rotation around "screen" y axis
	 * 		if positive, rotate rotate object clockwise
	 *
	 * beta:	rotation around "screen" x axis
	 * 		if positive, rotate object up
	 *
	 * gamma:	rotation around "screen" z axis
	 * 		if positive, rotate anti-clockwise (clockwise rotation of object)
	 */
	void rotate(	T alpha,
					T beta,
					T gamma
	)
	{
		CMatrix3<T> m;

		// alpha: rotate eye and right vector around up vector
		m.genRotation(-alpha*(M_PI/(T)180), up);
		eye = m*eye;
		right = m*right;

		// beta: rotate up around right vector
		m.genRotation(-beta*(M_PI/(T)180), right);
		eye = m*eye;
		up = m*up;

		// gamma: rotate up around right vector
		m.genRotation(-gamma*(M_PI/(T)180), eye);
		right = m*right;
		up = m*up;

		makeNumericStable();
	}

	/**
	 * rotate view around axis with angle alpha
	 */
	void rotate(	T alpha,			// rotation angle
					const CVector<3,T> &axis	// axis to rotate around
					)
	{
//		std::cout << axis << std::endl;
//		return;

		CMatrix4<float> rot;
		rot.genRotation(alpha*(M_PI/(T)180), axis);

		CMatrix3<float> invTranspMatrix(rot.getInverseTranspose());
//		CMatrix3<float> invTranspMatrix(rot);
//		CMatrix3<float> r(rotationMatrix);
//		CMatrix3<float> a(invTranspMatrix*r);
//		std::cout << a << std::endl;

		up = invTranspMatrix*up;
		right = invTranspMatrix*right;
//		eye = invTranspMatrix*eye;

		CVector<3,float> asdf = eye;
		asdf = asdf*(-1.0);
		eye = invTranspMatrix*asdf;
		eye = eye*(-1.0);
/*
		std::cout << std::endl;
		std::cout << "up: " << up << std::endl;
		std::cout << "right: " << right << std::endl;
		std::cout << "eye: " << eye << std::endl;
*/
		makeNumericStable();
	}
};


#endif
