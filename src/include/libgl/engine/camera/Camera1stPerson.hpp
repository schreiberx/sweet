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

#ifndef __I_CAMERA_1ST_PERSON_HPP__
#define __I_CAMERA_1ST_PERSON_HPP__

#include <libgl/engine/camera/CameraCommon.hpp>
#include <array>
#include <libmath/CGlSlMath.hpp>
#include <libmath/CQuaternion.hpp>

/**
 * this implements a 1st person free flight cammera.
 *
 * e. g. its possible to freely rotate, move to the right, left according
 * to the current cammera angle.
 */
template <typename T>
class Camera1stPerson :
	public CameraCommon<T>
{
	CQuaternion<T> rotation_quaternion;
	T angle_x, angle_y, angle_z;

public:
	CVector<3,T> position;


	/**
	 * reset to origin and zero rotation
	 */
	inline void reset()
	{
		for (int i : {0,1,2})
			position[i] = 0;

		rotation_quaternion.reset();
		CameraCommon<T>::view_matrix.loadIdentity();
		angle_x = angle_y = angle_z = 0;
	}

	Camera1stPerson()
	{
		reset();
	}

	/**
	 * return the current camera position
	 */
	inline const CVector<3,T>& getPosition()
	{
		return position;
	}

	/**
	 * move the player relative to the current view
	 */
	inline void moveRelative(const CVector<3,T> &movement)
	{
		position += movement*CMatrix3<T>(CameraCommon<T>::view_matrix);
	}

	/**
	 * set the position of the player
	 */
public:
	inline void setPosition(
			const CVector<3,T> &i_position
	)
	{
		position = i_position;
	}

	/**
	 * rotate the player using directly a quaternion
	 */
	inline void rotate(const CQuaternion<T> &quaternion)
	{
		rotation_quaternion *= quaternion;
	}

	/**
	 * reset the rotation
	 */
	inline void resetRotation()
	{
		rotation_quaternion.reset();
	}

	/**
	 * rotate the player around axes
	 */
	inline void rotate(T p_angle_x, T p_angle_y, T p_angle_z)
	{
		angle_x += p_angle_x;
		angle_y += p_angle_y;
		angle_z += p_angle_z;
	}

	inline void computeMatrices()
	{
		rotation_quaternion.setRotation(CVector<3,T>(0, 1, 0), angle_y);
		rotation_quaternion.rotatePost(CVector<3,T>(1, 0, 0), angle_x);

		CameraCommon<T>::view_matrix = CMatrix4<T>(rotation_quaternion.getRotationMatrix())*GLSL::translate(-position);
	}
};

#endif
