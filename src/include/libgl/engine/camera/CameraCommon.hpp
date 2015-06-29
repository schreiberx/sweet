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

#ifndef __I_CAMERA_HPP__
#define __I_CAMERA_HPP__

#include <libmath/CGlSlMath.hpp>


/**
 * camera interface which has to be implemented by the camera classes
 */
template <typename T>
class CameraCommon
{
public:
	/**
	 * the view matrix of this camera computed by computeMatrices()
	 */
	CMatrix4<T> view_matrix;

	/**
	 * the projection matrix of this camera computed by computeMatrices()
	 */
	CMatrix4<T> projection_matrix;

	CameraCommon()
	{
		view_matrix.loadIdentity();
		projection_matrix.loadIdentity();
	}

	/**
	 * setup the projection matrix with a frustum
	 *
	 * TODO: maybe this function should be implemented by the camera implementations.
	 */
	inline void frustum(	T left,
							T right,
							T bottom,
							T top,
							T nearval,
							T farval
					)
	{
		projection_matrix = GLSL::frustum(left, right, bottom, top, nearval, farval);
	}

	inline void setViewMatrix(
			CMatrix4<T> &i_view_matrix
	)
	{
		view_matrix = i_view_matrix;
	}

	/**
	 * return the current camera position
	 */
	virtual const CVector<3,T>& getPosition() = 0;


	/**
	 * this method has to be implemented by the implementations.
	 *
	 * if this method is called, the camera class is assumed to setup
	 * the view_matrix and projection_matrix.
	 */
	virtual void computeMatrices() = 0;

	virtual ~CameraCommon()
	{

	}
};

#endif //__I_OBJECT_HPP__
