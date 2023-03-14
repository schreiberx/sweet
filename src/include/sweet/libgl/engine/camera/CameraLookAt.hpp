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

#ifndef __I_CAMERA_LOOK_AT_HPP__
#define __I_CAMERA_LOOK_AT_HPP__

#include "CameraCommon.hpp"
#include <sweet/libmath/CGlSlMath.hpp>

/**
 * this implements a simple "lookat" cammera
 */
class CCameraLookAt : public CameraCommon
{
	CVector<3,float> position;
public:
	/**
	 * setup the camera to look at a some point
	 * \param eye		the current view position
	 * \param center	the distant center point to look at
	 * \param up		the up vector of the current view position to create a unique view
	 */
	inline void lookAt(	const GLSL::vec3 &eye,
						const GLSL::vec3 &center,
						const GLSL::vec3 &up
				)
	{
		view_matrix = GLSL::lookAt(eye, center, up);
		position = center;
	}

	inline void computeMatrices()
	{

	}

	/**
	 * return the current camera position
	 */
	inline const CVector<3,float> getPosition()	const
	{
		return position;
	}
};

#endif
