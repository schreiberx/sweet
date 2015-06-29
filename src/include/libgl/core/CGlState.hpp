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




#ifndef CGL_STATE_HPP_
#define CGL_STATE_HPP_

#include "libgl/incgl3.h"

/**
 * \brief enable an OpenGL state
 *
 * this class can be used to avoid resetting a state after the usage.
 *
 * initializing this class within a {} "sub program" frees the class at the "}" curly bracket.
 * this makes sure, that the old state is reset to it's old state.
 */
class CGlStateEnable
{
	GLuint flag;		///< old state
	bool not_modified;

public:

	/**
	 * enable an OpenGL state if not already enabled
	 */
	inline CGlStateEnable(GLuint p_flag)
	{
		not_modified = (glIsEnabled(p_flag) == GL_TRUE);
		if (not_modified)	return;
		flag = p_flag;
		glEnable(flag);
	}

	/**
	 * restore old state
	 */
	inline ~CGlStateEnable()
	{
		if (not_modified)	return;
		glDisable(flag);
	}
};


/**
 * \brief disable an OpenGL state
 *
 * (opposite of CGlStateEnable)
 */
class CGlStateDisable
{
	GLuint flag;
	bool not_modified;

public:

	/**
	 * disable an OpenGL state if not already disabled
	 */
	inline CGlStateDisable(GLuint p_flag)
	{
		not_modified = (glIsEnabled(p_flag) == GL_FALSE);
		if (not_modified)	return;
		flag = p_flag;
		glDisable(flag);
	}

	/**
	 * restore old state
	 */
	inline ~CGlStateDisable()
	{
		if (not_modified)	return;
		glEnable(flag);
	}
};

#endif /* CGLSTATE_HPP_ */
