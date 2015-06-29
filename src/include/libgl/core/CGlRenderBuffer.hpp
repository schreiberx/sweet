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


#ifndef C_GL_RENDERBUFFER_HPP
#define C_GL_RENDERBUFFER_HPP
/*
 * version 0.1
 *
 * changes:
 * 2010-01-30:
 *      created
 *
 */

#include "libgl/incgl3.h"
#include "libgl/core/CGlError.hpp"
#include <iostream>

/**
 * \brief	handler for opengl render buffers
 */
class CGlRenderbuffer
{
public:
	CError error;		///< error handler
	int width;			///< width of texture
	int height;			///< height of texture

	GLuint renderbufferid;	///< OpenGL render buffer handler
	GLenum target;			///< texture target, usually GL_TEXTURE_2D
	GLenum internalformat;	///< internalformat

	/**
	 * initialize renderbuffer
	 */
	inline CGlRenderbuffer(GLenum p_target = GL_RENDERBUFFER, GLenum p_internalformat = GL_RGBA)
	{
		glGenRenderbuffers(1, &renderbufferid);

		target = p_target;
		internalformat = p_internalformat;
		width = 0;
		height = 0;
	}


	/**
	 * resize renderbuffer
	 *
	 * call bind() before resizing!
	 */
	inline bool resize(int p_width, int p_height)
	{
		width = p_width;
		height = p_height;

		glRenderbufferStorage(target, internalformat, width, height);

		return true;
	}

	/**
	 * bind texture to texture unit
	 */
	inline void bind()
	{
		glBindRenderbuffer(target, renderbufferid);
	}

	/**
	 * unbind texture target
	 */
	inline void unbind()
	{
		glBindRenderbuffer(target, 0);
	}

	/**
	 * return renderbuffer id
	 */
	inline GLuint operator()()
	{
		return renderbufferid;
	}

	inline ~CGlRenderbuffer()
	{
		glDeleteRenderbuffers(1, &renderbufferid);
	}
};

#endif
