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


#ifndef CGL_BUFFER_HPP
#define CGL_BUFFER_HPP

#include <libgl/core/GlError.hpp>
#include "libgl/incgl3.h"

/**
 * \brief draw a box in opengl3 sufficient for volume rendering
 */
class CGlBuffer
{
public:
	GLuint buffer;	///< OpenGL buffer id
	GLenum target;	///< default target

	/**
	 * create a buffer for specific target
	 */
	inline CGlBuffer(GLenum p_target = GL_ARRAY_BUFFER)
	{
		target = p_target;
		glGenBuffers(1, &buffer);
		CGlErrorCheck();
	}

	/**
	 * bind buffer
	 */
	inline void bind()
	{
		glBindBuffer(target, buffer);
		CGlErrorCheck();
	}


	/**
	 * bind buffer with target specification
	 */
	inline void bind(	GLenum p_target,
						GLuint index)
	{
		glBindBufferBase(p_target, index, buffer);
		CGlErrorCheck();
	}

	/**
	 * unbind buffer
	 */
	inline void unbind()
	{
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		CGlErrorCheck();
	}

	/**
	 * unbind buffer
	 */
	inline void unbind(	GLenum p_target,
						GLuint index)
	{
/* binding to zero buffer not described in opengl spec
		glBindBufferBase(p_target, index, 0);
*/
	}

	/**
	 * setup buffer with data
	 */
	inline void data(	GLsizeiptr size,		///< size of data
						const void *data,		///< pointer to data
						GLenum usage = GL_STATIC_DRAW	///< usage of data
	)
	{
		glBufferData(target, size, data, usage);
		CGlErrorCheck();
	}



	/**
	 * setup a part of the buffer with data - the buffer isn't resized
	 */
	inline void subData(	GLintptr offset,		///< offset within existing buffer
							GLsizeiptr size,		///< size of data
							const void *data		///< pointer to data
	)
	{
		glBufferSubData(target, offset, size, data);
		CGlErrorCheck();
	}


	/**
	 * resize buffer
	 */
	inline void resize(	GLsizeiptr size,				///< size of data
						GLenum usage = GL_STATIC_DRAW	///< usage of data
	)
	{
		glBufferData(target, size, nullptr, usage);
		CGlErrorCheck();
	}


	/**
	 * setup buffer with data and special target
	 */
	inline void dataWithTarget(
							GLenum p_target,		///< different target than default one
							GLsizeiptr size,		///< size of data
							const void *data,		///< pointer to data
							GLenum usage = GL_STATIC_DRAW	///< usage of data
	)
	{
		glBufferData(p_target, size, data, usage);
		CGlErrorCheck();
	}

	/**
	 * default deconstructor
	 */
	inline ~CGlBuffer()
	{
		glDeleteBuffers(1, &buffer);
		CGlErrorCheck();
	}
};

#endif
