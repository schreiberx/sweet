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


#ifndef CGL_VERTEX_ARRAY_OBJECT_HPP
#define CGL_VERTEX_ARRAY_OBJECT_HPP

#include <sweet/libgl/core/GlError.hpp>
#include <sweet/libgl/incgl3.h>

/**
 * handle vertex array objects
 */
class GlVertexArrayObject
{
public:
	GLuint vbo_array;	///< OpenGL buffer id

	/**
	 * create a buffer for specific target
	 */
	inline GlVertexArrayObject()
	{
		glGenVertexArrays(1, &vbo_array);
		CGlErrorCheck();
	}

	/**
	 * bind buffer
	 */
	inline void bind()
	{
		glBindVertexArray(vbo_array);
		CGlErrorCheck();
	}


	/**
	 * unbind buffer
	 */
	inline void unbind()
	{
		glBindVertexArray(0);
		CGlErrorCheck();
	}


	/**
	 * default deconstructor
	 */
	inline ~GlVertexArrayObject()
	{
		glDeleteVertexArrays(1, &vbo_array);
		CGlErrorCheck();
	}
};

#endif
