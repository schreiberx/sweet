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


#ifndef C_GL_DRAW_QUAD_HPP
#define C_GL_DRAW_QUAD_HPP

#include <libgl/core/GlBuffer.hpp>
#include <libgl/core/GlProgram.hpp>
#include <libgl/core/GlTexture.hpp>
#include <libgl/core/GlVertexArrayObject.hpp>
#include "libgl/incgl3.h"
#include "libmath/CGlSlMath.hpp"

/**
 * \brief	render a simple textured quad
 */
class GlDrawQuad
{

public:
	GlVertexArrayObject vao;
	GlBuffer buffer;

	GlDrawQuad()
	{
		static const float vertices[4][3] = {
					{-1.0, -1.0, -1.0},
					{ 1.0, -1.0, -1.0},
					{-1.0,  1.0, -1.0},
					{ 1.0,  1.0, -1.0},
				};

		static const float normals[4][3] = {
					{ 0.0,  0.0, -1.0},
					{ 0.0,  0.0, -1.0},
					{ 0.0,  0.0, -1.0},
					{ 0.0,  0.0, -1.0},
				};

		static const float texcoords[4][2] = {
					{0.0, 0.0},
					{1.0, 0.0},
					{0.0, 1.0},
					{1.0, 1.0},
				};

		vao.bind();

			buffer.bind();
			buffer.resize(sizeof(vertices)+sizeof(normals)+sizeof(texcoords));

			buffer.subData(0, sizeof(vertices), vertices);
			buffer.subData(sizeof(vertices), sizeof(normals), normals);
			buffer.subData(sizeof(vertices) + sizeof(normals), sizeof(texcoords), texcoords);

			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(0);

			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(vertices)));
			glEnableVertexAttribArray(1);

			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, (void*)(sizeof(vertices)+sizeof(normals)));
			glEnableVertexAttribArray(2);

		vao.unbind();

		CGlErrorCheck();
	}



	/*
	 * render the quad
	 */
	void renderWithoutProgram()
	{
		vao.bind();
			glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		vao.unbind();
		CGlErrorCheck();
	}


	~GlDrawQuad()
	{
	}
};

#endif
