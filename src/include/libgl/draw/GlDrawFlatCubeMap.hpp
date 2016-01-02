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


#ifndef C_GL_DRAW_FLAT_CUBE_MAP_HPP
#define C_GL_DRAW_FLAT_CUBE_MAP_HPP

#include "libgl/incgl3.h"
#include "libgl/core/GlProgram.hpp"
#include "libgl/core/GlTexture.hpp"
#include "libgl/CGlCubeMap.hpp"

#include "libmath/CGlSlMath.hpp"


/**
 * \brief render flat cube map texture for debugging purposes
 */
class CGlDrawFlatCubeMap
{
	CGlBuffer buffer;
	GlVertexArrayObject vao;

	GlProgram program;
	GlUniform pvm_matrix_uniform;

public:
	CError error;		///< error handler

	/**
	 * default constructor
	 */
	CGlDrawFlatCubeMap()
	{
		program.initVertFragShadersFromDirectory("draw/flat_cube_map");
		if (program.error())
		{
			std::cerr << program.error.getString() << std::endl;
			return;
		}
		program.bindAttribLocation(0, "vertex_position");

		program.link();
		if (program.error())
		{
			std::string infoLog;
			program.getInfoLog(infoLog);
			std::cerr << "info Log: linking: " << infoLog << std::endl;
			return;
		}

		program.bindAttribLocation(0, "vertex_position");
		program.bindAttribLocation(1, "vertex_texture_coord");

		program.setupUniform(pvm_matrix_uniform, "pvm_matrix");

		static const float vertices[4][4] = {
					{-1.0, -1.0, 0.0, 1.0},
					{ 1.0, -1.0, 0.0, 1.0},
					{-1.0,  1.0, 0.0, 1.0},
					{ 1.0,  1.0, 0.0, 1.0},
				};

		const float N = -1;
		const float P = 1;
		static const float texcoords[6][4][3] = {
				{	{P, N, N},	{P, N, P},	{P, P, N},	{P, P, P}	},	// X pos
				{	{N, N, P},	{N, N, N},	{N, P, P},	{N, P, N}	},	// X neg
				{	{N, P, N},	{P, P, N},	{N, P, P},	{P, P, P}	},	// Y pos
				{	{N, N, P},	{P, N, P},	{N, N, N},	{P, N, N}	},	// Y neg
				{	{P, N, P},	{N, N, P},	{P, P, P},	{N, P, P}	},	// Z pos
				{	{N, N, N},	{P, N, N},	{N, P, N},	{P, P, N}	},	// Z neg
		};

		vao.bind();

			buffer.bind();
			buffer.resize(sizeof(vertices)+sizeof(texcoords));

			buffer.subData(0, sizeof(vertices), vertices);
			glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
			glEnableVertexAttribArray(0);

			buffer.subData(sizeof(vertices), sizeof(texcoords), texcoords);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)sizeof(vertices));
			glEnableVertexAttribArray(1);

		vao.unbind();
	}


	/**
	 * render flat cube map
	 */
	void render(	CGlCubeMap &cGlCubeMap,	///< handler to cube map
					float scale = 1.0f		///< scale flat cube map with this factor
	)
	{

		cGlCubeMap.texture_cube_map.bind();

		CGlViewport cGlViewport;
		cGlViewport.saveState();

		GLSL::mat4 pmv = GLSL::ortho(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
		program.use();
		pvm_matrix_uniform.set(pmv);

		vao.bind();
		CGlErrorCheck();

		pmv = GLSL::ortho(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);

// line in texture coordinates (offset = vertices)
#define TEX_LINE(l)		(sizeof(GLfloat)*4*4 + (l)*4*3*sizeof(GLfloat))

		buffer.bind();
		// LEFT
		glViewport(0*scale, cGlCubeMap.depth*scale, cGlCubeMap.depth*scale, cGlCubeMap.height*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(1));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		// FRONT
		glViewport(cGlCubeMap.depth*scale, cGlCubeMap.depth*scale, cGlCubeMap.width*scale, cGlCubeMap.height*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(5));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		// RIGHT
		glViewport(cGlCubeMap.depth*scale+cGlCubeMap.width*scale, cGlCubeMap.depth*scale, cGlCubeMap.depth*scale, cGlCubeMap.height*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(0));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		// BACK
		glViewport(cGlCubeMap.depth*2*scale+cGlCubeMap.width*scale, cGlCubeMap.depth*scale, cGlCubeMap.width*scale, cGlCubeMap.height*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(4));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		// BOTTOM
		glViewport(cGlCubeMap.depth*scale, 0, cGlCubeMap.width*scale, cGlCubeMap.depth*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(3));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		// TOP
		glViewport(cGlCubeMap.depth*scale, cGlCubeMap.depth*scale + cGlCubeMap.height*scale, cGlCubeMap.width*scale, cGlCubeMap.depth*scale);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)TEX_LINE(2));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		cGlViewport.restoreState();

		vao.unbind();

		program.disable();

		CGlErrorCheck();

		glFlush();
	}
};

#endif
