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


/*
 * CGlDrawFboQuad.hpp
 *
 *  Created on: Jan 30, 2010
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef CGLDRAWFBOQUAD_HPP_
#define CGLDRAWFBOQUAD_HPP_

#include "libgl/core/GlVertexArrayObject.hpp"

/**
 * draw an opengl quad with given viewport size
 *
 * this class is useful for fbo operations when quads
 * have to be drawn to use GL shaders for computations
 */
class CGlDrawFboQuad
{
private:
	GlVertexArrayObject vao_quad_vertices;
	GlBuffer quad_vertices_buffer;
	CGlViewport viewport;

public:
	inline CGlDrawFboQuad()
	{
		const GLfloat quad_vertices[4][4] = {
					{-1.0, -1.0, 0.0, 1.0},
					{ 1.0, -1.0, 0.0, 1.0},
					{-1.0,  1.0, 0.0, 1.0},
					{ 1.0,  1.0, 0.0, 1.0},
				};

		vao_quad_vertices.bind();

		quad_vertices_buffer.bind();
		quad_vertices_buffer.data(sizeof(quad_vertices), quad_vertices);

		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
		glEnableVertexAttribArray(0);

		vao_quad_vertices.unbind();
	}

	/**
	 * step 1) save the viewport state
	 */
	inline void draw_save_state()
	{
		viewport.saveState();
	}

	/**
	 * step 2) resize the viewport
	 */
	inline void draw_viewport_resize(GLsizei width, GLsizei height)
	{
		viewport.setSize(width, height);
	}


	/**
	 * step 3) initialization routines to draw a quad
	 */
	inline void draw_init_array()
	{
		vao_quad_vertices.bind();
	}

	/**
	 * step 4) draw a quad covering the whole viewport
	 */
	inline void draw()
	{
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	}

	/**
	 * step 5) finish drawing the quad
	 */
	inline void draw_finish_array()
	{
		vao_quad_vertices.unbind();
	}

	/**
	 * step 6) restore the previous viewport state
	 */
	inline void draw_restore_state()
	{
		viewport.restoreState();
	}

	/**
	 * convenient routine to draw unoptimized quad.
	 *
	 * use the other class functions if more than one quad has to be drawn.
	 */
	inline void draw_all(GLsizei width, GLsizei height)
	{
		draw_save_state();
		draw_init_array();
		draw_viewport_resize(width, height);
		draw();
		draw_finish_array();
		draw_restore_state();
	}
};

#endif /* CGLDRAWFBOQUAD_HPP_ */
