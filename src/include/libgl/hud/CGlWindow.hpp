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
 * CGlWindow.hpp
 *
 *  Created on: Mar 22, 2010
 *      Author: martin
 */

#ifndef CGLWINDOW_HPP_
#define CGLWINDOW_HPP_

#include "libmath/CGlSlMath.hpp"
#include "libgl/core/CGlViewport.hpp"
#include "libgl/core/CGlState.hpp"
#include "libgl/draw/CGlDrawQuad.hpp"

class CGlWindow
{
	CGlViewport viewport;
	CGlProgram program;
	CGlUniform background_color_uniform;

	CGlDrawQuad draw_quad;

	GLSL::vec4 background_color;

public:
	CError error;

	GLSL::ivec2 pos;		///< position of window
	GLSL::ivec2 size;		///< size of window

	bool inWindow(int mouse_x, int mouse_y)
	{
		return (mouse_x >= pos[0] && mouse_y >= pos[1] && mouse_x < size[0]+pos[0] && mouse_y < size[1]+pos[1]);
	}

	CGlWindow()
	{
		// setup shader program
		program.initVertFragShadersFromDirectory("window");
		if (program.error())
		{
			std::cerr << program.error.getString() << std::endl;
			return;
		}

		program.link();
		if (program.error())
		{
			std::cerr << "info Log: linking: " << program.getInfoLog() << std::endl;
			return;
		}

		program.use();
		program.bindAttribLocation(0, "vertex_attrib");
		program.setupUniform(background_color_uniform, "background_color");
		program.disable();

		setBackgroundColor(GLSL::vec4(74.0/256.0, 96.0/256.0, 150.0/256.0, 0.7));
	}

	/**
	 * set the window background color
	 */
	void setBackgroundColor(const GLSL::vec4 &p_background_color)
	{
		program.use();
		background_color_uniform.set(p_background_color);
		program.disable();
	}

	void setSize(const GLSL::ivec2 &p_size)
	{
		size = p_size;
	}
	void setPosition(const GLSL::ivec2 &p_pos)
	{
		pos = p_pos;
	}

	/**
	 * this method has to be called, when something should be rendered to the window
	 */
	void startRendering()
	{
		viewport.saveState();
		viewport.set(pos[0], pos[1], size[0], size[1]);

		{
			CGlStateDisable depth_test(GL_DEPTH_TEST);
			CGlStateEnable blend(GL_BLEND);

			// enable blending
			glBlendEquation(GL_FUNC_ADD);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			program.use();
			draw_quad.render();
			program.disable();
		}
	}

	/**
	 * this method has to be called, when rendering to the window is finished
	 */
	void finishRendering()
	{
		viewport.restoreState();
	}
};

#endif /* CGLWINDOW_HPP_ */
