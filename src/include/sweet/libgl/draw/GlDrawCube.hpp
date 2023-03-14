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



#ifndef CGL_DRAW_CUBE_HPP__
#define CGL_DRAW_CUBE_HPP__


#include <sweet/libgl/draw/GlDrawSphere.hpp>
#include <sweet/libgl/core/GlBuffer.hpp>
#include <sweet/libgl/core/GlProgram.hpp>
#include <sweet/libgl/core/GlVertexArrayObject.hpp>
#include <sweet/libmath/CGlSlMath.hpp>



/**
 * \brief	draw a box sufficient for volume rendering
 */
class GlDrawCube
{
	GlBuffer vertex_buffer;
	GlBuffer normal_buffer;
	GlBuffer texture_coord_buffer;

	GlVertexArrayObject vao;


public:
	GlDrawCube()	:
		vertex_buffer(GL_ARRAY_BUFFER),
		normal_buffer(GL_ARRAY_BUFFER)//,
//		texture_coord_buffer(GL_ARRAY_BUFFER)
	{
		/**
		 * initialize buffers
		 */
		/*
		 * vertices for cube drawn counterclockwise
		 * use quads to draw surfaces
		 */
		static const GLfloat vertices[36][3] = {

				// front
				{-1,1,1},
				{-1,-1,1},
				{1,1,1},
				{1,-1,1},
				{1,1,1},
				{-1,-1,1},

				// back
				{1,-1,-1},
				{-1,-1,-1},
				{1,1,-1},
				{-1,1,-1},
				{1,1,-1},
				{-1,-1,-1},


				// left
				{-1,-1,-1},
				{-1,-1,1},
				{-1,1,1},
				{-1,1,1},
				{-1,1,-1},
				{-1,-1,-1},

				// right
				{1,-1,-1},
				{1,1,-1},
				{1,1,1},
				{1,1,1},
				{1,-1,1},
				{1,-1,-1},


				// bottom
				{-1,-1,1},
				{-1,-1,-1},
				{1,-1,1},
				{1,-1,-1},
				{1,-1,1},
				{-1,-1,-1},

				// top
				{-1,1,-1},
				{-1,1,1},
				{1,1,1},
				{1,1,1},
				{1,1,-1},
				{-1,1,-1},
		};


		static const GLfloat normals[36][3] = {
				{0, 0, 1},
				{0, 0, 1},
				{0, 0, 1},
				{0, 0, 1},
				{0, 0, 1},
				{0, 0, 1},

				{0, 0, -1},
				{0, 0, -1},
				{0, 0, -1},
				{0, 0, -1},
				{0, 0, -1},
				{0, 0, -1},

				{-1, 0, 0},
				{-1, 0, 0},
				{-1, 0, 0},
				{-1, 0, 0},
				{-1, 0, 0},
				{-1, 0, 0},

				{1, 0, 0},
				{1, 0, 0},
				{1, 0, 0},
				{1, 0, 0},
				{1, 0, 0},
				{1, 0, 0},

				{0, -1, 0},
				{0, -1, 0},
				{0, -1, 0},
				{0, -1, 0},
				{0, -1, 0},
				{0, -1, 0},

				{0, 1, 0},
				{0, 1, 0},
				{0, 1, 0},
				{0, 1, 0},
				{0, 1, 0},
				{0, 1, 0},
		};
#if 0
		static const float texture_coords[24][3] = {
				{-1,-1,1},
				{-1,-1,-1},
				{-1,1,1},
				{-1,1,-1},
				{1,-1,1},
				{1,-1,1},
				{-1,-1,1},
				{-1,-1,-1},
				{-1,1,1},
				{-1,1,-1},
				{1,-1,1},
				{1,-1,1},
			};
#endif

		vao.bind();

			vertex_buffer.bind();
			vertex_buffer.data(sizeof(vertices), vertices);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
			glEnableVertexAttribArray(0);

			normal_buffer.bind();
			normal_buffer.data(sizeof(normals), normals);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
			glEnableVertexAttribArray(1);

//			texture_coord_buffer.bind();
//			texture_coord_buffer.data(sizeof(texture_coords), texture_coords);
//			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
//			glEnableVertexAttribArray(2);

		vao.unbind();
	}



	/**
	 * render the volume box without using a program
	 */
	void renderWithoutProgram()
	{
		vao.bind();
		glDrawArrays(GL_TRIANGLES, 0, 36);
		vao.unbind();

		CGlErrorCheck();
	}

};

#endif
