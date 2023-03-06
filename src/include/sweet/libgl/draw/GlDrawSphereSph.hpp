/*
 * Copyright 2016 Martin Schreiber
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


#ifndef CGL_DRAW_SPHERE_SPH_HPP
#define CGL_DRAW_SPHERE_SPH_HPP

#include <sweet/libgl/core/GlBuffer.hpp>
#include <sweet/libgl/core/GlProgram.hpp>
#include <sweet/libgl/core/GlVertexArrayObject.hpp>
#include <sweet/libmath/CGlSlMath.hpp>

#if !SWEET_USE_SPHERE_SPECTRAL_SPACE
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not activated"
#endif

#include <sweet/core/sphere/SphereData_Config.hpp>

/**
 * \brief	create and render a sphere with polygons and
 * the latitudes placed with SPH ALP nodal points
 */
class GlDrawSphereSph
{
	GlBuffer vertex_buffer;
	GlBuffer texture_coord_buffer;
	GlBuffer index_buffer;
	GlVertexArrayObject vao;

	std::size_t triangle_strip_indices_count;


public:
	/**
	 * initialize OpenGL buffers to render a sphere with the given parameters
	 */
	void initBuffers(
			const sweet::SphereDataConfig *i_sphereDataConfig
	)
	{
		int hsegments = i_sphereDataConfig->physical_num_lon;
		int vsegments = i_sphereDataConfig->physical_num_lat;

		triangle_strip_indices_count = (hsegments+1)*2*vsegments;
		GLuint vertices_count = hsegments*(vsegments+1);

		GLfloat *vertices = new GLfloat[vertices_count*3];
		GLfloat *texture_coords = new GLfloat[vertices_count*2];
		GLuint *triangle_strip_indices = new GLuint[triangle_strip_indices_count];


		// create triangle indices
		GLuint *t = triangle_strip_indices;
		int count = 0;
		for (int dbeta = 0; dbeta < vsegments; dbeta++)
		{
			for (int dalpha = 0; dalpha < hsegments; dalpha++)
			{
				*t = dalpha+dbeta*hsegments;
				t++;
				*t = dalpha+(dbeta+1)*hsegments;
				t++;

				count += 2;
			}

			*t = dbeta*hsegments;
			t++;
			*t = (dbeta+1)*hsegments;
			t++;
			count += 2;
		}

		// create vertex data
		GLfloat *v = vertices;
		GLfloat *tc = texture_coords;
		int vcount = 0;
		for (int dbeta = 0; dbeta <= vsegments; dbeta++)
		{
			float y;

			if (dbeta == 0)
			{
				// first cell edge
				y = 0;
			}
			else if (dbeta == vsegments)
			{
				// last cell edge
				y = M_PI;
			}
			else
			{
				// compute average between two Gaussian latitude sampling points
				y = M_PI*0.5 - 0.5*(i_sphereDataConfig->lat[dbeta-1]+i_sphereDataConfig->lat[dbeta]);
			}


			for (int dalpha = 0; dalpha < hsegments; dalpha++)
			{
				float r = ((float)dalpha/(float)hsegments) * 2.0 * M_PI + M_PI;

				float scale = std::sin(y);
				v[0] = std::sin(r)*scale;
				v[1] = -std::cos(y);
				v[2] = -std::cos(r)*scale;
				v += 3;

				tc[0] = ((float)dalpha)/(float)hsegments;
				tc[1] = ((float)dbeta)/(float)vsegments;

				tc[0] = -tc[0];

#if 1
//				tc[0] += 0.5;
//				tc[0] += 0.25;
				tc[0] += 0.75;
				while (tc[0] > 1.0)
					tc[0]--;
				while (tc[0] < 0.0)
					tc[0]++;
#endif
				tc += 2;

				vcount++;
			}
		}

		GLuint *pt = triangle_strip_indices;
		for (std::size_t i = 0; i < triangle_strip_indices_count; i++)
		{
			if (*pt >= vertices_count)
				std::cerr << "ERROR" << std::endl;
			pt++;
		}

		vao.bind();

			vertex_buffer.bind();
			vertex_buffer.data(vertices_count*3*sizeof(GLfloat), vertices);

			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
			glEnableVertexAttribArray(0);

			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
			glEnableVertexAttribArray(1);

			texture_coord_buffer.bind();
			texture_coord_buffer.data(vertices_count*2*sizeof(GLfloat), texture_coords);
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
			glEnableVertexAttribArray(2);

			index_buffer.bind();
			index_buffer.data((triangle_strip_indices_count)*sizeof(GLuint), triangle_strip_indices);

		vao.unbind();

		delete[] vertices;
		delete[] texture_coords;
		delete[] triangle_strip_indices;
	}



	GlDrawSphereSph()	:
		vertex_buffer(GL_ARRAY_BUFFER),
		texture_coord_buffer(GL_ARRAY_BUFFER),
		index_buffer(GL_ELEMENT_ARRAY_BUFFER),
		triangle_strip_indices_count(0)
	{
	}



	/**
	 * initialize buffers to render sphere with given parameters
	 */
	void initSphere(
			const sweet::SphereDataConfig *i_sphereDataConfig
	)
	{
		initBuffers(i_sphereDataConfig);
	}



	/**
	 * render the sphere without using a GLSL program
	 */
	void renderWithoutProgram()
	{
		vao.bind();
			glDrawElements(GL_TRIANGLE_STRIP, triangle_strip_indices_count, GL_UNSIGNED_INT, 0);
		vao.unbind();

		CGlErrorCheck();
	}
};

#endif
