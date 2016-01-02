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


#ifndef CGL_CUBE_MAP_HPP
#define CGL_CUBE_MAP_HPP

#include "libgl/core/GlTexture.hpp"
#include "libgl/core/GlFbo.hpp"
#include "libgl/core/GlError.hpp"
#include "libgl/core/GlViewport.hpp"

#include "libmath/CGlSlMath.hpp"

/**
 * \brief callback interface to draw scene for cube map creation
 */
class CGlCubeMapCallbackClass
{
	virtual void createCallback(GLSL::mat4 &projection_matrix, GLSL::mat4 &view_matrix) = 0;

	virtual ~CGlCubeMapCallbackClass() {}
};

/**
 * \brief create cubemap sufficient for mirror textures or other fancy stuff
 */
class GlCubeMap
{

	GlViewport cGlViewport;

	GlFbo fbo_pos_x;		///< framebuffer to render the view aimed to positive x direction
	GlFbo fbo_neg_x;
	GlFbo fbo_pos_y;
	GlFbo fbo_neg_y;
	GlFbo fbo_pos_z;
	GlFbo fbo_neg_z;

public:
	CError error;	///< error handler

	GLuint width;	///< width of cubemap
	GLuint height;	///< height of cubemap
	GLuint depth;	///< depth of cubemap

	GlTexture texture_cube_map;	///< handler to cubemap texture
	GlTexture texture_depth;		///< handler to depth map texture


	/**
	 * resize cubemap to given parameters
	 *
	 * it's important to note, that cubemaps seem to be restricted to the same width, height and depth,
	 * but this class offers more flexibility which cannot be used with current drivers
	 */
	void resize(	GLuint p_width,
					GLuint p_height,
					GLuint p_depth
	)
	{
		width = p_width;
		height = p_height;
		depth = p_depth;

		GLuint depth_width = (depth > width ? depth : width);
		GLuint depth_height = (depth > height ? depth : height);

		texture_depth.bind();
		texture_depth.resize(depth_width, depth_height);
		texture_depth.unbind();

		/**
		 * initialize ALL cube map textures before binding to fbo
		 */
		texture_cube_map.bind();
		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_POSITIVE_X, depth, height);
		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, depth, height);

		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, width, depth);
		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, width, depth);

		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, width, height);
		texture_cube_map.resizeCubeMap(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, width, height);
		texture_cube_map.unbind();

		CGlErrorCheck();

// TODO: maybe we want to use only one depth buffer
#define prepare(alignment, textarget, a, b)										\
		fbo_##alignment.bind();													\
		fbo_##alignment.bindFramebufferTextureWithTarget(texture_cube_map, textarget);	\
		fbo_##alignment.bindDepthTexture(texture_depth);				\
		fbo_##alignment.unbind();												\
		if (fbo_##alignment.error())											\
		{																		\
			std::cerr << fbo_##alignment.error.getString();							\
			return;																\
		}																		\
		CGlErrorCheck();

		prepare(pos_x, GL_TEXTURE_CUBE_MAP_POSITIVE_X, depth, height);
		prepare(neg_x, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, depth, height);

		prepare(pos_y, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, width, depth);
		prepare(neg_y, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, width, depth);

		prepare(pos_z, GL_TEXTURE_CUBE_MAP_POSITIVE_Z, width, height);
		prepare(neg_z, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, width, height);
#undef prepare
	}


	/**
	 * default cube map constructor
	 */
	GlCubeMap()	:
		width(0),
		height(0),
		depth(0),

		texture_cube_map(GL_TEXTURE_CUBE_MAP),
		texture_depth(GL_TEXTURE_2D, GL_DEPTH_COMPONENT32, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE)
	{
	}

	/**
	 * create cube map
	 */
	void create(	void (*callback)(GLSL::mat4 &projection_matrix, GLSL::mat4 &view_matrix, void *user_data),	///< pointer to class with callback handlers
					const GLSL::vec3 &translate,	///< translation vector of model_matrix (last column)
					void *user_data = nullptr,			///< pointer to user data forwarded to callback function
					GLfloat near_plane = 1.0f,		///< near perspective plane
					GLfloat far_plane = 10.0f		///< far perspective plane
					)
	{
		// IMPORTANT!! we mirror the scene with the frustum values

		GLSL::mat4 projection_matrix = GLSL::frustum(-near_plane, near_plane, -near_plane, near_plane, near_plane, far_plane);
		GLSL::mat4 view_matrix = GLSL::rotate(180.0f, 0.0f, 0.0f, 1.0f);
		GLSL::mat4 p_view_matrix;

		cGlViewport.saveState();
		glClearColor(1, 1, 1, 0);

		// FRONT
		fbo_neg_z.bind();
		p_view_matrix = view_matrix;//*GLSL::rotate(180.0f, 0.0f, 0.0f, 1.0f);
		p_view_matrix *= GLSL::translate(translate);
		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		// BACK
		fbo_pos_z.bind();
		p_view_matrix = view_matrix*GLSL::rotate(180.0f, 0.0f, 1.0f, 0.0f);
		p_view_matrix *= GLSL::translate(translate);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		// LEFT
		fbo_neg_x.bind();
		p_view_matrix = view_matrix*GLSL::rotate(-90.0f, 0.0f, 1.0f, 0.0f);
		p_view_matrix *= GLSL::translate(translate);
		glViewport(0, 0, depth, height);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		// RIGHT
		fbo_pos_x.bind();
		p_view_matrix = view_matrix*GLSL::rotate(90.0f, 0.0f, 1.0f, 0.0f);
		p_view_matrix *= GLSL::translate(translate);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		// TOP
		fbo_pos_y.bind();
		p_view_matrix = view_matrix*GLSL::rotate(180.0f, 0.0f, 0.0f, 1.0f);
		p_view_matrix *= GLSL::rotate(-90.0f, 1.0f, 0.0f, 0.0f);
		p_view_matrix *= GLSL::translate(translate);
		glViewport(0, 0, depth, width);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		// BOTTOM
		fbo_neg_y.bind();
		p_view_matrix = view_matrix*GLSL::rotate(180.0f, 0.0f, 0.0f, 1.0f);
		p_view_matrix *= GLSL::rotate(90.0f, 1.0f, 0.0f, 0.0f);
		p_view_matrix *= GLSL::translate(translate);
		glViewport(0, 0, depth, width);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		callback(projection_matrix, p_view_matrix, user_data);

		fbo_neg_y.unbind();

		CGlErrorCheck();

		cGlViewport.restoreState();
	}

	~GlCubeMap()
	{

	}
};


#endif
