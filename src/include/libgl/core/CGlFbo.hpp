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




#ifndef CGL_FBO_HPP
#define CGL_FBO_HPP

#include "libgl/incgl3.h"
#include "libgl/core/CGlError.hpp"
#include "libgl/core/CGlTexture.hpp"


/**
 * \brief class to handle framebuffer objects
 */
class CGlFbo
{
	GLuint fbo;		///< framebuffer object itself

public:
	CError error;	///< error handler

	/**
	 * create new fbo
	 * \param cTexture	texture to bind to framebuffer
	 */
	CGlFbo(CGlTexture &cTexture)
	{
		glGenFramebuffers(1, &fbo);					// create fbo

		bind();
		bindTexture(cTexture);
		unbind();
		CGlErrorCheck();
	}


	/**
	 * create new fbo
	 */
	CGlFbo()
	{
		glGenFramebuffers(1, &fbo);		// create fbo
		CGlErrorCheck();
	}

	/**
	 * free framebuffer OpenGL handler
	 */
	~CGlFbo()
	{
		glDeleteFramebuffers(1, &fbo);
	}


	/**
	 * bind texture to framebuffer and texture unit
	 * \param texture		texture to bind to fbo
	 * \param textureUnit	bind texture to attachment 'textureUnit'
	 * \param level			mipmap level
	 */
	void bindTexture(const CGlTexture &texture, GLint textureUnit = 0, GLint level = 0)
	{
		static const GLuint fbo_attachments[] = {
				GL_COLOR_ATTACHMENT0,
				GL_COLOR_ATTACHMENT1,
				GL_COLOR_ATTACHMENT2,
				GL_COLOR_ATTACHMENT3,
				GL_COLOR_ATTACHMENT4,
				GL_COLOR_ATTACHMENT5,
				GL_COLOR_ATTACHMENT6,
				GL_COLOR_ATTACHMENT7,
				GL_COLOR_ATTACHMENT8,
				GL_COLOR_ATTACHMENT9,
				GL_COLOR_ATTACHMENT10,
				GL_COLOR_ATTACHMENT11,
				GL_COLOR_ATTACHMENT12,
				GL_COLOR_ATTACHMENT13,
				GL_COLOR_ATTACHMENT14,
				GL_COLOR_ATTACHMENT15,
			};

		glFramebufferTexture2D(GL_FRAMEBUFFER, fbo_attachments[textureUnit], texture.target, texture.textureid, level);
		checkStatus();
		CGlErrorCheck();
	}



	/**
	 * bind a layer of a 3d or 2d array texture to a framebuffer
	 */
	void bindTextureLayer(	const CGlTexture &texture,	///< texture to bind to fbo
							GLint textureUnit,			///< bind texture to attachment 'textureUnit'
							GLint level,				///< mipmap level
							GLint layer					///< layer of 2d texture array or 3d texture
					)
	{
		static const GLuint fbo_attachments[] = {
				GL_COLOR_ATTACHMENT0,
				GL_COLOR_ATTACHMENT1,
				GL_COLOR_ATTACHMENT2,
				GL_COLOR_ATTACHMENT3,
				GL_COLOR_ATTACHMENT4,
				GL_COLOR_ATTACHMENT5,
				GL_COLOR_ATTACHMENT6,
				GL_COLOR_ATTACHMENT7,
				GL_COLOR_ATTACHMENT8,
				GL_COLOR_ATTACHMENT9,
				GL_COLOR_ATTACHMENT10,
				GL_COLOR_ATTACHMENT11,
				GL_COLOR_ATTACHMENT12,
				GL_COLOR_ATTACHMENT13,
				GL_COLOR_ATTACHMENT14,
				GL_COLOR_ATTACHMENT15,
			};

		glFramebufferTextureLayer(GL_FRAMEBUFFER, fbo_attachments[textureUnit], texture.textureid, level, layer);
		checkStatus();
		CGlErrorCheck();
	}

	/**
	 * bind a texture as depth buffer to the framebuffer
	 */
	void bindDepthTexture(	const CGlTexture &texture,
							GLint textureUnit = 0,
							GLint level = 0)
	{
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, texture.target, texture.textureid, level);
		checkStatus();
		CGlErrorCheck();
	}


	/**
	 * bind a depth layer of a 3d or 2d array texture to a framebuffer
	 */
	void bindDepthTextureLayer(	const CGlTexture &texture,	///< texture to bind to fbo
								GLint textureUnit,			///< bind texture to attachment 'textureUnit'
								GLint level,				///< mipmap level
								GLint layer					///< layer of 2d texture array or 3d texture
	)
	{
		glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, texture.textureid, level, layer);
		checkStatus();
		CGlErrorCheck();
	}
///////////////////////////////////////////

	/**
	 * bind texture to framebuffer and target
	 * \param texture		texture to bind to fbo
	 * \param textarget		target to bind fbo
	 */
	void bindFramebufferTextureWithTarget(CGlTexture &texture, GLenum textarget)
	{
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, textarget, texture.textureid, 0);
		checkStatus();
		CGlErrorCheck();
	}

	/**
	 * bind framebuffer, so we can work with it or initialize it
	 */
	void bind()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
		CGlErrorCheck();
	}

	/**
	 * unbind fbo
	 */
	void unbind()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		CGlErrorCheck();
	}

	/**
	 * check status and pipe message to 'error' in case of invalid fbo
	 */
	bool checkStatus()
	{
		GLint status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		switch (status)
		{
			case GL_FRAMEBUFFER_COMPLETE:
				return true;

			case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
				std::cerr << "FBO: incomplete attachment" << std::endl;
				return false;

			case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
				std::cerr << "FBO: incomplete missing attachment" << std::endl;
				return false;

			case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
				std::cerr << "FBO: incomplete draw buffer" << std::endl;
				return false;

			case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
				std::cerr << "FBO: incomplete read buffer" << std::endl;
				return false;

			case GL_FRAMEBUFFER_UNSUPPORTED:
				std::cerr << "FBO: unsupported" << std::endl;
				return false;

			default:
				std::cerr << "FBO: unknown error " << status << std::endl;
				return false;
		}
	}
};


#endif
