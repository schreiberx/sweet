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


#ifndef C_GL_TEXTURE_HPP
#define C_GL_TEXTURE_HPP
/*
 * version 0.4
 *
 * changes:
 * 2009-12-17:
 *      updates to handle GL_TEXTURE_3D
 *
 * changes:
 * 2009-10-25:
 *      OpenGl 3.2 compatibility
 *
 * 2007-09-27:
 * 	added proxy texture to detect out of memory
 */

#include "libgl/incgl3.h"
#include <iostream>

// DON'T REMOVE!! HAS TO BE INCLUDED BEFORE SDL!!!
#include <cmath>
extern "C" {
	#include <SDL.h>
//	#include <SDL_image.h>
}


#include <libgl/core/GlError.hpp>



/**
 * \brief	general texture handler for 2d, 3d, rectangle and cubemap texture
 */
class GlTexture
{
public:
	int width;			///< width of texture
	int height;			///< height of texture
	int depth;			///< depth of texture (for 3d textures)

	GLuint textureid;	///< OpenGL texture handler

	GLint int_format;	///< internal texture format
	GLenum ext_format;	///< external texture format to load data from
	GLenum ext_type;	///< external texture format type

	GLenum target;		///< texture target, usually GL_TEXTURE_2D


	/**
	 * setup the texture types, parameters, internal formats, etc.
	 *
	 * this is only relevant for further usage, no OpenGL command is called and no conversions are done
	 */
	inline void setTextureParameters(
						GLenum p_target = GL_TEXTURE_2D,
						GLint p_int_format = GL_RGBA,
						GLenum p_ext_format = GL_RGBA,
						GLenum p_ext_type = GL_UNSIGNED_BYTE
					)
	{
		target = p_target;
		int_format = p_int_format;
		ext_format = p_ext_format;
		ext_type = p_ext_type;
	}



	/**
	 * initialize texture
	 */
	inline GlTexture(	GLenum p_target = GL_TEXTURE_2D,
						GLint p_int_format = GL_RGBA,
						GLenum p_ext_format = GL_RGBA,
						GLenum p_ext_type = GL_UNSIGNED_BYTE
					)
	{
		glGenTextures(1, &textureid);

		width = 0;
		height = 0;
		depth = 0;

		setTextureParameters(p_target, p_int_format, p_ext_format, p_ext_type);

		bind();

		setLinearInterpolation();

		setParam(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		setParam(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		setParam(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	}



	/**
	 * resize of initialize texture with the given size
	 *
	 * call BIND() before resizing!
	 */
	inline bool resize(int p_width, int p_height)
	{
		/**
		 * free texture data
		 */

#if 0
		if (width != 0)
			glTexImage2D(target, 0, int_format, 0, 0, 0, ext_format, ext_type, 0);

// proxy textures dont seem to work with current nvidia drivers
		GLint x, y;
		if (target == GL_TEXTURE_2D)
		{
			// try to create proxy texture
			glTexImage2D(GL_PROXY_TEXTURE_2D, 0, int_format, width, height, 0, ext_format, ext_type, 0);

			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);
		}
		else if (target == GL_TEXTURE_RECTANGLE)
		{
			// TODO
			std::cerr << "TODO" << std::endl;
			return false;

			// try to create proxy texture
			glTexImage2D(GL_PROXY_TEXTURE_2D, 0, int_format, width, height, 0, ext_format, ext_type, 0);

			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);
		}
		else
		{
			std::cerr << "ERROR: unsupported target texture" << std::endl;
			return false;
		}

		if (x == 0 && y == 0)
		{
			std::cerr << "ERROR: failed to create texture (maybe not enough memory)" << std::endl;
			return false;
		}
#endif

		width = p_width;
		height = p_height;

		glTexImage2D(target, 0, int_format, width, height, 0, ext_format, ext_type, 0);

		return true;
	}

	/**
	 * resize of initialize texture with the given size
	 *
	 * call BIND() before resizing!
	 */
	inline bool resizeMipMap(GLint p_width, GLint p_height, GLint level)
	{
		if (level == 0)
		{
			width = p_width;
			height = p_height;
		}

		glTexImage2D(target, level, int_format, p_width, p_height, 0, ext_format, ext_type, 0);
		return true;
	}

	/**
	 * resize cubemap
	 * target: GL_TEXTURE_CUBE_MAP_[POSITIVE/NEGATIVE]_[X/Y/Z]
	 */
	inline bool resizeCubeMap(GLenum target, GLuint p_width, GLuint p_height)
	{
		/**
		 * free texture data
		 */
		glTexImage2D(target, 0, int_format, 0, 0, 0, ext_format, ext_type, 0);

#if 0
			// try to create proxy texture
		glTexImage2D(GL_PROXY_TEXTURE_2D, 0, int_format, width, height, 0, ext_format, ext_type, 0);
		glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &x);
		glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &y);

		if (x == 0 && y == 0)
		{
			std::cerr << "ERROR: failed to create texture (maybe not enough memory)" << std::endl;
			return false;
		}
#endif

		width = p_width;
		height = p_height;

		glTexImage2D(target, 0, int_format, width, height, 0, ext_format, ext_type, 0);

		return true;
	}

	/**
	 * only update the width and height values in this class without modifying the texture itself
	 *
	 * \param p_width	new width value
	 * \param p_height	new height value
	 */
	inline void setSizeParameter(int p_width, int p_height)
	{
		width = p_width;
		height = p_height;
	}

	/**
	 * resize of initialize texture with the given size
	 */
	inline bool resize(int p_width, int p_height, int p_depth)
	{
		/**
		 * free texture data
		 */

#if 0
		if (width != 0)
			glTexImage3D(target, 0, int_format, 0, 0, 0, 0, ext_format, ext_type, 0);

// proxy textures dont seem to work with current nvidia drivers
		GLint x, y;
		if (target == GL_TEXTURE_3D)
		{
			// try to create proxy texture
			glTexImage3D(GL_PROXY_TEXTURE_3D, 0, int_format, width, height, 0, ext_format, ext_type, 0);

			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
		}
		else if (target == GL_TEXTURE_RECTANGLE)
		{
			// TODO
			std::cerr << "TODO" << std::endl;
			return false;

			// try to create proxy texture
			glTexImage2D(GL_PROXY_TEXTURE_3D, 0, int_format, width, height, 0, ext_format, ext_type, 0);

			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
		}
		else if (target == GL_TEXTURE_CUBE_MAP)
		{
			// TODO
			std::cerr << "TODO" << std::endl;
			return false;

			// try to create proxy texture
			glTexImage2D(GL_PROXY_TEXTURE_3D, 0, int_format, width, height, 0, ext_format, ext_type, 0);

			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_WIDTH, &x);
			glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D, 0, GL_TEXTURE_HEIGHT, &y);
		}
		else
		{
			std::cerr << "ERROR: unsupported target texture" << std::endl;
			return false;
		}

		if (x == 0 && y == 0)
		{
			std::cerr << "ERROR: failed to create texture (maybe not enough memory)" << std::endl;
			return false;
		}
#endif

		width = p_width;
		height = p_height;
		depth = p_depth;

		glTexImage3D(target, 0, int_format, width, height, depth, 0, ext_format, ext_type, 0);
		CGlErrorCheck();

		return true;
	}


	/**
	 * initialize texture from file specifying target
	 * this is useful for cubemap
	 */
	inline bool loadFromFile(
			const std::string &filename,
			GLuint target
	)
	{
#if 1
		std::cerr << "NOT SUPPORTED" << std::endl;
		exit(0);
#else
		SDL_Surface *surface = IMG_Load(filename.c_str());

		// could not load filename
		if (!surface)
		{
			std::cerr << "Cannot load image \"" << filename << "\""<< SDL_GetError() << std::endl;
			return false;
		}

		// work out what format to tell glTexImage2D to use...
		if (surface->format->BytesPerPixel == 3)
		{
			// RGB 24bit
			int_format = GL_RGB;
		}
		else if (surface->format->BytesPerPixel == 4)
		{
			// RGBA 32bit
			int_format = GL_RGBA;
		}
		else
		{
			SDL_FreeSurface(surface);
			std::cerr << "Image has to be 3 or 4 bytes per pixel \"" << filename << "\""<< std::endl;
			return false;
		}

		width = surface->w;
		height = surface->h;

		bind();

		resize(width, height);

		ext_format = int_format;
		ext_type = GL_UNSIGNED_BYTE;

		// this reads from the sdl surface and puts it into an opengl texture
		glTexImage2D(target, 0, int_format, surface->w, surface->h, 0, ext_format, ext_type, surface->pixels);
		SDL_FreeSurface(surface);

		setLinearInterpolation();
		setParam(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		setParam(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
#endif
		return true;
	}

	/**
	 * initialize texture with image
	 */
	inline bool loadFromFile(const std::string &filename)
	{
		return loadFromFile(filename, GL_TEXTURE_2D);
	}

	/**
	 * create texture from file
	 */
	inline GlTexture(const char *filename)
	{
		loadFromFile(filename);
	}

	/**
	 * bind texture to texture unit
	 */
	inline void bind()
	{
		glBindTexture(target, textureid);
	}

	/**
	 * unbind texture target
	 */
	inline void unbind()
	{
		glBindTexture(target, 0);
	}

	/**
	 * set texture GLfloat parameter
	 */
	inline void setParam(GLenum pname, GLfloat param)
	{
		glTexParameterf(target, pname, param);
	}

	/**
	 * set texture GLint parameter
	 */
	inline void setParam(GLenum pname, GLint param)
	{
		glTexParameteri(target, pname, param);
	}


	/**
	 * set texture GLfloat[4] parameter
	 */
	inline void setParamfv(GLenum pname, const GLfloat *param)
	{
		glTexParameterfv(target, pname, param);
	}


	/**
	 * set texture GLint[4] parameter
	 */
	inline void setParamiv(GLenum pname, const GLint *param)
	{
		glTexParameteriv(target, pname, param);
	}

	/**
	 * set texture GLuint[4] parameter
	 */
	inline void setParamIuiv(GLenum pname, const GLuint *param)
	{
		glTexParameterIuiv(target, pname, param);
	}

	/**
	 * set texture GLenum parameter
	 */
	inline GLfloat getParami(GLenum pname)
	{
		GLfloat param;
		glGetTexParameterfv(target, pname, &param);
		return param;
	}

	/**
	 * return texture GLenum parameter
	 */
	inline GLint getParam(GLenum pname)
	{
		GLint param;
		glGetTexParameteriv(target, pname, &param);
		return param;
	}

	/**
	 * return texture id
	 */
	inline GLuint getTextureId()
	{
		return textureid;
	}

	/**
	 * return internal texture format
	 */
	inline GLenum getFormat()
	{
		return int_format;
	}

	/**
	 * copy the data pointed by data to the texture
	 *
	 * \param data	pointer to new texture data (float array with only one [red] component)
	 */
	inline void setData(void *data)
	{
		bind();
		switch(target)
		{
			case GL_TEXTURE_2D:
			case GL_TEXTURE_RECTANGLE:
				glTexImage2D(target, 0, int_format, width, height, 0, ext_format, ext_type, data);
				break;

			case GL_TEXTURE_3D:
				glTexImage3D(GL_TEXTURE_3D, 0, int_format, width, height, depth, 0, ext_format, ext_type, data);
				break;

			default:
				std::cerr << "unsupported target (in setData())" << std::endl;
				break;
		}
	}

	/**
	 * copy the data pointed by data to the texture
	 * \param data	pointer to new texture data (float array with only one [red] component)
	 * \param level	mipmap level
	 */
	inline void setDataMipMap(void *data, GLint level = 0)
	{
		GLint nwidth = width >> level;
		GLint nheight = height >> level;
		if (nwidth == 0)	nwidth = 1;
		if (nheight == 0)	nheight = 1;

		bind();
		switch(target)
		{
			case GL_TEXTURE_2D:
			case GL_TEXTURE_RECTANGLE:
				glTexImage2D(target, level, int_format, nwidth, nheight, 0, ext_format, ext_type, data);
				break;

			case GL_TEXTURE_3D:
				glTexImage3D(GL_TEXTURE_3D, level, int_format, nwidth, nheight, depth, 0, ext_format, ext_type, data);
				break;

			default:
				std::cerr << "unsupported target (in setData())" << std::endl;
				break;
		}
	}

	/**
	 * copy the texture data to the memory location referred by data
	 * 
	 * \param data	pointer to store the texture data
	 * \param level	mipmap level
	 */
	inline void getData(void *data, GLint level = 0)
	{
		bind();
		glGetTexImage(target, level, ext_format, ext_type, data);
	}


	/**
	 * setup quickly nearest neighbor filtering
	 */
	inline void setNearestNeighborInterpolation()
	{
		bind();
		setParam(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		setParam(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	}

	/**
	 * setup quickly nearest neighbor filtering
	 */
	inline void setLinearInterpolation()
	{
		bind();
		setParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		setParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}

	inline ~GlTexture()
	{
		glDeleteTextures(1, &textureid);
	}

	/**
	 * print red channel of a texture with internal type unsigned integer
	 * \param level	mipmap level to read data
	 */
	void printR32UI(GLint level = 0)
	{
		int nwidth = this->width >> level;
		int nheight = this->height >> level;

		if (nwidth == 0)	nwidth = 1;
		if (nheight == 0)	nheight = 1;

		std::cout << nwidth << " x " << nheight << std::endl;;

		GLuint *data = new GLuint[nwidth*nheight];
		bind();
		glGetTexImage(GL_TEXTURE_2D, level, GL_RED_INTEGER, GL_UNSIGNED_INT, data);
		glFinish();

		for(int y = 0; y < nheight; y++)
		{
			for(int x = 0; x < nwidth; x++)
				std::cout << data[x+nwidth*y] << " ";
			std::cout << std::endl;
		}

		unbind();

		delete[] data;
	}

	/**
	 * print red and green channel of a texture with internal type unsigned integer
	 * \param level	mipmap level to read data
	 */
	void printRG32UI(GLint level = 0)
	{
		int nwidth = this->width >> level;
		int nheight = this->height >> level;

		if (nwidth == 0)	nwidth = 1;
		if (nheight == 0)	nheight = 1;

		std::cout << nwidth << " x " << nheight << std::endl;;

		GLuint *data = new GLuint[nwidth*nheight*2];
		bind();
		glGetTexImage(GL_TEXTURE_2D, level, GL_RG_INTEGER, GL_UNSIGNED_INT, data);
		glFinish();

		for(int y = 0; y < nheight; y++)
		{
			for(int x = 0; x < nwidth; x++)
				std::cout << "(" << data[(x+nwidth*y)*2] << "," << data[(x+nwidth*y)*2+1] << ") ";
			std::cout << std::endl;
		}

		unbind();

		delete[] data;
	}

	/**
	 * print rgba channels of a texture with internal type unsigned integer
	 * \param level	mipmap level to read data
	 */
	void printRGBA32UI(GLint level = 0)
	{
		int nwidth = this->width >> level;
		int nheight = this->height >> level;

		if (nwidth == 0)	nwidth = 1;
		if (nheight == 0)	nheight = 1;

		std::cout << nwidth << " x " << nheight << std::endl;;

		GLuint *data = new GLuint[nwidth*nheight*4];
		bind();
		glGetTexImage(GL_TEXTURE_2D, level, GL_RGBA_INTEGER, GL_UNSIGNED_INT, data);
		glFinish();

		for(int y = 0; y < nheight; y++)
		{
			for(int x = 0; x < nwidth; x++)
				std::cout << "(" << data[(x+nwidth*y)*4] << "," << data[(x+nwidth*y)*4+1] << "," << data[(x+nwidth*y)*4+2] << "," << data[(x+nwidth*y)*4+3] << ") ";
			std::cout << std::endl;
		}

		unbind();

		delete[] data;
	}


	/**
	 * print rgba channels of a texture with internal type float
	 * \param level	mipmap level to read data
	 */
	void printRGBA32F(GLint level = 0)
	{
		int nwidth = this->width >> level;
		int nheight = this->height >> level;

		if (nwidth == 0)	nwidth = 1;
		if (nheight == 0)	nheight = 1;

		std::cout << nwidth << " x " << nheight << std::endl;;

		GLfloat *data = new GLfloat[nwidth*nheight*4];
		bind();
		glGetTexImage(GL_TEXTURE_2D, level, GL_RGBA, GL_FLOAT, data);
		glFinish();

		for(int y = 0; y < nheight; y++)
		{
			for(int x = 0; x < nwidth; x++)
				std::cout << "(" << data[(x+nwidth*y)*4] << "," << data[(x+nwidth*y)*4+1] << "," << data[(x+nwidth*y)*4+2] << "," << data[(x+nwidth*y)*4+3] << ") ";
			std::cout << std::endl;
		}

		unbind();

		delete[] data;
	}
};

#endif
