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
 * CGlFreeType.hpp
 *
 *  Created on: Mar 21, 2010
 *      Author: martin
 */

#ifndef CGLFREETYPE_HPP_
#define CGLFREETYPE_HPP_

#include <algorithm>
#include <ft2build.h>
#include FT_FREETYPE_H



#include <cmath>

#include "libgl/core/CGlTexture.hpp"
#include "libgl/core/CGlBuffer.hpp"
#include "libgl/core/CGlVertexArrayObject.hpp"
#include "libgl/core/CGlProgram.hpp"
#include "libgl/core/CGlState.hpp"

#include <string.h>

class CGlFreeType
{
private:
	int advance_x[256];

	CGlTexture font_texture;
	CGlProgram program;

	CGlUniform glyph_size_uniform;
	CGlUniform inv_viewport_size_uniform;
	CGlUniform color_uniform;

	int newline_x_position;	///< x position for newline

	int *vertex_attrib;		///< allocated memory for vertex attributes
	CGlBuffer buffer;
	CGlVertexArrayObject vao;

	bool valid;

public:
	int char_width = 0;			///< width of each char in texture
	int char_height = 0;		///< height of each char in texture
	int char_origin_left = 0;	///< base (0,0) point of each char in texture - this is the distance of the char origin to the top border in texels
	int char_origin_top = 0;	///< base (0,0) point of each char in texture

	std::array<float,2> stylus_pos;	///< current text drawing position
	int font_size;		///< font size in pixels

	int max_length;		///< maximum length of string

	CGlFreeType(
			float p_font_size = 12,
			const char *p_font_file = nullptr
	)	:
			font_texture(GL_TEXTURE_2D, GL_RED, GL_RED, GL_UNSIGNED_BYTE),
			font_size(p_font_size),
			max_length(1024)
	{
		vertex_attrib = new GLint[4*4*max_length];	// allow a maximum of 1024 chars to be rendered at once

		// setup shader program
		if (!program.initVertFragShadersFromDirectory("free_type"))
			return;

		if (!program.link())
		{
			std::cerr << "info Log: linking: " << program.getInfoLog() << std::endl;
			return;
		}

		program.use();
		program.bindAttribLocation(0, "vertex_attrib");

		program.setupUniform(glyph_size_uniform, "glyph_size");
		program.setupUniform(inv_viewport_size_uniform, "inv_viewport_size");
		program.setupUniform(color_uniform, "color");
		program.disable();

		setColor(GLSL::vec3(1,1,1));

		newline_x_position = 0;

		vao.bind();
			buffer.bind();
			buffer.resize(max_length*4*4*sizeof(GLint), GL_DYNAMIC_DRAW);
			glVertexAttribIPointer(0, 4, GL_INT, 0, 0);
			glEnableVertexAttribArray(0);
		vao.unbind();

		valid = false;

		CGlErrorCheck();
	}

	~CGlFreeType()
	{
		delete[] vertex_attrib;

		CGlErrorCheck();
	}

	bool trySetupFont(
			FT_Library &p_ft_library,
			const char *p_font_file,
			FT_Face &ret_ft_face)
	{
		FT_Error ft_error;			///< handle to error handling
		ft_error = FT_New_Face(p_ft_library, p_font_file, 0, &ret_ft_face);
		if (ft_error)
		{
			std::cerr << "unable to load font file " << p_font_file << std::endl;
			return false;
		}
		return true;
	}


	void loadFont(
			int i_font_size = 14,
			bool i_bold = false,
			int i_char_spacing = 0
	)
	{
		FT_Error ft_error;			///< handle to error handling

		font_size = i_font_size;

		FT_Library ft_library;		///< handle to freetype library
		ft_error = FT_Init_FreeType(&ft_library);

		if (ft_error)
		{
			std::cerr << "an error occurred during freetype library initialization" << std::endl;
			return;
		}

		FT_Face	ft_face;			///< handle to face object
		bool no_bold_font_found = false;
		if (i_bold)
		{
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/ubuntu-font-family/Ubuntu-B.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/msttcorefonts/arialbd.ttf", ft_face))
			{
				no_bold_font_found = true;
			}

		}
		else
		{
			no_bold_font_found = false;
		}

		if (no_bold_font_found)
		{
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/ubuntu-font-family/Ubuntu-R.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/freefont/FreeSans.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/msttcorefonts/arial.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/corefonts/arial.ttf", ft_face))
			if (!trySetupFont(ft_library, "/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans.ttf", ft_face))
			{
				std::cerr << "ERROR: Could not find valid font file... aborting" << std::endl;
				return;
			}
		}

		ft_error = FT_Set_Pixel_Sizes(ft_face, 0, font_size);
		if (ft_error)
		{
			std::cerr << "unable to set font size" << std::endl;
			FT_Done_Face(ft_face);
			return;
		}


		/*
		 * compute borders and origin point
		 */
		FT_GlyphSlot ft_slot = ft_face->glyph;

		char_width = 0;
		char_height = 0;

		char_origin_left = 0;	// char origin within texture
		char_origin_top = 0;	// char origin within texture

		// ft_slot->bitmap_top -> bitmap rows from top to glyph origin
		// ft_slot->bitmap_left -> bitmap colums from left to glyph origin

		int dright = 0, dtop = 0;	// distance to bitmap borders from origin

		for (int i = 0; i < 256; i++)
		{
			unsigned int glyph_index = FT_Get_Char_Index(ft_face, (unsigned char)i);

			// load glyph to slot (render it)
			ft_error = FT_Load_Glyph(ft_face, glyph_index, FT_LOAD_RENDER);
			if (ft_error)
			{
				std::cerr << "unable to load glyph" << std::endl;
				FT_Done_Face(ft_face);
				return;
			}

			FT_Bitmap bitmap = ft_slot->bitmap;

			dright = std::max((int)(bitmap.width - ft_slot->bitmap_left), dright);
			char_origin_left = std::max(ft_slot->bitmap_left, char_origin_left);

			dtop = std::max((int)(bitmap.rows - ft_slot->bitmap_top), dtop);
			char_origin_top = std::max(ft_slot->bitmap_top, char_origin_top);
		}

		char_width = dright + char_origin_left;
		char_height = dtop + char_origin_top;

		font_texture.bind();
		font_texture.resize(char_width*16, char_height*16);

		char *texture = new char[font_texture.width*font_texture.height];

		for (int i = 0; i < font_texture.width*font_texture.height; i++)
			texture[i] = 0;

		// current position in texture to write glyph
		int tex_base_x = 0;
		int tex_base_y = 0;

		for (int i = 0; i < 256; i++)
		{
			unsigned int glyph_index = FT_Get_Char_Index(ft_face, (unsigned char)(i)%256);

			// load glyph to slot (render it)
			ft_error = FT_Load_Glyph(ft_face, glyph_index, FT_LOAD_RENDER);

			FT_Bitmap bitmap = ft_slot->bitmap;

			for (std::size_t by = 0; by < bitmap.rows; by++)
			{
				int ty = tex_base_y + (char_height-1) - (by+char_origin_top-ft_slot->bitmap_top);

				for (std::size_t bx = 0; bx < bitmap.width; bx++)
				{
					int tx = tex_base_x + char_origin_left + ft_slot->bitmap_left + bx;

					texture[tx + ty*font_texture.width] = bitmap.buffer[bx + by*bitmap.width];
				}
			}

			advance_x[i] = ft_slot->advance.x/64 + i_char_spacing;

			tex_base_x += char_width;
			if (tex_base_x >= font_texture.width)
			{
				tex_base_x = 0;
				tex_base_y += char_height;
			}
		}

		font_texture.setData(texture);
		font_texture.unbind();
		delete [] texture;

		FT_Done_Face(ft_face);
		FT_Done_FreeType(ft_library);


		/*
		 * GL shader stuff
		 */
		program.use();
		glyph_size_uniform.set(GLSL::ivec2(char_width, char_height));
		program.disable();

		setPosition(0,font_size);

		valid = true;

		CGlErrorCheck();
	}

public:
	/**
	 * this function has to be called when the rendering viewport size is changed
	 */
	void viewportChanged(
			GLfloat i_width,
			GLfloat i_height
	)
	{
		program.use();
		inv_viewport_size_uniform.set(GLSL::vec2(1.0/(float)i_width, 1.0/(float)i_height));
		program.disable();
	}


	void setNewlineXPosition(int p_newline_x_position)
	{
		newline_x_position = p_newline_x_position;
	}


	void setColor(const GLSL::vec3 &color)
	{
		program.use();
		color_uniform.set(color);
		program.disable();
	}


	void setPosition(const std::array<float,2> &i_stylus_pos)
	{
		stylus_pos = i_stylus_pos;
	}


	void setPosition(GLfloat i_x, GLfloat i_y)
	{
		stylus_pos[0] = i_x;
		stylus_pos[1] = i_y;
	}


	void renderString(
		const char* i_text
	)
	{
		if (!valid)
		{
			std::cerr << "GLFreeType not yet valid, load font first" << std::endl;
			return;
		}

		int length = std::min<size_t>(strlen(i_text), 1024);

		// prepare vertex attributes
		int attrib_counter = 0;

		int *v = vertex_attrib;

		stylus_pos[0] -= char_origin_left;
		for (int i = 0; i < length; i++)
		{
			char c = i_text[i];
			if (c == '\n')
			{
				stylus_pos[0] = newline_x_position - char_origin_left;
				stylus_pos[1] -= font_size;
				continue;
			}

			int px = stylus_pos[0];
			int py = stylus_pos[1] - (char_height - char_origin_top);
			v[0] = px;	v[1] = py;	v[2] = c;	v[3] = 0;	v+=4;
			v[0] = px;	v[1] = py;	v[2] = c;	v[3] = 1;	v+=4;
			v[0] = px;	v[1] = py;	v[2] = c;	v[3] = 2;	v+=4;
			v[0] = px;	v[1] = py;	v[2] = c;	v[3] = 3;	v+=4;

			stylus_pos[0] += advance_x[(int)c];
			attrib_counter +=4;
		}
		stylus_pos[0] += char_origin_left;

		CGlStateDisable depth_test(GL_DEPTH_TEST);

		// enable blending
		glEnablei(GL_BLEND, 0);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		CGlErrorCheck();

		vao.bind();

			buffer.bind();
			buffer.subData(0, attrib_counter*4*sizeof(GLint), vertex_attrib);
			glVertexAttribIPointer(0, 4, GL_INT, 0, 0);
			glEnableVertexAttribArray(0);

				program.use();
				font_texture.bind();
					glDrawArrays(GL_TRIANGLE_STRIP, 0, attrib_counter);
				font_texture.unbind();
				program.disable();

		vao.unbind();

		glDisablei(GL_BLEND, 0);

		CGlErrorCheck();
	}

	int getTextLength(const char* text)
	{
		int length = 0;

		int chars = std::min<size_t>(::strlen(text), 1024);
		for (int i = 0; i < chars; i++)
			length += advance_x[(int)text[i]];

		return length;
	}
};


#endif /* CGLFREETYPE_HPP_ */
