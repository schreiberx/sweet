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


#ifndef C_GL_PROGRAM_HPP
#define C_GL_PROGRAM_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlUniform.hpp>
#include <libmath/CGlSlMath.hpp>
#include <libgl/tools/File.hpp>
#include <list>
#include <string>

#include "libgl/shaders/CDefaultShaderDir.hpp"


#include "libgl/core/GlShader.hpp"

/**
 * \brief GLSL program abstraction class
 */
class GlProgram
{
public:
	GLuint program;					///< OpenGL program ID
	bool error;

	std::list<CGlShader> shaders;	///< list to attached shaders

	std::string shaders_dir;	///< for debug purposes
	std::string prefix_string;	///< string to place before the source code


protected:
	/**
	 * set a string which is used as a prefix for every sourcefile
	 */
	void setSourcePrefix(const std::string &p_prefix_string)
	{
		prefix_string = p_prefix_string;
	}

	/**
	 * initialize and allocate OpenGL program handler
	 */
	void init()
	{
		error = false;
		program = glCreateProgram();
		CGlErrorCheck();
	}

public:
	GlProgram()
	{
		init();
	}

	/**
	 * attach shader to OpenGL program
	 */
	void attachShader(CGlShader &shader)
	{
		glAttachShader(program, shader());
		CGlErrorCheck();
	}

	/**
	 * detach and free shader if valid
	 */
	void detachAndFreeShaders()
	{
		while (!shaders.empty())
		{
			CGlShader &shader = shaders.front();
			glDetachShader(program, shader());
			shaders.pop_front();
		}
	}

	/**
	 * load shader from file and attach shader to program
	 */
	bool attachShader(const char *shader_file, GLuint type)
	{
		// insert new shader at begin of texture
		CGlShader tmp_shader;
		shaders.push_front(tmp_shader);

		CGlShader &new_shader = shaders.front();
		new_shader.init(type);
		if (!new_shader.loadSource(shader_file, prefix_string))
		{
			std::cerr << "Load error in shader_file " << shader_file << std::endl;
			std::cerr << new_shader.getInfoLog() << std::endl;
			error = true;
			// call deconstructor and remove
			shaders.pop_front();
			return false;
		}

		if (!new_shader.compile())
		{
			std::cerr << "Compile error in shader_file " << shader_file << std::endl;
			std::cerr << new_shader.getInfoLog() << std::endl;
			error = true;

			// call deconstructor and remove
			shaders.pop_front();
			return false;
		}

		glAttachShader(program, new_shader());

		return true;
	}

	/**
	 * attach vertex shader loaded from shader_file
	 */
	bool attachVertShader(const char *shader_file)
	{
		return attachShader(shader_file, GL_VERTEX_SHADER);
	}

	/**
	 * attach geometry shader loaded from shader_file
	 */
	bool attachGeomShader(const char *shader_file)
	{
		return attachShader(shader_file, GL_GEOMETRY_SHADER);
	}

	/**
	 * attach fragment shader loaded from shader_file
	 */
	bool attachFragShader(const char *shader_file)
	{
		return attachShader(shader_file, GL_FRAGMENT_SHADER);
	}

	/**
	 * initialize and attach vertex and fragment shaders from file
	 */
	bool initVertFragShaders(
				const char *vert_shader_file,
				const char *frag_shader_file
	)
	{
		std::string infoLog;

		detachAndFreeShaders();

		if (!attachVertShader(vert_shader_file))
			return false;

		if (!attachFragShader(frag_shader_file))
			return false;

		return true;
	}


	/**
	 * initialize and attach geometry, vertex and fragment shaders from file
	 */
	bool initGeomVertFragShaders(
				const char *vert_shader_file,
				const char *geom_shader_file,
				const char *frag_shader_file
	)
	{
		std::string infoLog;
		detachAndFreeShaders();


		if (!attachVertShader(vert_shader_file))
			return false;

		if (!attachGeomShader(geom_shader_file))
			return false;

		if (!attachFragShader(frag_shader_file))
			return false;

		return true;
	}

	/**
	 * initialize and attach vertex and fragment shaders from directory assuming that
	 * vertex shader is given in [directory]/vertex_shader.glsl and
	 * fragment shader is given in [directory]/fragment_shader.glsl
	 */
	bool initVertFragShadersFromDirectory(const char *directory)
	{
		shaders_dir = directory;

		std::string dir;

		dir = SHADER_GLSL_DEFAULT_DIR;
		dir += directory;
		dir += "/";

		std::string vert_shader = dir;
		vert_shader += "vertex_shader.glsl";

		std::string fragment_shader = dir;
		fragment_shader += "fragment_shader.glsl";

		return initVertFragShaders(vert_shader.c_str(), fragment_shader.c_str());
	}


	/**
	 * initialize and attach geometry, vertex and fragment shaders from directory assuming that
	 * vertex shader is given in [directory]/vertex_shader.glsl and
	 * geometry shader is given in [directory]/geometry_shader.glsl and
	 * fragment shader is given in [directory]/fragment_shader.glsl
	 */
	bool initGeomVertFragShadersFromDirectory(const char *directory)
	{
		shaders_dir = directory;

		std::string dir;
		dir = "shaders/";
		dir += directory;
		dir += "/";

		std::string geom_shader = dir;
		geom_shader += "geometry_shader.glsl";

		std::string vert_shader = dir;
		vert_shader += "vertex_shader.glsl";

		std::string fragment_shader = dir;
		fragment_shader += "fragment_shader.glsl";

		return initGeomVertFragShaders(vert_shader.c_str(), geom_shader.c_str(), fragment_shader.c_str());
	}




	/**
	 * return true, if information log is valid and store information log to parameter
	 */
	bool getInfoLog(std::string &infoLog)
	{
		GLint length;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);

		if (length == 0)
		{
			// no info log available
			infoLog = "Info log has zero length";
			return false;
		}

		GLchar *info_log_buf = new GLchar[length];

		// returned string is already zero terminated
		glGetProgramInfoLog(program, length, nullptr, info_log_buf);

		infoLog = info_log_buf;
		return true;
	}

	/**
	 * return string for info log
	 */
	std::string getInfoLog()
	{

		GLint length = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);

		if (length == 0)
			return "Info log has zero length";

		GLchar *info_log_buf = new GLchar[length];

		// returned string is already zero terminated
		glGetProgramInfoLog(program, length, nullptr, info_log_buf);

		return info_log_buf;
	}

	/**
	 * link compiled shaders
	 */
	bool link()
	{
		glLinkProgram(program);
		CGlErrorCheck();

		GLint status;
		glGetProgramiv(program, GL_LINK_STATUS, &status);

		if (status == GL_TRUE)
			return true;

		std::cerr << getInfoLog() << std::endl;
		error = true;

		for (CGlShader &s : shaders)
		{
			std::cerr << " Shader filename: " << s.filename << std::endl;
			error = true;
		}
		return false;
	}

	/**
	 * bind attribute (input) location
	 */
	inline void bindAttribLocation(GLuint index, const GLchar *name)
	{
		glBindAttribLocation(program, index, name);
		CGlErrorCheck();
	}

	/**
	 * bind fragment data (output) location
	 */
	void bindFragDataLocation(GLuint color_number, const GLchar *name)
	{
		glBindFragDataLocation(program, color_number, name);
		CGlErrorCheck();
	}

	/**
	 * set uniform to the value
	 */
	inline void setUniform1i(const GLchar *name, GLint value)
	{
		GlUniform uniform;
		setupUniform(uniform, name);
		uniform.set1i(value);

		CGlErrorCheckWithMessage(name);
	}

	/**
	 * set uniform to the 3 components of value
	 */
	inline void setUniform3fv(const GLchar *name, GLfloat *value)
	{
		GlUniform uniform;
		setupUniform(uniform, name);
		uniform.set3fv(value);

		CGlErrorCheckWithMessage(name);
	}

	/**
	 * return the uniform location for a given uniform
	 */
	inline GLint getUniformLocation(const GLchar *p_name)
	{
		return glGetUniformLocation(program, p_name);
	}


	/**
	 * initialize uniform from program
	 */
	inline void setupUniform(
			GlUniform &p_uniform,	///< uniform to setup
			const GLchar *p_name	///< name of uniform variable in program
	)
	{
		p_uniform.name = p_name;
		p_uniform.location = glGetUniformLocation(program, p_name);

		if (p_uniform.location == -1)
		{
			/**
			 * just output some useful information
			 *
			 * if location is set to -1, setting a value for such a uniform variable
			 * produces no error
			 */
#if DEBUG
			std::cout << "!!! info: uniform location \"" << p_name << "\" not found";
			if (shaders_dir != "")
				std::cout << " for shaders \"" << shaders_dir << "\"";
			std::cout << "!!!" << std::endl;
#endif
		}
		CGlErrorCheck();
	}


	/**
	 * use the program
	 */
	inline 	void use()
	{
		glUseProgram(program);
		CGlErrorCheck();
	}

	/**
	 * disable usage of program
	 */
	static inline 	void disable()
	{
#if DEBUG==1
		glUseProgram(0);
		CGlErrorCheck();
#endif
	}

	/**
	 * return state of program (true, if program can be used)
	 */
	inline bool validate()
	{
		glValidateProgram(program);

		GLint status;
		glGetProgramiv(1, GL_VALIDATE_STATUS, &status);

		return status == GL_TRUE;
	}

	/**
	 * free and delete OpenGL program
	 */
	inline void free()
	{
		if (program == 0)
			return;
		glDeleteProgram(program);
		program = 0;
	}

	inline ~GlProgram()
	{
		free();
		CGlErrorCheck();
	}
};

/**
 * \brief	convenient function to activate usage of programs within {} blocks
 */
class GlProgramUse
{
public:
	/**
	 * activate usage of program in current program block {}
	 */
	inline GlProgramUse(GlProgram &p_program)
	{
		p_program.use();
	}

	/**
	 * disable program
	 */
	inline ~GlProgramUse()
	{
		GlProgram::disable();
	}
};

#endif
