/*
 * CGlUniform.hpp
 *
 *  Created on: Sep 6, 2011
 *      Author: schreibm
 */

#ifndef CGLUNIFORM_HPP_
#define CGLUNIFORM_HPP_

#include <array>

//#include "libgl/core/CGlProgram.hpp"

/**
 * \brief	handle uniforms for GLSL programs
 */
class CGlUniform
{
public:
	GLint location;		///< location id for uniform
	std::string name;	///< name of uniform for debugging purposes

	/**
	 * initialize uniform with OpenGL Uniform location Id
	 */
	inline CGlUniform(GLint p_uniformLocation)
	{
		location = p_uniformLocation;
	}

	/**
	 * initialize empty uniform
	 */
	inline CGlUniform()
	{
		location = -1;
	}


	/**
	 * copy operator to copy uniform ids
	 */
	inline void operator=(GLint p_uniformLocation)
	{
		location = p_uniformLocation;
	}

	/**
	 * access opengl uniform id with operator()
	 */
	inline GLint operator()()
	{
		return location;
	}

	/**
	 * set mat4 uniform
	 */
	inline void set(const GLSL::mat4 &matrix)
	{
		glUniformMatrix4fv(location, 1, GL_TRUE, GLSL::value_ptr(matrix));
		CGlErrorCheck();
	}

	/**
	 * set mat3 uniform
	 */
	inline void set(const GLSL::mat3 &matrix)
	{
		glUniformMatrix3fv(location, 1, GL_TRUE, GLSL::value_ptr(matrix));
		CGlErrorCheck();
	}

	/**
	 * set vec2 uniform
	 */
	inline void set(const GLSL::vec2 &vector)
	{
		glUniform2fv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}

	/**
	 * set vec3 uniform
	 */
	inline void set(const GLSL::vec3 &vector)
	{
		glUniform3fv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}

	/**
	 * set vec3 uniform
	 */
	inline void set(const std::array<float,3> &i_vector)
	{
		glUniform3fv(location, 1, i_vector.data());
		CGlErrorCheck();
	}
	/**
	 * set vec2 uniform
	 */
	inline void set(const std::array<float,2> &i_vector)
	{
		glUniform2fv(location, 1, i_vector.data());
		CGlErrorCheck();
	}

	/**
	 * set vec4 uniform
	 */
	inline void set(const GLSL::vec4 &vector)
	{
		glUniform4fv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}


	/**
	 * set ivec2 uniform
	 */
	inline void set(const GLSL::ivec2 &vector)
	{
		glUniform2iv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}

	/**
	 * set ivec3 uniform
	 */
	inline void set(const GLSL::ivec3 &vector)
	{
		glUniform3iv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}

	/**
	 * set ivec4 uniform
	 */
	inline void set(const GLSL::ivec4 &vector)
	{
		glUniform4iv(location, 1, GLSL::value_ptr(vector));
		CGlErrorCheck();
	}

	/**
	 * set boolean uniform
	 */
	inline void set1b(bool value)
	{
		glUniform1i(location, value);
		CGlErrorCheck();
	}

	/**
	 * set GLfloat uniform
	 */
	inline void set1f(const GLfloat value)
	{
		glUniform1f(location, value);
		CGlErrorCheck();
	}

	/**
	 * set GLint uniform
	 */
	inline void set1i(const GLint value)
	{
		glUniform1i(location, value);
		CGlErrorCheck();
	}

	/**
	 * set GLfloat[2] uniform
	 */
	inline void set2fv(const GLfloat *value)
	{
		glUniform2fv(location, 1, value);
		CGlErrorCheck();
	}

	/**
	 * set GLint[2] uniform
	 */
	inline void set2iv(const GLint *value)
	{
		glUniform2iv(location, 1, value);
		CGlErrorCheck();
	}

	/**
	 * set GLfloat[3] uniform
	 */
	inline void set3fv(const GLfloat *value)
	{
		glUniform3fv(location, 1, value);
		CGlErrorCheck();
	}

	/**
	 * set GLint[3] uniform
	 */
	inline void set3iv(const GLint *value)
	{
		glUniform3iv(location, 1, value);
		CGlErrorCheck();
	}
};



#endif /* CGLUNIFORM_HPP_ */
