/*
 * CGlShader.hpp
 *
 *  Created on: Sep 6, 2011
 *      Author: schreibm
 */

#ifndef CGLSHADER_HPP_
#define CGLSHADER_HPP_

/**
 * \brief	load & compile GLSL vertex, geometry and fragment shaders
 */
class CGlShader
{
	GLuint shader;


public:
	std::string filename;


	/**
	 * initialize with given type
	 */
	void init(GLenum type)
	{
		freeIfValid();

		shader = glCreateShader(type);
		CGlErrorCheck();
	}

	/**
	 * constructor with type
	 */
	CGlShader(GLenum type)
	{
		init(type);
	}



	/**
	 * initialize with existing shader
	 */
	CGlShader(const CGlShader &p_shader)
	{
		if (p_shader() != 0)
			std::cerr << "Shaders may not be initialized for copy constructor!!!" << std::endl;
		shader = 0;
	}



	/**
	 * default constructor
	 */
	CGlShader()
	{
		shader = 0;
	}



	/**
	 * access OpenGL shader id
	 */
	inline GLuint operator()()	const
	{
		return shader;
	}



	/**
	 * load source
	 */
	bool loadSource(	const std::string &i_filename,		///< filename with source
						const std::string &prefix_string	///< prefix to place before source (e. g. definitions)
	)
	{
		filename = i_filename;

		std::string content;
		std::string errorLog;

		FileContents file;
		if (!file.getTextFileContent(i_filename, content, errorLog))
		{
			return false;
		}

		content = prefix_string + content;

		const GLchar *content_str = content.c_str();
		glShaderSource(shader, 1, &content_str, nullptr);

		CGlErrorCheck();
		return true;
	}



	/**
	 * load source from filename
	 */
	bool loadSource(
			const std::string &i_filename
	)
	{
		filename = i_filename;

		std::string content;
		std::string errorLog;

		FileContents file;
		if (!file.getTextFileContent(i_filename, content, errorLog))
		{
			return false;
		}

		const GLchar *content_str = content.c_str();
		glShaderSource(shader, 1, &content_str, nullptr);

		CGlErrorCheck();
		return true;
	}



	/**
	 * compile the shader
	 */
	bool compile()
	{
		glCompileShader(shader);
		CGlErrorCheck();

		GLint status;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &status);

		if (status == GL_TRUE)
			return true;

		std::cerr << getInfoLog();
		return false;
	}



	/**
	 * return the information log for the compilation
	 */
	bool getInfoLog(std::string &infoLog)
	{
		GLint length;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);

		if (length == 0)
		{
			// no info log available
			infoLog = "";
			return false;
		}

		GLchar *info_log_buf = new GLchar[length];

		// returned string is already zero terminated
		glGetShaderInfoLog(shader, length, nullptr, info_log_buf);

		infoLog = info_log_buf;
		return true;
	}



	/**
	 * return the information log for the compilation
	 */
	std::string getInfoLog()
	{
		GLint length = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);

		if (length == 0)
			return "Info log has zero length";

		GLchar *info_log_buf = new GLchar[length];

		// returned string is already zero terminated
		glGetShaderInfoLog(shader, length, nullptr, info_log_buf);

		std::string ret_string = info_log_buf;
		return ret_string;
	}



	/**
	 * free the shader
	 */
	void free()
	{
		if (shader != 0)
		{
			glDeleteShader(shader);
			CGlErrorCheck();
			shader = 0;
		}
		else
		{
			std::cerr << "Warning: ~CGlShader: shader was not initialized!" << std::endl;
		}
	}

	/**
	 * return true, if the shader was loaded or created
	 */
	bool loaded()
	{
		return shader != 0;
	}

	/**
	 * free and delete shader data
	 */
	void freeIfValid()
	{
		if (shader != 0)
		{
			glDeleteShader(shader);
			CGlErrorCheck();
			shader = 0;
		}
	}

	~CGlShader()
	{
		freeIfValid();
	}
};


#endif /* CGLSHADER_HPP_ */
