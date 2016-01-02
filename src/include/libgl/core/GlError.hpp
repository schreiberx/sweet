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

#ifndef C_GL_ERRORS_HPP
#define C_GL_ERRORS_HPP

#include "libgl/incgl3.h"
#include <iostream>
#include <stdio.h>

/**
 * \brief Handle & output OpenGL errors
 *
 * the errors are printed immediately to std::cerr without buffering!
 */
class CGlError
{
public:
	CGlError()
	{
	}

/*
#if CHECK_GL_ERROR_OFF
	inline static void check(const char *file, const char *function, int line)
	{
	}
#else
*/
	/**
	 * check for OpenGL error
	 *
	 * this function should not be called directly. Use CGlErrorCheck() for OpenGL error checks
	 */
	static void check(	const char *file,			///< source filename
						const char *function,		///< name of function
						int line,					///< line in source file
						const char *message = nullptr	///< message to append error message
	)
	{
		GLenum errNr = glGetError();
		if (errNr == GL_NO_ERROR)
			return;

		const char *errorMessage = 0;
		char raw_message[100];


		switch(errNr)
		{
			case GL_INVALID_ENUM:
				errorMessage = "INVALID_ENUM";
				break;

			case GL_INVALID_VALUE:
				errorMessage = "INVALID_VALUE";
				break;
			
			case GL_INVALID_OPERATION:
				errorMessage = "INVALID_OPERATION";
				break;

			case GL_INVALID_FRAMEBUFFER_OPERATION:
				errorMessage = "INVALID_FRAMEBUFFER_OPERATION";
				break;

			case GL_OUT_OF_MEMORY:
				errorMessage = "OUT_OF_MEMORY";
				break;

			default:
				sprintf(raw_message, "UNKNOWN %i", errNr);
				errorMessage = raw_message;
				break;
		}

		std::string string_message;
		if (message != nullptr)
			string_message = message;

//		std::cerr << "\033[0;40;31m" << "GL ERROR  in   " << file << "(" << line << ")" << ": " << function << "() -> " << errorMessage << " [" << string_message << "\033[0m" << std::endl;
		std::cerr << "GL ERROR  in   " << file << "(" << line << ")" << ": " << function << "() -> " << errorMessage << " [" << string_message << std::endl;

#if 1
		std::cerr << "[creating null pointer exception to debug program...]" << std::endl;
		int *null_exception = 0;
		*null_exception = 0;
#endif
	}
//#endif

	/**
	 * manually trigger an OpenGL error
	 */
	static void trigger(const char *file, const char *function, int line, const char *errorMessage)
	{
		std::cerr << "\033[0;40;31m" << "GL ERROR TRIGGERED in   " << file << "(" << line << ")" << ": " << function << "() -> " << errorMessage << "\033[0m" << std::endl;
	}
//#endif
};


#if C3S_DEBUG_MODE
	/**
	 * check for gl error and output the line and function where the error occured
	 */
	#define CGlErrorCheck()	CGlError::check(__FILE__, __FUNCTION__, __LINE__)

	/**
	 * check for gl error and output the line and function where the error occured
	 *
	 * append 'message' to the OpenGL error message
	 */
	#define CGlErrorCheckWithMessage(message)	CGlError::check(__FILE__, __FUNCTION__, __LINE__, message)

	/**
	 * check for gl error and output the line and function where the error occured
	 */
	#define CGlErrorTrigger(message)	CGlError::trigger(__FILE__, __FUNCTION__, __LINE__, message)

#else

	#define CGlErrorCheck()	{}
	#define CGlErrorCheckWithMessage(message)	{}
	#define CGlErrorTrigger(message)	{}

#endif


#endif
