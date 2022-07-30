/*
 * FatalError.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SWEETERROR_HPP_
#define SRC_INCLUDE_SWEET_SWEETERROR_HPP_


#include <cassert>
#include <string>
#include <iostream>
#include <signal.h>
#include <stdlib.h>

// For ptrace
#include <stdio.h>
#include <execinfo.h>



/************************************************************
 * SWEETError in release and debug mode
 ************************************************************/

class SWEETError_
{
	// From https://stackoverflow.com/questions/3899870/print-call-stack-in-c-or-c/26529030
	/* Paste this on the file you want to debug. */
	void print_trace(void) {
		char **strings;
		size_t i, size;
		enum Constexpr { MAX_SIZE = 1024 };
		void *array[MAX_SIZE];
		size = backtrace(array, MAX_SIZE);
		strings = backtrace_symbols(array, size);
		for (i = 0; i < size; i++)
		printf("%s\n", strings[i]);
		puts("");
		free(strings);
	}


public:
//	[[ noreturn ]]
	SWEETError_(
			const std::string &i_error_type,
			const std::string &i_error_message,
			const char* i_filename,
			int i_line_no,
			const char* i_func,
			bool stop_after_error = true
		)
	{
		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;

		if (i_error_message != "")
			std::cerr << " " << i_error_type << ": " << i_error_message << std::endl;
		else
			std::cerr << " " << i_error_type << std::endl;

		std::cerr << "********************************************" << std::endl;
		std::cerr << " +     File: " << i_filename << std::endl;
		std::cerr << " + Line Nr.: " << i_line_no << std::endl;
		std::cerr << " + Function: " << i_func << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		if (stop_after_error)
		{
			print_trace();
			assert(false);	// Trigger every possible error we can trigger to suport debuggers
			raise(SIGABRT);
			exit(-1);
		}
	}
};


// Errors which should never happen
#define SWEETErrorInternal(msg)	SWEETError_("Internal SWEET ERROR", msg, __FILE__, __LINE__, __func__)

// Errors which should never happen
#define SWEETErrorTODO(msg)		SWEETError_("TODO", msg, __FILE__, __LINE__, __func__)

// Regular errors such as wrong time integration method, negative resolution, etc.
#define SWEETError(msg)			SWEETError_("ERROR", msg, __FILE__, __LINE__, __func__)

// Regular errors such as wrong time integration method, negative resolution, etc.
#define SWEETError_nostop(msg)			SWEETError_("ERROR", msg, __FILE__, __LINE__, __func__, false)




/************************************************************
 * SWEET Debug assertions active during debug mode
 ************************************************************/

#ifdef NDEBUG

	#define SWEETDebugAssert_msg(assertion, msg)
	#define SWEETDebugAssert(assertion)

#else

	class SWEETDebugAssert_
	{
	public:
		SWEETDebugAssert_(
				bool i_assertion,
				const std::string &i_error_message,
				const char* i_filename,
				int i_line_no,
				const char* i_func
		)
		{
			if (i_assertion)
				return;

			SWEETError_("ASSERTION ERROR", i_error_message, i_filename, i_line_no, i_func);
		}
	};

	#define SWEETDebugAssert_msg(assertion, msg)	SWEETDebugAssert_(assertion, msg, __FILE__, __LINE__, __func__)

	#define SWEETDebugAssert(assertion)			SWEETDebugAssert_(assertion, "", __FILE__, __LINE__, __func__)

#endif



/************************************************************
 * SWEET assertions which also work in release mode
 ************************************************************/

class SWEETAssert_
{
public:
	SWEETAssert_(
			bool i_assertion,
			const std::string &i_error_message,
			const char* i_filename,
			int i_line_no,
			const char* i_func
	)
	{
		if (i_assertion)
			return;

		SWEETError_("ASSERTION ERROR", i_error_message, i_filename, i_line_no, i_func);
	}
};

#define SWEETAssert(assertion, msg)	SWEETAssert_(assertion, msg, __FILE__, __LINE__, __func__)


#endif /* SRC_INCLUDE_SWEET_SWEETERROR_HPP_ */
