/*
 * FatalError.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_FATALERROR_HPP_
#define SRC_INCLUDE_SWEET_FATALERROR_HPP_


#include <cassert>
#include <string>
#include <iostream>
#include <signal.h>
#include <stdlib.h>

// For ptrace
#include <stdio.h>
#include <execinfo.h>


class FatalError
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
	// [[ noreturn ]] FatalError(const std::string &i_error)
	FatalError(const std::string &i_error)
	{
		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << "ERROR: " << i_error << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		print_trace();
		assert(false);
		raise(SIGABRT);
		exit(-1);
	}
};



class AssertFatalError
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
#ifdef NDEBUG
	AssertFatalError(bool i_assertion, const std::string &i_error)	{}
#elif !SWEET_DEBUG
	AssertFatalError(bool i_assertion, const std::string &i_error)	{}
#else
	AssertFatalError(bool i_assertion, const std::string &i_error)
	{
		if (i_assertion)
			return;

		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << "ASSERT ERROR: " << i_error << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		print_trace();
		assert(false);
		exit(-1);
	}
#endif

};


#endif /* SRC_INCLUDE_SWEET_FATALERROR_HPP_ */
