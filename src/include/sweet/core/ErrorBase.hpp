/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_ERRORBASE_HPP_
#define SRC_INCLUDE_SWEET_ERRORBASE_HPP_

#include <string>
#include <ostream>
#include <sweet/core/SWEETError.hpp>
#include <sweet/core/Backtrace.hpp>

/*
 * Do an error check.
 * If there's an error: forward error and return false
 */
#define ERROR_CHECK_WITH_RETURN_BOOLEAN(classWithError)	\
	{ if ((classWithError).error.exists()) return error.forwardWithPositiveReturn((classWithError).error); }

/*
 * Do an error check.
 * If there's an error: forward error and return with EXIT_FAILURE
 */
#define ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(classWithError)	\
	{ if ((classWithError).error.exists()) { (classWithError).error.print(); return EXIT_FAILURE; } }


/*
 * Do an error check.
 * If there's an error: forward error and return
 * If there's no error: continue
 */
#define ERROR_CHECK_WITH_RETURN(classWithError) \
	{ if ((classWithError).error.exists()) { error.forward((classWithError).error); return; } }


/*
 * If there's an error: forward error and return false
 * If there's no error: return true
 */
#define ERROR_FORWARD_WITH_RETURN_BOOLEAN(classWithError)	\
	{ return error.forwardWithPositiveReturn((classWithError).error); }

/*
 * This is useful for contructors
 * If there's an error: simply forward
 */
#define ERROR_FORWARD(classWithError)	\
	{ error.forward((classWithError).error); }


#if SWEET_DEBUG
	#define SWEET_ASSERT(x)	\
		if (!(x)) \
			{ sweet::Backtrace b; std::cerr << b.getGDBBacktrace() << std::endl; };

	#define SWEET_ASSERT_MSG(x, msg)	\
		if (!(x)) \
			{ std::cerr << msg << std::endl; sweet::Backtrace b; std::cerr << b.getGDBBacktrace() << std::endl; };
#else
	#define SWEET_ASSERT(x)
	#define SWEET_ASSERT_MSG(x, msg)
#endif


namespace sweet
{

class ErrorBase
{
	bool _hasError;
	std::string _errorMessage;


private:
	/*
	 * Check for environment variable SWEET_ERROR_WITH_STACKTRACE and if it exists,
	 * add stack trace to error information
	 */
	bool _withStacktrace()
	{
		static bool envLoaded = false;
		static bool envExists = false;

		if (!envLoaded)
		{
			char *val = std::getenv("SWEET_ERROR_WITH_STACKTRACE");
			envExists = (val != nullptr);
			envLoaded = true;
		}
		return envExists;
	}


public:
	ErrorBase()	:
		_hasError(false)
	{
	}

	/*
	 * Set an error
	 *
	 * \return always *false* to be able to use a one-liner error message
	 */
	bool set(const std::string &i_errorMessage)
	{
		if (_hasError)
		{
			_errorMessage += " | " + i_errorMessage;
		}
		else
		{
			_hasError = true;
			_errorMessage = i_errorMessage;
		}


		if (_withStacktrace())
		{
			std::string gdb = Backtrace::getGDBBacktrace();

			if (gdb != "")
			{
				_errorMessage += "\n";
				_errorMessage += "Stacktrace (from GDB):\n";
				_errorMessage += Backtrace::getGDBBacktrace();
			}
		}
		return false;
	}

	bool forward(ErrorBase &i_error)
	{
		if (!i_error._hasError)
			return false;

		_hasError = i_error._hasError;
		_errorMessage = i_error._errorMessage;

		i_error.reset();

		return true;
	}

	/*
	 * Forward errors
	 *
	 * \return **true** if there's no error
	 */
	bool forwardWithPositiveReturn(ErrorBase &i_error)
	{
		if (!i_error._hasError)
			return true;

		_hasError = i_error._hasError;
		_errorMessage = i_error._errorMessage;

		i_error.reset();

		return false;
	}

	bool exists() const
	{
		return _hasError;
	}

	void assertNoError() const
	{
		if (_hasError)
		{
			std::cerr << "Internal problem detected" << std::endl;
			SWEETError("Use backtrace to debug this problem!");
			exit(1);
		}
	}

	void reset()
	{
		_hasError = false;
		_errorMessage = "";
	}

	std::string get()
	{
		if (!_hasError)
		{
			std::cerr << "Error message requested, but has no error!" << std::endl;
			SWEETError("Use backtrace to debug this problem!");
			exit(1);
		}

		std::string tmp = _errorMessage;
		reset();
		return tmp;
	}

public:
	void print(std::ostream &io_os = std::cerr)
	{
		io_os << "ERROR: " << get() << std::endl;
	}

	~ErrorBase()
	{
		if (_hasError)
		{
			std::cerr << "ERROR was not processed!" << std::endl;
			std::cerr << "************************************************************" << std::endl;
			std::cerr << _errorMessage << std::endl;
			std::cerr << "************************************************************" << std::endl;
			exit(1);
		}
	}
};

}

#endif
