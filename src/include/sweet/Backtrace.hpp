/*
 * Backtrace.hpp
 *
 *  Created on: Feb 20, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_BACKTRACE_HPP_
#define SRC_INCLUDE_SWEET_BACKTRACE_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/prctl.h>

namespace sweet
{

class Backtrace
{
public:
	/*
	 * We avoid using ErrorBase to avoid circular dependency
	 */
	std::string error;

	/*
	 * Source: https://stackoverflow.com/questions/4636456/how-to-get-a-stack-trace-for-c-using-gcc-with-line-number-information
	 *
	 * Get stacktrace
	 */
public:
	static
	void printGDBTrace()
	{
	    char pid_buf[30];
	    sprintf(pid_buf, "%d", getpid());

		// get name to program executable (which we want to debug)
		char program_exec_buf[512];
		std::size_t len = readlink("/proc/self/exe", program_exec_buf, 511);
		program_exec_buf[len] = '\0';

		prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY, 0, 0, 0);
		int child_pid = fork();
		if (!child_pid) {
			//dup2(2,1); // redirect output to stderr - edit: unnecessary?
			execl("/usr/bin/gdb", "gdb", "--batch", "-n", "-ex", "thread", "-ex", "bt", program_exec_buf, pid_buf, NULL);
			abort(); /* If gdb failed to start */
		} else {
			waitpid(child_pid,NULL,0);
		}
	}


	/*
	 * Use GDB to output detailed stacktrace which is returned in a string.
	 */
public:
	static
	std::string getGDBBacktrace()
	{
		/*
		 * Get PID of current process we want to debug
		 */
	    char pid_buf[30];
	    int pid = getpid();
	    sprintf(pid_buf, "%d", pid);

		/*
		 * Get path to program executable (which we want to debug)
		 */
		char program_exec_buf[1024];
		int len = readlink("/proc/self/exe", program_exec_buf, sizeof(program_exec_buf)-1);
		if (len == -1)
			SWEETError("Failed to read link");
		program_exec_buf[len] = '\0';

	    /*
	     * Create unidirectional pipeline
	     * pipefd[0]: read
	     * pipefd[1]: write
	     */
	    int pipefd_cout[2];
	    int result_pipe_cout = pipe(pipefd_cout);
	    if (result_pipe_cout == -1)
	    	SWEETError("Failed to create cout pipe");

	    int pipefd_cerr[2];
	    int result_pipe_cerr = pipe(pipefd_cerr);
	    if (result_pipe_cerr == -1)
	    	SWEETError("Failed to create cerr pipe");

	    /*
	     * Allow ptrace for forked processes
	     */
		prctl(PR_SET_PTRACER, PR_SET_PTRACER_ANY, 0, 0, 0);

	    /*
	     * Create fork
	     */
		int child_pid = fork();
		if (child_pid == -1)
	    	SWEETError("Failed to create fork");

		std::string pipe_buffer;
		if (child_pid == 0)
		{
			/*
			 * Redirect stdout and stderr to pipeline
			 */
		    dup2(pipefd_cout[1], STDOUT_FILENO);
		    dup2(pipefd_cerr[1], STDERR_FILENO);

		    /*
		     * Close reading ends of pipe
		     */
		    close(pipefd_cout[0]);
		    close(pipefd_cerr[0]);

			/*
			 * Run GDB
			 */
			execl(
					"/usr/bin/gdb",
					"gdb",
					"--batch",
					"-n",
					"-ex",
					"thread",
					"-ex",
					"bt",
					program_exec_buf,
					pid_buf,
					NULL
				);

			/*
			 * Close writing ends
			 */
		    close(pipefd_cout[1]);
		    close(pipefd_cerr[1]);

			abort();	// directly stop program without calling any atexit functions
		}
		else
		{
		    /*
		     * Close writing ends of pipe
		     */
		    close(pipefd_cout[1]);
		    close(pipefd_cerr[1]);

			/*
			 * Process std::cout
			 */
			{
				FILE* fcout = fdopen(pipefd_cout[0], "r");
				char buf[1024];

				while (true)
				{
					// automatic '\0' termination
					char *s = fgets(buf, sizeof(buf), fcout);

					if (s == nullptr)
						break;

					pipe_buffer += buf;
				}
			}

			/*
			 * Process std::cerr
			 */
			{
				FILE* fcout = fdopen(pipefd_cerr[0], "r");
				char buf[1024];

				while (true)
				{
					// automatic '\0' termination
					char *s = fgets(buf, sizeof(buf), fcout);

					if (s == nullptr)
						break;

					pipe_buffer += buf;
				}
			}

		    /*
		     * Close reading ends of pipe
		     */
		    close(pipefd_cout[0]);
		    close(pipefd_cerr[0]);

		}

		waitpid(child_pid, nullptr, 0);
		return pipe_buffer;
	}


	/*
	 * From https://stackoverflow.com/questions/3899870/print-call-stack-in-c-or-c/26529030
	 *
	 * This doensn't include any line numbers, etc. :-(
	 */
public:
	static
	std::string getBacktrace()
	{
		const int MAX_SIZE = 1024;
		void *array[1024];
		int size = backtrace(array, MAX_SIZE);
		char **strings = backtrace_symbols(array, size);

		std::string retval;
		for (int i = 0; i < size; i++)
		{
			retval += strings[i];
			retval += "\n";
		}

		free(strings);
		return retval;
	}
};

}

#endif /* SRC_INCLUDE_SWEET_BACKTRACE_HPP_ */
