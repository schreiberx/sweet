/*
 * Get information about OMP parallelization
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <atomic>

#include <sched.h>
#include <omp.h>
#include <time.h>

#include <sys/syscall.h>
#include <unistd.h>

#if SWEET_MPI
	#include <mpi.h>
#endif



/*
 * Partly based on this code snipped:
 *
 * https://web.archive.org/web/20150209010321/http://docs.cray.com/books/S-2496-4101/html-S-2496-4101/cnlexamples.html
 */
/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
	char *ptr = str;
	int i, j, entry_made = 0;
	for (i = 0; i < CPU_SETSIZE; i++) {
		if (CPU_ISSET(i, mask)) {
			int run = 0;
			entry_made = 1;
			for (j = i + 1; j < CPU_SETSIZE; j++) {
				if (CPU_ISSET(j, mask)) run++;
				else break;
			}
			if (!run)
				sprintf(ptr, "%d,", i);
			else if (run == 1) {
				sprintf(ptr, "%d,%d,", i, i + 1);
				i++;
			} else {
				sprintf(ptr, "%d-%d,", i, i + run);
				i += run;
			}
			while (*ptr != 0) ptr++;
		}
	}
	ptr -= entry_made;
	*ptr = 0;
	return(str);
}



void hline_star()
{
	std::cout << "********************************************************************************" << std::endl;
}



void hline_dash()
{
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
}



void schedInfo(
		std::vector<std::ostringstream> &ss,
		std::map<int, int> &tid_to_worker_id,
		const std::string &i_str_iter
)
{
	char hostname_buf[1024];
	gethostname(hostname_buf, sizeof(hostname_buf));

	int mpi_rank;

	#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	#else
		mpi_rank = 0;
	#endif

	int thread_num = omp_get_thread_num();

#if 0
	// not supported on old systems
	pid_t tid = gettid();
#else
	pid_t tid = syscall(SYS_gettid);
#endif
	int worker_id = tid_to_worker_id[tid];

	cpu_set_t coremask;
	sched_getaffinity(0, sizeof(coremask), &coremask);

	char cores_buf[1024];
	cpuset_to_cstr(&coremask, cores_buf);

#pragma omp critical
	{
		ss[worker_id] << "host=" << hostname_buf;
		ss[worker_id] << ", ";
		ss[worker_id] << "rank=" << mpi_rank;
		ss[worker_id] << ", ";
		ss[worker_id] << "id=" << worker_id;
		ss[worker_id] << ", ";
		ss[worker_id] << "thread_num=" << thread_num;
		ss[worker_id] << ", ";
		ss[worker_id] << "cores=" << cores_buf;

		if (i_str_iter != "")
		{
			ss[worker_id] << ", ";
			ss[worker_id] << i_str_iter;
		}
	}
}

int main(int argc, char *argv[])
{
#if SWEET_MPI
	MPI_Init(&argc, &argv);
#endif

	int max_threads = omp_get_max_threads();

	std::cout << "max_threads: " << max_threads << std::endl;

	//for (int i = 0; i < 3; i++)
	{
		std::cout << std::endl;
		for (int j = 0; j < 10; j++)
			std::cout << "WARNING ";
		std::cout << std::endl;

		std::cout << "This is a test program for OpenMP to test nested parallelism!" << std::endl;
		std::cout << "If your OpenMP runtime doesn't support all features, it may deadlock!" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "You need to activate OpenMP nesting" << std::endl;
		std::cout << "	OMP_MAX_ACTIVE_LEVELS=2" << std::endl;
		std::cout << "" << std::endl;

		for (int j = 0; j < 10; j++)
			std::cout << "WARNING ";
		std::cout << std::endl;

		std::cout << std::endl;
	}

	if (omp_get_max_active_levels() == 0)
	{
		std::cerr << "ERROR: Nesting is not active" << std::endl;
		std::cerr << "You need to activate OpenMP nesting" << std::endl;
		std::cerr << "	OMP_NESTED=true" << std::endl;
		return EXIT_FAILURE;
	}

	if (max_threads % 2 != 0)
	{
		std::cout << "Max threads needs to be an even number" << std::endl;
		exit(1);
	}

	std::map<int, int> tid_to_worker_id;

	std::atomic<int> counter(0);

	#pragma omp parallel for schedule(static,1)
	for (int j = 0; j < max_threads; j++)
	{
		pid_t tid = gettid();

#pragma omp critical
		tid_to_worker_id[tid] = j;

		counter++;
		while (counter != max_threads);
	}


	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);

		std::cout << "#pragma omp parallel" << std::endl;

		std::atomic<int> counter(0);
		#pragma omp parallel
		{
			int i = omp_get_thread_num();
			std::ostringstream ss_iter;
			ss_iter << "(i=" << i << ")";
			schedInfo(ss, tid_to_worker_id, ss_iter.str());

			counter++;
			while (counter != max_threads);
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}

	hline_dash();

	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);

		std::cout << "#pragma omp parallel num_threads(max_threads/2)" << std::endl;
		hline_dash();

		std::atomic<int> counter(0);
		#pragma omp parallel num_threads(max_threads/2)
		{
			int i = omp_get_thread_num();

			std::ostringstream ss_iter;
			ss_iter << "(i=" << i << ")";
			schedInfo(ss, tid_to_worker_id, ss_iter.str());

			counter++;

			while (counter != max_threads/2);
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}
	hline_dash();

	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);

		std::cout << "#pragma omp parallel num_threads(max_threads/2)" << std::endl;
		std::cout << "#pragma omp for schedule(static,1)" << std::endl;
		hline_dash();

		std::atomic<int> counter(0);

		/*
		 * VERY IMPORTANT:
		 * Use num_threads only on "omp parallel!"
		 */
		#pragma omp parallel num_threads(max_threads/2)
		/*
		 * After setting num_threads, iterate over for loop as usual
		 */
		#pragma omp for schedule(static,1)
		for (int i = 0; i < max_threads/2; i++)
		{
			std::ostringstream ss_iter;
			ss_iter << "(i=" << i << ")";
			schedInfo(ss, tid_to_worker_id, ss_iter.str());

			counter++;
			while (counter != max_threads/2);
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}

	hline_dash();

	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);


		std::cout << "#pragma omp parallel num_threads(max_threads/2)" << std::endl;
		std::cout << "	#pragma omp parallel num_threads(2)" << std::endl;
		hline_dash();

		std::atomic<int> counter(0);
		#pragma omp parallel num_threads(max_threads/2)
		{

			int i = omp_get_thread_num();

			#pragma omp parallel num_threads(2)
			{
				int j = omp_get_thread_num();

				std::ostringstream ss_iter;
				ss_iter << "(i=" << i << ", j=" << j << ")";
				schedInfo(ss, tid_to_worker_id, ss_iter.str());

				counter++;
				while (counter != max_threads);
			}
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}

	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);

		std::cout << "#pragma omp parallel num_threads(max_threads/2)" << std::endl;
		std::cout << "#pragma omp for schedule(static,1)" << std::endl;
		std::cout << "	#pragma omp parallel num_threads(2)" << std::endl;
		std::cout << "	#pragma omp for schedule(static,1)" << std::endl;
		hline_dash();

		std::atomic<int> counter(0);

		/*
		 * VERY IMPORTANT:
		 * Use num_threads only on "omp parallel!"
		 */
		#pragma omp parallel num_threads(max_threads/2)
		#pragma omp for schedule(static,1)
		for (int j = 0; j < max_threads/2; j++)
		{
			#pragma omp parallel num_threads(2)
			#pragma omp for schedule(static,1)
			for (int k = 0; k < 2; k++)
			{
				std::ostringstream ss_iter;
				ss_iter << "(j=" << j << ", k=" << k << ")";
				schedInfo(ss, tid_to_worker_id, ss_iter.str());

				counter++;
				while (counter != max_threads);
			}
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}

#if SWEET_MPI
	MPI_Finalize();
#endif
	return(0);
}

