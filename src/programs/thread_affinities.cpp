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
		std::map<int, int> &tid_to_worker_id
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

	pid_t tid = gettid();
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
	}
}

int main(int argc, char *argv[])
{
#if SWEET_MPI
	MPI_Init(&argc, &argv);
#endif

	int max_threads = omp_get_max_threads();

	std::map<int, int> tid_to_worker_id;

	#pragma omp parallel for
	for (int j = 0; j < max_threads; j++)
	{
		pid_t tid = gettid();

#pragma omp critical
		tid_to_worker_id[tid] = j;
	}


	if (true)
	{
		std::vector<std::ostringstream> ss;
		ss.resize(max_threads);

		std::cout << "#pragma omp parallel" << std::endl;

		std::atomic<int> counter(0);
		#pragma omp parallel
		{
			schedInfo(ss, tid_to_worker_id);

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
			schedInfo(ss, tid_to_worker_id);

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

		std::cout << "#pragma omp parallel for" << std::endl;
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
			schedInfo(ss, tid_to_worker_id);

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

		std::atomic<int> counter(0);

		std::cout << "#pragma omp parallel num_threads(max_threads/2)" << std::endl;
		std::cout << "	#pragma omp parallel num_threads(2)" << std::endl;
		hline_dash();

		#pragma omp parallel num_threads(max_threads/2)
		{
			std::atomic<int> counter2(0);

			#pragma omp parallel num_threads(2)
			{
				schedInfo(ss, tid_to_worker_id);

				counter2++;

				while (counter2 != 2);
			}

			counter++;

			while (counter != max_threads/2);
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
			std::atomic<int> counter2(0);

			#pragma omp parallel num_threads(2)
			#pragma omp for schedule(static,1)
			for (int k = 0; k < 2; k++)
			{
				schedInfo(ss, tid_to_worker_id);

				counter2++;

				while (counter2 != 2);
			}

			counter++;

			while (counter != max_threads/2);
		}

		for (int i = 0; i < max_threads; i++)
			std::cout << " + worker " << i << ": " << ss[i].str() << std::endl;

	}

#if SWEET_MPI
	MPI_Finalize();
#endif
	return(0);
}

