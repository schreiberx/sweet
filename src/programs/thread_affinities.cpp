/*
 * This script originates from Brian Dobbins (NCAR) who
 * might have got it from another person.
 *
 * It is not directly related to SWEET, but helps to identify affinities
 * if running it on supercomputers.
 */

/*
 * http://docs.cray.com/books/S-2496-4101/html-S-2496-4101/cnlexamples.html
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#if SWEET_MPI
	#include <mpi.h>
#endif
#include <omp.h>
#include <time.h>

#include <stdlib.h>     /* srand, rand */


#define SEED 35791246
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
int main(int argc, char *argv[])
{
	int rank, thread;
	cpu_set_t coremask;
	int niter = 100000;            //number of iterations per FOR loop
	double x,y;                     //x,y value for the random coordinate
	int i;                          //loop counter
	int count=0;                //Count holds all the number of how many good coordinates
	double z;                       //Used to check if x^2+y^2<=1
	char clbuf[7 * CPU_SETSIZE], hnbuf[64];
	#if SWEET_MPI
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
		rank = 0;
	#endif
	srand(SEED);
	memset(clbuf, 0, sizeof(clbuf));
	memset(hnbuf, 0, sizeof(hnbuf));
	(void)gethostname(hnbuf, sizeof(hnbuf));

	#pragma omp parallel firstprivate(x, y, z, i) private(thread, coremask, clbuf)
	{
		thread = omp_get_thread_num();
		clock_t t;
		t= clock();
		/* Borrowed from https://www.olcf.ornl.gov/tutorials/monte-carlo-pi/ By Jake Wynn */
		//Let's do some work generating random nubmers to see if core affinity impacts execution time.
		if((rank==0)||(rank==1)) // put this labor on ranks 0 and 1.
		{
			srandom((int)time(NULL) ^ omp_get_thread_num());    //Give random() a seed value
			for (i=0; i<niter; ++i)              //main loop
			{
				x = (double)random()/1989.98;      //gets a random x coordinate
				y = (double)random()/1974.9171;      //gets a random y coordinate
				z = ((x*x)+(y*y));          //Checks to see if number is inside unit circle
				++count;            //if it is, consider it a valid random point
			}
		}
		(void)sched_getaffinity(0, sizeof(coremask), &coremask);
		cpuset_to_cstr(&coremask, clbuf);

		#pragma omp barrier

		t = clock() - t;
		printf("Rank %d, thread %d, on %s. core = %s,(%f seconds).\n",
				rank, thread, hnbuf, clbuf,((float)t)/CLOCKS_PER_SEC);
	}
#if SWEET_MPI
	MPI_Finalize();
#endif
	return(0);
}

