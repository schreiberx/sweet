/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */


//#include <sched.h>
#include <numa.h>
#include <omp.h>
#include <iostream>


#include <sched.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>


class Affinity
{
public:
	void print()
	{

		cpu_set_t *cpusetp;
		size_t size;
		int num_cpus, cpu;
/*
		if (argc < 2) {
			fprintf(stderr, "Usage: %s <num-cpus>\n", argv[0]);
			exit(EXIT_FAILURE);
		}
		num_cpus = atoi(argv[1]);
*/
		num_cpus = 4;
		cpusetp = CPU_ALLOC(num_cpus);
		if (cpusetp == NULL) {
			perror("CPU_ALLOC");
			exit(EXIT_FAILURE);
		}
		size = CPU_ALLOC_SIZE(num_cpus);
		CPU_ZERO_S(size, cpusetp);
		for (cpu = 0; cpu < num_cpus; cpu += 2)
			CPU_SET_S(cpu, size, cpusetp);
		printf("CPU_COUNT() of set:    %d\n", CPU_COUNT_S(size, cpusetp));
		CPU_FREE(cpusetp);

		int num_configured_nodes = numa_num_configured_nodes();
		std::cout << "num_configured_nodes: " << num_configured_nodes << std::endl;

		int num_configured_cpus = numa_num_configured_cpus();
		std::cout << "num_configured_cpus: " << num_configured_cpus << std::endl;

/*
		struct cpuCPU_ALLOC()
		struct bitmask mask;
		int retval = sched_getaffinity(0, &mask);
		if (retval != 0)
		{
			std::cout << "ERROR: numa_sched_getaffinity" << std::endl;
			exit(1);
		}
		std::cout << "sched_getaffinity: ";
		unsigned int mask_size = numa_bitmask_nbytes(&mask);
		std::cout << "mask_size: " << mask_size << std::endl;
		for (int i = 0; i < mask_size; i++)
			std::cout << (numa_bitmask_isbitset(&mask, i) == 1 ? '1' : '0');
		std::cout << std::endl;
*/
		//int numa_sched_getaffinity(pid_t pid, struct bitmask *mask);
		
/*
		std::vector<cpu_set_t> mask;
		mask.resize()
		size_t cpu_set_t mask


int sched_getaffinity(pid_t pid, size_t cpusetsize,
                      cpu_set_t *mask);
		unsigned max_threads = omp_get_max_threads();
		std::cout << "max_threads: " << max_threads << std::endl;

		for (int j = 0; j < max_threads; j++)
		{
			#pragma omp parallel for schedule(static,1)
			for (int i = 0; i < max_threads; i++)
			{
				getcpu(&cpu, &node, nullptr);
				std::cout << "Thread " << i << " runs on " << cpu << " which is on node " << node << std::endl;
			}
		}
*/
	}
};


int main(int argc, char *argv[])
{
	Affinity a;
	a.print();

	return 0;
}

