#include <omp.h>
#include <sweet/NUMABlockAlloc.hpp>
#include <vector>
#include <iostream>



int main(int argc, char *argv[])
{
	std::cout << "MAIN" << std::endl;
	NUMABlockAlloc::setup();

	int num_threads = omp_get_max_threads();

	////////////////////////////////////////////////////////////
	std::vector<double*> data_a_1024;
	data_a_1024.resize(num_threads);

	std::cout << "FOR START" << std::endl;
#pragma omp parallel for
	for (int i = 0; i < num_threads; i++)
		data_a_1024[i] = NUMABlockAlloc::alloc<double>(1024);			// alloc a
	std::cout << "FOR END" << std::endl;

	////////////////////////////////////////////////////////////
	std::vector<double*> data_b_1024;
	data_b_1024.resize(num_threads);

#pragma omp parallel for
	for (int i = 0; i < num_threads; i++)
		data_b_1024[i] = NUMABlockAlloc::alloc<double>(1024);			// alloc b

	////////////////////////////////////////////////////////////
	std::vector<double*> data_c_2024;
	data_c_2024.resize(num_threads*2);

#pragma omp parallel for
	for (int i = 0; i < num_threads*2; i++)
		data_c_2024[i] = NUMABlockAlloc::alloc<double>(2024);			// alloc c

	////////////////////////////////////////////////////////////

#pragma omp parallel for
	for (int i = 0; i < num_threads*2; i++)
		NUMABlockAlloc::free(data_c_2024[i], 2024);						// free c

	////////////////////////////////////////////////////////////

#pragma omp parallel for
	for (int i = 0; i < num_threads; i++)
		NUMABlockAlloc::free(data_b_1024[i], 1024);						// free b

	////////////////////////////////////////////////////////////

#pragma omp parallel for
	for (int i = 0; i < num_threads; i++)
		NUMABlockAlloc::free(data_a_1024[i], 1024);						// free a


	return 0;
}
