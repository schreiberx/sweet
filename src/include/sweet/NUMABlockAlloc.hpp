/*
 * NUMABlockAlloc.hpp
 *
 *  Created on: 14 Sep 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_NUMABLOCKALLOC_HPP_
#define SRC_INCLUDE_SWEET_NUMABLOCKALLOC_HPP_

/**
 * define granularity of allocation
 *
 * 0: default allocator
 *
 * 1: NUMA domains
 *    -> Allocate one memory block chain per numa domain
 *    -> requires additional synchronization (critical regions)
 *
 * 2: Threads
 *    -> Allocate one memory block chain per thread!
 *    -> requires no synchronization
 *
 * 3: Non-NUMA:
 *    -> Allocate one memory block chain
 *    -> requires synchronization
 *    -> This is useful for XeonPhi KNC
 *
 */
#if !(NUMA_BLOCK_ALLOCATOR_TYPE >= 0 && NUMA_BLOCK_ALLOCATOR_TYPE <= 3)
#	error	"Please specify allocator type via NUMA_BLOCK_ALLOCATOR_TYPE"
#endif

#if NUMA_BLOCK_ALLOCATOR_TYPE == 1 || NUMA_BLOCK_ALLOCATOR_TYPE == 2
#	include <numa.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
#	include <omp.h>
#endif

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>





/**
 * This class implements a memory manager which caches the allocation of large memory blocks.
 *
 * The idea is to avoid freeing blocks directly.
 */
class NUMABlockAlloc
{
	/**
	 * Number of allocation domains
	 *
	 * Either this is the number of NUMA domains or the number of threads,
	 * depending on NUMA_BLOCK_ALLOCATOR_TYPE
	 */
	int num_alloc_domains = 1;

	/**
	 * verbosity
	 */
	int verbosity = 1;

	/**
	 * List of memory blocks of same size
	 */
private:
	class MemBlocksSameSize
	{
	public:
		/**
		 * size of memory blocks
		 */
		std::size_t block_size;

		/**
		 * Array of memory blocks
		 */
		std::vector<void*> free_blocks;
	};


	/**
	 * List of varying memory blocks of same size
	 */
	class DomainMemBlocks
	{
	public:
		/**
		 * Array of memory blocks
		 */
		std::vector<MemBlocksSameSize> block_groups;
	};


	/**
	 * List over all domains which are a
	 * -> list over all Block sizes which are a
	 * ---> list over all free blocks with that size
	 */
	std::vector<DomainMemBlocks> domain_block_groups;


	/**
	 * Hardware page size to use for memory alignment.
	 * (avoid overlapping pages for different NUMA nodes)
	 */
//	long page_size;

	/**
	 * setup already executed?
	 */
	bool setup_done;


private:
	inline
	static
	int& getThreadLocalDomainIdRef()
	{
#if 0
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		static int domain_id = omp_get_thread_num();
#else
		static int domain_id = 0;
#endif

		return domain_id;

#else

		/**
		 * Domain node for current thread
		 */
		// WARNING: This class has to be compiled always with
		static thread_local int domain_id;
		return domain_id;
#endif
	}


public:
	NUMABlockAlloc()	:
		setup_done(false)
	{
		p_setup();
	}


public:
	static
	void setup()
	{
		// access singleton to call constructor
		getSingletonRef();
	}

private:
	void p_setup()
	{
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		if (omp_in_parallel())
		{
			std::cerr << "ERROR: NUMAMemManager may not be initialized within parallel region!" << std::endl;
			std::cerr << "       Call NUMAMemManager::setup() at program start" << std::endl;
			exit(1);
		}
#endif

		if (setup_done)
			return;

		const char* env_verbosity = getenv("NUMA_BLOCK_ALLOC_VERBOSITY");
		if (env_verbosity == nullptr)
			verbosity = 0;
		else
			verbosity = atoi(env_verbosity);


#if  NUMA_BLOCK_ALLOCATOR_TYPE == 0

		if (verbosity > 0)
			std::cout << "NUMA block alloc: Using default system's allocator" << std::endl;

		num_alloc_domains = 1;
		getThreadLocalDomainIdRef() = 0;

#elif  NUMA_BLOCK_ALLOCATOR_TYPE == 1

		if (verbosity > 0)
			std::cout << "NUMA block alloc: Using NUMA node granularity" << std::endl;

		/*
		 * NUMA granularity
		 */
		num_alloc_domains = numa_num_configured_nodes();
		if (verbosity > 0)
			std::cout << "num_alloc_domains: " << num_alloc_domains << std::endl;

		// set NUMA id in case that master thread has a different id than the first thread
		int cpuid = sched_getcpu();
		getThreadLocalDomainIdRef() = numa_node_of_cpu(cpuid);


#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
#pragma omp parallel
		{
			int cpuid = sched_getcpu();
			getThreadLocalDomainIdRef() = numa_node_of_cpu(cpuid);
		}
#else
		getThreadLocalDomainIdRef() = 0;
#endif

#elif NUMA_BLOCK_ALLOCATOR_TYPE == 2

		if (verbosity > 0)
			std::cout << "NUMA block alloc: Using allocator based on thread granularity" << std::endl;

		/*
		 * Thread granularity, use this also per default
		 */
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		num_alloc_domains = omp_get_max_threads();
#else
		num_alloc_domains = 1;
#endif

		if (verbosity > 0)
			std::cout << "num_alloc_domains: " << num_alloc_domains << std::endl;

		// set NUMA id in case that master thread has a different id than the first thread
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		getThreadLocalDomainIdRef() = omp_get_thread_num();

#pragma omp parallel
		getThreadLocalDomainIdRef() = omp_get_thread_num();

#else
		getThreadLocalDomainIdRef() = 0;
#endif


#elif  NUMA_BLOCK_ALLOCATOR_TYPE == 3


		if (verbosity > 0)
			std::cout << "NUMA block alloc: Using non-numa single memory block chain" << std::endl;

		num_alloc_domains = 1;
		getThreadLocalDomainIdRef() = 0;

#else

#	error "Invalid NUMA_BLOCK_ALLOCATOR_TYPE"

#endif

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
		if (verbosity > 0)
		{
			#pragma omp parallel
			{
				#pragma omp critical
				{
					std::cout << "	thread id " << omp_get_thread_num() << " is assigned to memory allocator domain " << getThreadLocalDomainIdRef() << std::endl;
				}
			}
		}
#endif

		domain_block_groups.resize(num_alloc_domains);

#if 0
		// TODO: care about first-touch policy
		for (auto& n : domain_block_groups)
		{
			std::size_t S = num_alloc_domains*10;

			// preallocate S different size of blocks which should be sufficient
			n.block_groups.reserve(S);
		}
#endif

		setup_done = true;
	}



private:
	~NUMABlockAlloc()
	{
		if (verbosity > 1)
			std::cout << "NUMABlockAlloc EXIT" << std::endl;

		for (auto& n : domain_block_groups)
		{
			for (auto& g : n.block_groups)
			{
//				std::cout << "cleaning up " << g.free_blocks.size() << " blocks of size " << g.block_size << std::endl;

				for (auto& b : g.free_blocks)
				{
#if NUMA_BLOCK_ALLOCATOR_TYPE == 0 || NUMA_BLOCK_ALLOCATOR_TYPE == 3
					::free(b);
#else
					numa_free(b, g.block_size);
#endif
				}
			}
		}
	}


public:
	static
	inline
	NUMABlockAlloc& getSingletonRef()
	{
		/*
		 * Compilation details:
		 * this is compiled to a code such as
		 * '
		 * 		cmp memManager, 0
		 * 		jne return
		 * 		# memManager = NUMABlockAlloc()
		 * 		...
		 * 	return:
		 * 		ret ...
		 * '
		 */
		static NUMABlockAlloc memManager;
		return memManager;
	}



	/**
	 * return a list of blocks with the same size
	 *
	 * If the list does not exist, insert it
	 */
	static
	std::vector<void*>& getBlocksSameSize(
			std::size_t i_size				///< size of blocks
	)
	{
		NUMABlockAlloc &n = NUMABlockAlloc::getSingletonRef();

		assert(n.getThreadLocalDomainIdRef() < (int)n.domain_block_groups.size());

		std::vector<MemBlocksSameSize>& block_groups = n.domain_block_groups[n.getThreadLocalDomainIdRef()].block_groups;

		// iterate over blocks available for this NUMA domain
		for (auto& block_group : block_groups)
		{
			// check for matching block size
			if (block_group.block_size == i_size)
				return block_group.free_blocks;
		}

		// add new vector of blocks at the end with given size
		block_groups.emplace_back();
		block_groups.back().block_size = i_size;

		return block_groups.back().free_blocks;
	}



public:
	template <typename T=void>
	static
	inline
	T *alloc(
			std::size_t i_size		///< size of block
	)
	{
		T *data = nullptr;

#if NUMA_BLOCK_ALLOCATOR_TYPE == 1 || NUMA_BLOCK_ALLOCATOR_TYPE == 2

#	if NUMA_BLOCK_ALLOCATOR_TYPE == 1
		// dummy call here to initialize this class as part of the singleton out of critical region
		getSingletonRef();

#	if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
#		pragma omp critical
#	endif

#	endif
		{
			std::vector<void*>& block_list = getBlocksSameSize(i_size);

			if (block_list.size() > 0)
			{
				data = (T*)block_list.back();
				block_list.pop_back();
			}
		}

		if (data != nullptr)
			return data;

		return (T*)numa_alloc(i_size);

#elif NUMA_BLOCK_ALLOCATOR_TYPE == 3

#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
#	pragma omp critical
#endif
		{
			std::vector<void*>& block_list = getBlocksSameSize(i_size);

			if (block_list.size() > 0)
			{
				data = (T*)block_list.back();
				block_list.pop_back();
			}
		}

		if (data != nullptr)
			return data;

		int retval = posix_memalign((void**)&data, 4096, i_size);
		if (retval != 0)
		{
			std::cerr << "Unable to allocate memory" << std::endl;
			assert(false);
			exit(-1);
		}

		return data;

#else

		// allocate a new element to the list of blocks given in block_list

		// posix_memalign is thread safe
		// http://www.qnx.com/developers/docs/6.3.0SP3/neutrino/lib_ref/p/posix_memalign.html
		int retval = posix_memalign((void**)&data, 4096, i_size);
		if (retval != 0)
		{
			std::cerr << "Unable to allocate memory" << std::endl;
			assert(false);
			exit(-1);
		}

		return data;

#endif
	}



public:
	inline
	static
	void free(
			void *i_data,
			std::size_t i_size
	)
	{
		if (i_data == nullptr)
			return;

#if NUMA_BLOCK_ALLOCATOR_TYPE == 0

		::free(i_data);

#else

	#if NUMA_BLOCK_ALLOCATOR_TYPE == 1 || NUMA_BLOCK_ALLOCATOR_TYPE == 3
#if SWEET_THREADING || SWEET_REXI_THREAD_PARALLEL_SUM
	#pragma omp critical
#endif
	#endif
		{
			std::vector<void*>& block_list = getBlocksSameSize(i_size);
			block_list.push_back(i_data);
		}

#endif
	}
};


#endif /* SRC_INCLUDE_SWEET_NUMABLOCKALLOC_HPP_ */
