/*
 * MemBlockAlloc.hpp
 *
 *  Created on: 14 Sep 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * Changelog:
 *   - 2021-12-23: Made fully configurable via environment variable
 *   - 2022-01-08: Various updates to help finding bugs, more information if used with help
 *
 */
#ifndef INCLUDE_MEMBLOCKALLOC_NEW_HPP_
#define INCLUDE_MEMBLOCKALLOC_NEW_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>

#include <stdexcept>
#include <cstdlib>
#include <cassert>

#include "StringSplit.hpp"


/**
 * Prefix string for console output
 */
#define MEMBLOCKALLOC_PREFIX	"[MEMBLOCKALLOC] "


/**
 * Support allocator with NUMA granularity
 */
#ifndef MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
	#define MEMBLOCKALLOC_ENABLE_NUMA_ALLOC 1
#endif

#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
	#include <numa.h>
#endif



/*
 * SWEET specific part
 */
#if SWEET_THREADING_SPACE == 1 || SWEET_THREADING_TIME_REXI == 1
	#define MEMBLOCKALLOC_ENABLE_OMP 1
#endif

#if SWEET_THREADING_SPACE == 0 && SWEET_THREADING_TIME_REXI == 0
	#define MEMBLOCKALLOC_ENABLE_OMP 0
#endif

#ifndef MEMBLOCKALLOC_ENABLE_OMP
	#define MEMBLOCKALLOC_ENABLE_OMP 1
#endif

/*
 * Debug mode with extra sanity checks
 */
#ifndef MEMBLOCKALLOC_DEBUG
	#define MEMBLOCKALLOC_DEBUG 0
#endif


/*
 * SWEET specific part
 */
#if MEMBLOCKALLOC_ENABLE_OMP
	#include <sweet/core/openmp_helper.hpp>
	#include <omp.h>
#endif


#define MEMBLOCKALLOC_MODE__SYSTEM 0
#define MEMBLOCKALLOC_MODE__ONE 1
#define MEMBLOCKALLOC_MODE__PERTHREAD 2
#define MEMBLOCKALLOC_MODE__PERNUMA 3

#if SWEET_DEBUG == 1
	#undef MEMBLOCKALLOC_DEBUG
	#define MEMBLOCKALLOC_DEBUG 1
#endif


/**
 * This class implements a memory manager with various features:
 *  - caches the allocation of large memory blocks
 *  - first touch policy (optional)
 *  - highly performing if only the same block sizes are required
 */
class MemBlockAlloc
{
	static const char *get_helptext()
	{
		return
		 " Use environment variable\n"
		 " 		MEMBLOCKALLOC=[option1=...][,option2=...][...]\n"
		 " to configure.\n"
		 "\n"
		 " Options:\n"
		 "\n"
		 " 	verbose=[int]\n"
		 " 		0: Disable verbose mode\n"
		 " 		1: Enable verbose mode during initialization\n"
		 " 		9: Print information on allocating / releasing memory\n"
		 "\n"
		 " 	firsttouch=[int]\n"
		 " 		0: Disabled\n"
		 " 		1: Enabled, threaded with 'omp parallel for'\n"
		 " 		2: Enabled, nonthreaded\n"
		 "\n"
		 "	alloc={system,one"
#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
				",pernuma"
				",perthread"
#endif
		 "}:\n"
		 "  	'system': Use system's allocator!\n"
		 "\n"
		 "		'one': Allocate one memory block chain\n"
		 "    		- Performance: Requires finer granular synchronization\n"
		 "    		- Hint: This was useful for XeonPhi KNC (not KNL)\n"
		 "\n"

#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
		 "  	'perthread': Allocate one memory block chain per thread!\n"
		 "    		- Performance: Requires no synchronization\n"
		 "\n"
		 "  	'pernuma': Allocate one memory block chain per NUMA domain\n"
		 "    		- Performance: Requires additional synchronization (critical regions)\n"
		 "\n"
#endif
		;
	};

	/**
	 * Block allocation mode
	 */
	int mem_block_allocation_mode = MEMBLOCKALLOC_MODE__ONE;

	/**
	 * First touch policy
	 * 0: disabled
	 * 1: enabled, threaded "omp parallel for"
	 * 2: enabled, no threading
	 */
	#if MEMBLOCKALLOC_ENABLE_OMP
		int first_touch_policy = 1;
	#else
		int first_touch_policy = 2;
	#endif

	/**
	 * verbosity
	 */
	int verbosity_level = 1;

	/**
	 * Number of allocation domains
	 *
	 * Either this is the number of NUMA domains or the number of threads,
	 * depending on NUMA_BLOCK_ALLOCATOR_TYPE
	 */
	int _num_block_chain_domains = 1;


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
	std::vector<DomainMemBlocks> _domain_block_groups;


	/**
	 * setup already executed?
	 */
	bool _setup_done = false;



private:
	inline
	static
	int& getThreadLocalDomainIdRef()
	{
		/**
		 * Domain node for current thread
		 */
		// WARNING: This class has to be compiled always with -fopenmp activated
		static thread_local int domain_id;
		return domain_id;
	}

	static
	void fatal_error(const std::string &str)
	{
		std::cerr << MEMBLOCKALLOC_PREFIX "ERROR: " << str << std::endl;
		assert(false);
		throw std::runtime_error(str);
		std::exit(1);
	}

	void print_configuration()
	{
		std::cout << MEMBLOCKALLOC_PREFIX "*** MemBlockAlloc VERBOSE information ***" << std::endl;
		std::cout << MEMBLOCKALLOC_PREFIX " + verbosity_level: " << verbosity_level << std::endl;
		std::cout << MEMBLOCKALLOC_PREFIX " + first_touch_policy: " << first_touch_policy << std::endl;
		std::cout << MEMBLOCKALLOC_PREFIX " + mem_block_allocation_mode: ";

		if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__SYSTEM)
			std::cout << "system";
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__ONE)
			std::cout << "one";
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERTHREAD)
			std::cout << "perthread";
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA)
			std::cout << "pernuma";
		else
			fatal_error("Internal error (invalid mode enum)");

		std::cout << std::endl;
	}

	/*
	 * Parse environment variable (if it exists) and set parameters
	 */
	void parse_argv()
	{
		char *val = std::getenv("MEMBLOCKALLOC");

		if (val == nullptr)
		{
			/*
			 * Nothing to do
			 */
			return;
		}

		const std::string envstring(val);

		/**
		 * Split comma separated parameters, e.g.,
		 * 	 verbose=1,alloc=perthread
		 *
		 * List with items
		 *   verbose=1
		 *   alloc=perthread
		 */
		std::vector<std::string> params = StringSplit::split(envstring, ",");

		for (auto iter = params.begin(); iter != params.end(); iter++)
		{
			const std::string &param = *iter;


			std::vector<std::string> split_params = StringSplit::split(param, "=");

			if (split_params.size() == 0)
				fatal_error(std::string("Error with parameter")+param);

			if (split_params[0] == "help")
			{
				print_configuration();
				std::cout << std::endl;
				std::cout << get_helptext() << std::endl;
				std::cout << std::endl;
				std::exit(1);
			}
			else if (split_params[0] == "verbose")
			{
				/*
				 * Parse, e.g.,
				 * 	verbose=13
				 */
				if (split_params.size() > 2)
					fatal_error(std::string("Parameter missing for ")+split_params[0]);

				if (split_params.size() == 1)
					verbosity_level = 99;
				else
					verbosity_level = std::atoi(split_params[1].c_str());
			}
			else if (split_params[0] == "firsttouch")
			{
				/*
				 * Parse, e.g.,
				 * 	firsttouch=1
				 */
				if (split_params.size() != 2)
					fatal_error(std::string("first touch option must have exactly one parameter, given ")+param);

				first_touch_policy = std::atoi(split_params[1].c_str());

				if (first_touch_policy != 0 && first_touch_policy != 1 && first_touch_policy != 2)
					fatal_error(std::string("first touch policy must be set to 0, 1 or 2"));

			}
			else if (split_params[0] == "alloc")
			{
				/*
				 * Parse, e.g.,
				 * 	alloc=perthread
				 */
				if (split_params.size() != 2)
					fatal_error(std::string("Only one parameter required for 'alloc='"));

				if (split_params[1] == "system")
				{
					mem_block_allocation_mode = MEMBLOCKALLOC_MODE__SYSTEM;
				}
				else if (split_params[1] == "one")
				{
					mem_block_allocation_mode = MEMBLOCKALLOC_MODE__ONE;
				}
#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
				else if (split_params[1] == "perthread")
				{
					mem_block_allocation_mode = MEMBLOCKALLOC_MODE__PERTHREAD;
				}
				else if (split_params[1] == "pernuma")
				{
					mem_block_allocation_mode = MEMBLOCKALLOC_MODE__PERNUMA;
				}
#endif
				else
				{
					fatal_error(std::string("Unknown parameter '") + split_params[1] + ("' for parameter alloc=..."));
				}
			}
			else
			{
				fatal_error(std::string("Unknown option '") + split_params[0] + "'");
			}
		}

		if (verbosity_level >= 1)
			print_configuration();

		return;
	}


public:
	/**
	 * Constructor.
	 *
	 * To be used directly at the beginning of 'main'
	 */
	MemBlockAlloc(int secret_code)	:
		_setup_done(false)
	{
		if (_setup_done)
			fatal_error("Setup in MemBlockAlloc called twice!");

		if (secret_code != 1337)
			fatal_error("Secret code mismatch, do NOT call MemBlockAlloc(...) on your own!");

		parse_argv();

		if (verbosity_level > 1)
		{
			std::cout << MEMBLOCKALLOC_PREFIX << "MemBlockAlloc() called (constructor, should be called only once)" << std::endl;
		}

#if MEMBLOCKALLOC_ENABLE_OMP == 0
		if (	mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA		||
			mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERTHREAD	)
		{
			std::cerr << MEMBLOCKALLOC_PREFIX "WARNING: Using thread-supported memory allocator, but without OMP enabled." << std::endl;
		}
#endif

		if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__SYSTEM)
		{
			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "NUMA block alloc: Using default system's posix_memalign/free allocator" << std::endl;

			/*
			 * this is the system's default allocator
			 */
			_num_block_chain_domains = 1;
			getThreadLocalDomainIdRef() = 0;
		}
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__ONE)
		{
			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "NUMA block alloc: Using just a single memory block chain" << std::endl;

			_num_block_chain_domains = 1;
			getThreadLocalDomainIdRef() = 0;
		}
#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERTHREAD)
		{
			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "NUMA block alloc: Using allocator based on thread granularity" << std::endl;

			/*
			 * Thread granularity, use this also per default
			 */
			#if MEMBLOCKALLOC_ENABLE_OMP
				_num_block_chain_domains = omp_get_max_threads();
			#else
				_num_block_chain_domains = 1;
			#endif

			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "num_block_chain_domains: " << _num_block_chain_domains << std::endl;

			#if MEMBLOCKALLOC_ENABLE_OMP
				getThreadLocalDomainIdRef() = omp_get_thread_num();

				#pragma omp parallel
				{
					getThreadLocalDomainIdRef() = omp_get_thread_num();
				}
			#else
				getThreadLocalDomainIdRef() = 0;
			#endif

		}
		else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA)
		{
			/*
			 * Allocator which works per NUMA domain
			 */
			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "NUMA block alloc: Using NUMA node granularity" << std::endl;

			/*
			 * NUMA granularity
			 */
			_num_block_chain_domains = numa_num_configured_nodes();
			if (verbosity_level > 0)
				std::cout << MEMBLOCKALLOC_PREFIX "_num_block_chain_domains: " << _num_block_chain_domains << std::endl;

			// set NUMA id in case that master thread has a different id than the first thread
			int cpuid = sched_getcpu();
			getThreadLocalDomainIdRef() = numa_node_of_cpu(cpuid);

			#if MEMBLOCKALLOC_ENABLE_OMP
				#pragma omp parallel
				{
					int cpuid = sched_getcpu();
					getThreadLocalDomainIdRef() = numa_node_of_cpu(cpuid);
				}
			#else
				getThreadLocalDomainIdRef() = 0;
			#endif
		}
#endif
		else
		{
			fatal_error("Internal error (allocation mode enum)");
		}

		#if MEMBLOCKALLOC_ENABLE_OMP
			if (verbosity_level >= 10)
			{
				#pragma omp parallel
				{
					#pragma omp critical
					{
						std::cout << MEMBLOCKALLOC_PREFIX " + thread_id " << omp_get_thread_num() << " is assigned to memory allocator domain " << getThreadLocalDomainIdRef() << std::endl;
					}
				}
			}
		#else
			if (verbosity_level >= 10)
			{
				std::cout << MEMBLOCKALLOC_PREFIX " + thread_id 0 is assigned to memory allocator domain " << getThreadLocalDomainIdRef() << std::endl;
			}
		#endif

		_domain_block_groups.resize(_num_block_chain_domains);

		_setup_done = true;
	}


	static
	inline
	void init()
	{
		getSingletonRef();
	}

public:
	static
	inline
	MemBlockAlloc& getSingletonRef()
	{
		/*
		 * Compilation details:
		 *
		 * this is compiled to a code such as
		 * '
		 * 		cmp memManager, 0
		 * 		jne return
		 * 		# memManager = NUMABlockAlloc()
		 * 		...
		 * 	return:
		 * 		ret ...
		 * '
		 * In other words, it's initialized on the fly during the first access
		 */
		static MemBlockAlloc memBlockAlloc(1337);
		return memBlockAlloc;
	}



	/**
	 * return a list of blocks with the same size
	 *
	 * If the list does not exist, add an empty one
	 *
	 * *** NOT THREAD SAFE ***
	 */
	static
	std::vector<void*>& getBlockListSameSize(
			std::size_t i_size				///< size of blocks
	)
	{
		MemBlockAlloc &n = MemBlockAlloc::getSingletonRef();

		assert(n.getThreadLocalDomainIdRef() < (int)n._domain_block_groups.size());

		std::vector<MemBlocksSameSize>& block_groups = n._domain_block_groups[n.getThreadLocalDomainIdRef()].block_groups;

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


	/**
	 * return a block of the given size or 0 if no block is already allocated
	 *
	 * *** NOT THREAD SAFE ***
	 */
	static
	void* getBlockSameSize(
			std::size_t i_size				///< size of blocks
	)
	{
		std::vector<void*>& block_list = getBlockListSameSize(i_size);

		void *data = nullptr;
		if (block_list.size() > 0)
		{
			data = (void*)block_list.back();
			block_list.pop_back();
		}

		return data;
	}


	/**
	 * Explicitly write data to the areas instead of relying on
	 * the program to apply a first touch policy
	 */
	template <typename T=void>
	static
	T *first_touch_init(
			T *i_data,
			std::size_t i_size
	)
	{
		int _first_touch_policy = getSingletonRef().first_touch_policy;

		char *data = (char*)i_data;
		if (_first_touch_policy == 0)
		{
			// nothing
		}
		else if (_first_touch_policy == 1)
		{
#if MEMBLOCKALLOC_ENABLE_OMP
			#pragma omp parallel for
#endif
			for (std::size_t i = 0; i < i_size; i++)
				data[i] = 0;
		}
		else if (_first_touch_policy == 2)
		{
			for (std::size_t i = 0; i < i_size; i++)
				data[i] = 0;
		}

		return i_data;
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

		int _mem_block_allocation_mode = getSingletonRef().mem_block_allocation_mode;

		if (_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__SYSTEM)
		{
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

			first_touch_init(data, i_size);
		}
		else if (_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__ONE)
		{
			#if MEMBLOCKALLOC_ENABLE_OMP
			#	pragma omp critical
			#endif
			{
				data = (T*)getBlockSameSize(i_size);
			}

			if (data == nullptr)
			{
				int retval = posix_memalign((void**)&data, 4096, i_size);
				if (retval != 0)
				{
					std::cerr << "Unable to allocate memory" << std::endl;
					assert(false);
					exit(-1);
				}

				first_touch_init(data, i_size);
			}
		}
#if MEMBLOCKALLOC_ENABLE_NUMA_ALLOC
		else if (_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERTHREAD)
		{
			data = (T*)getBlockSameSize(i_size);

			if (data == nullptr)
			{
				data = (T*)numa_alloc_local(i_size);

				first_touch_init(data, i_size);
			}
		}
		else if (_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA)
		{
			// use critical section since the block chain is shared by different threads of the same domain
			#if MEMBLOCKALLOC_ENABLE_OMP
			#	pragma omp critical
			#endif
			{
				data = (T*)getBlockSameSize(i_size);
			}

			if (data == nullptr)
			{
				// Allocate block
				data = (T*)numa_alloc_local(i_size);

				first_touch_init(data, i_size);
			}
		}
#endif
		else
		{
			fatal_error("ALLOC: mode not found");
			return nullptr;
		}

#if MEMBLOCKALLOC_DEBUG
		MemBlockAlloc &n = MemBlockAlloc::getSingletonRef();

		if (n.verbosity_level >= 100)
			std::cout << "ALLOC " << (long long)data << ", " << i_size << std::endl;
#endif
		return data;
	}



public:
	inline
	static
	void free(
			void *i_data,
			std::size_t i_size
	)
	{
		MemBlockAlloc &n = MemBlockAlloc::getSingletonRef();

		if (!n._setup_done)
		{
			fatal_error("free: _setup_done not set! Either not set up or deconstructor already called\nHint 1: Compile with MEMBLOCKALLOC_DEBUG=1\nHint 2: Execute by setting environment flag to MEMBLOCKALLOC=verbose=10 to get more information");
		}

#if MEMBLOCKALLOC_DEBUG

		if (n.verbosity_level >= 100)
			std::cout << "FREE " << (long long)i_data << ", " << i_size << std::endl;

		if (i_data == nullptr)
		{
			fatal_error("free: nullptr received, not allowed with MemBlock");
			return;
		}
#endif

		int _mem_block_allocation_mode = n.mem_block_allocation_mode;

		if (_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__SYSTEM)
		{
			::free(i_data);
		}
		else
		{
			if (	_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA ||
				_mem_block_allocation_mode == MEMBLOCKALLOC_MODE__ONE)
			{
				#if MEMBLOCKALLOC_ENABLE_OMP
					#pragma omp critical
				#endif
				{
					std::vector<void*>& block_list = getBlockListSameSize(i_size);
					block_list.push_back(i_data);
				}
			}
			else
			{
				std::vector<void*>& block_list = getBlockListSameSize(i_size);
				block_list.push_back(i_data);
			}
		}
	}

public:
	~MemBlockAlloc()
	{
		if (verbosity_level > 1)
		{
			std::cout << MEMBLOCKALLOC_PREFIX << "~MemBlockAlloc() called (deconstructor, should be called only once)" << std::endl;
		}

		getSingletonRef()._shutdown();

		if (verbosity_level > 1)
		{
			std::cout << MEMBLOCKALLOC_PREFIX << "~MemBlockAlloc() finished" << std::endl;
		}
	}


public:
	void _shutdown()
	{
		if (!_setup_done)
			fatal_error("Setup not executed, but shutdown requested");

		for (auto& n : _domain_block_groups)
		{
			for (auto& g : n.block_groups)
			{
//				std::cout << "cleaning up " << g.free_blocks.size() << " blocks of size " << g.block_size << std::endl;

				for (auto& b : g.free_blocks)
				{
					if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__SYSTEM)
					{
						::free(b);
					}
					else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__ONE)
					{
						::free(b);
					}
					else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERTHREAD)
					{
						numa_free(b, g.block_size);
					}
					else if (mem_block_allocation_mode == MEMBLOCKALLOC_MODE__PERNUMA)
					{
						numa_free(b, g.block_size);
					}
					else
					{
						fatal_error("Internal error (_shutdown)");
					}
				}

				// free all blocks
				g.free_blocks.clear();
			}
		}

		_setup_done = false;
	}

};


#endif
