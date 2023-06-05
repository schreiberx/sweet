/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <iostream>
#include <fftw3.h>
#include <omp.h>
#include <cstdlib>


class TestFFTPlans
{
public:
	bool importWisdom(int i_reuse_spectral_transformation_plans)
	{
		static const char *wisdom_file = "sweet_fftw";

		std::cout << "fftw_import_wisdom_from_filename(" << wisdom_file << ")" << std::endl;

		int wisdom_plan_loaded = fftw_import_wisdom_from_filename(wisdom_file);
		if (wisdom_plan_loaded == 0)
		{
			std::cerr << "Failed to load FFTW wisdom from file '" << wisdom_file << "'" << std::endl;
			if (i_reuse_spectral_transformation_plans == 2)
				exit(1);
		}

		std::cout << "WISDOM: " << fftw_export_wisdom_to_string() << std::endl;

		return true;
	}

	bool exportWisdom()
	{
		static const char *wisdom_file = "sweet_fftw";

		std::cout << "fftw_export_wisdom_to_filename(" << wisdom_file << ")" << std::endl;

		int wisdom_plan_stored = fftw_export_wisdom_to_filename(wisdom_file);
		if (wisdom_plan_stored == 0)
		{
			std::cerr << "Failed to store FFTW wisdom to file " << wisdom_file << std::endl;
			exit(1);
		}

		std::cout << "WISDOM: " << fftw_export_wisdom_to_string() << std::endl;

		return true;
	}

	bool printWisdom()
	{
		std::cout << "WISDOM: " << fftw_export_wisdom_to_string() << std::endl;

		return true;
	}


	int run(
			int i_reuse_spectral_transformation_plans,
			int i_res[2],
			int i_nthreads
	)
	{
#if SWEET_THREADING
		std::cout << "fftw_init_threads()" << std::endl;
		int retval = fftw_init_threads();
		if (retval == 0)
		{
			std::cerr << "ERROR: fftw_init_threads()" << std::endl;
			exit(1);
		}


		std::cout << "fftw_plan_with_nthreads(" << i_nthreads << ")" << std::endl;
		fftw_plan_with_nthreads(i_nthreads);
#endif

		importWisdom(i_reuse_spectral_transformation_plans);


		unsigned int num_cells = i_res[0]*i_res[1];


		unsigned flags = 0;

		if (i_reuse_spectral_transformation_plans == -1)
		{
			flags |= FFTW_ESTIMATE;
		}
		else
		{
			// estimation based on workload
			if (i_reuse_spectral_transformation_plans == 2)
			{
				std::cout << "Enforcing to use Wisdom" << std::endl;
				flags |= FFTW_WISDOM_ONLY;
			}
			else
			{
				if (num_cells < 64*64)
					//flags |= FFTW_EXHAUSTIVE;
					flags |= FFTW_MEASURE;
				else
					flags |= FFTW_PATIENT;
			}
		}

		// allocate more data than necessary for spectral space
		fftw_complex *data_spectral = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * i_res[0]*i_res[1]);

		// physical space data
		double *data_physical = (double*)fftw_malloc(sizeof(double)*2 * i_res[0]*i_res[1]);

		std::cout << "fftw_plan_dft_r2c_2d(...)" << std::endl;
		fftw_plan	fftw_plan_forward;
		fftw_plan_forward =
				fftw_plan_dft_r2c_2d(
					i_res[1],
					i_res[0],
					data_physical,
					(fftw_complex*)data_spectral,
					flags
				);

		printWisdom();

		if (fftw_plan_forward == nullptr)
		{
			std::cout << "Wisdom: " << fftw_export_wisdom_to_string() << std::endl;
			std::cerr << "Failed to get forward plan dft_r2c fftw" << std::endl;
			exit(-1);
		}

		exportWisdom();

#if SWEET_THREADING
		fftw_cleanup_threads();
#endif
		fftw_cleanup();

		return 0;
	}

};



#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>

int main(int i_argc, char **i_argv)
{
	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	sweet::ShackPlaneDataOps *shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

	shackProgArgDict.printShackData();


#if SWEET_THREADING
	int nthreads = omp_get_max_threads();
	std::cout << " + nthreads: " << nthreads << std::endl;
#else
	int nthreads = 0;
#endif

	TestFFTPlans t;

	t.run(
			shackPlaneDataOps->reuse_spectral_transformation_plans,
			shackPlaneDataOps->space_res_physical,
			nthreads
		);

	std::cout << "FIN" << std::endl;

	return 0;
}
