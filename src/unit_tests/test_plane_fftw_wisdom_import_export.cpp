#include <iostream>
#include <fftw3.h>
#include <omp.h>
#include <cstdlib>
#include <sweet/SimulationVariables.hpp>


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
			// estimation base don workload

			if (num_cells < 32*32)
				//flags |= FFTW_EXHAUSTIVE;
				num_cells |= FFTW_MEASURE;
			else if (num_cells < 128*128)
				num_cells |= FFTW_MEASURE;
			else
				num_cells |= FFTW_PATIENT;

			if (i_reuse_spectral_transformation_plans == 2)
			{
				std::cout << "Enforcing to use Wisdom" << std::endl;
				flags |= FFTW_WISDOM_ONLY;
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


int main(int i_argc, char **i_argv)
{
	SimulationVariables simVars;

	simVars.setupFromMainParameters(i_argc, i_argv, nullptr, true);
	simVars.outputConfig();

#if SWEET_THREADING
	int nthreads = omp_get_max_threads();
	std::cout << " + nthreads: " << nthreads << std::endl;
#else
	int nthreads = 0;
#endif

	TestFFTPlans t;

	t.run(simVars.misc.reuse_spectral_transformation_plans, simVars.disc.space_res_physical, nthreads);

	std::cout << "FIN" << std::endl;

	return 0;
}
