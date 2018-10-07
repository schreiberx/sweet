/*
 * REXI_SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_REXI_SIMULATION_VARIABLES_HPP_
#define SRC_REXI_SIMULATION_VARIABLES_HPP_

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <sweet/FatalError.hpp>



/**
 * REXI
 */
struct REXI_SimulationVariables
{
	/**
	 * Choose REXI solver method
	 */
	std::string rexi_method = "ci";


	/**
	 * Use only half of the poles for REXI
	 */
	bool use_half_poles = false;

	/**
	 * Extend modes for certain operations
	 */
	int use_sphere_extended_modes = 2;

	/**
	 * Normalize REXI for geostrophic balance
	 */
	bool normalization = true;

	/**
	 * REXI beta cutoff
	 * Remove REXI terms if the absolute value of the
	 * beta coefficients are below this value
	 */
	double beta_cutoff = 0.0;

	/**
	 * Add U0 (initial condition of current time step) to REXI terms
	 */
	bool rexi_terms_add_u0 = 0;

	/**
	 * Use REXI preallocation
	 */
	bool sphere_solver_preallocation = true;

	/**
	 * Use direct solution instead of REXI
	 */
	bool use_direct_solution = false;


	/***************************************************
	 * REXI Terry
	 */

	/**
	 * REXI parameter h
	 */
	double terry_h = -1.0;

	/**
	 * REXI parameter M: Number of Gaussian basis functions.
	 * This is not the number of rational functions!
	 */
	int terry_M = 128;

	/**
	 * REXI parameter L
	 */
	int terry_L = 0;

	/***************************************************
	 * REXI File
	 */

	std::string file_filename_alpha = "";
	std::string file_filename_beta = "";


#if 0
	/**
	 * integration range for which NGREXI should be valid
	 * [-test_range_abs;test_range_abs]
	 */
	double file_test_min = 0;
	double file_test_max = 0;

	/*
	 * Filename to directly load REXI coefficients
	 */
	std::string file_filename = "";

	std::string file_faf_dir = "./data/faf_data";

	/*
	 * Number of rational functions
	 */
	int file_N = 0;

	/*
	 * Spacing of rational functions
	 */
	double file_h = -1.0;

	/**
	 * Max double precision error within test region
	 */
	double file_max_error_double_precision = 1e-12;
#endif

	/***************************************************
	 * REXI CI
	 */

	/*
	 * Number of quadrature points
	 */
	int ci_n = 16;

	/*
	 * Primitive name
	 */
	std::string ci_primitive = "circle";

	/*
	 * Maximum real (positive) Eigenvalue enscribed by circle.
	 */
	double ci_max_real = -999;
	double ci_max_imag = -999;

	/*
	 * Size of primitive in real axis
	 */
	double ci_s_real = 5;

	/*
	 * Size of primitive in imag axis
	 */
	double ci_s_imag = 5;

	/*
	 * Shift
	 */
	double ci_mu = 0;


	/*
	 * Gaussian filter
	 *
	 * exp( -pow(dt_norm / (a * dt) * x, N) )
	 *
	 */
//	double ci_gaussian_filter = 0.0;

	/*
	 * Parameter to control the filter spread
	 */
	double ci_gaussian_filter_scale_a = 1.0;

	/*
	 * Timestep size to scale the spectrum to the interval [-i;i]
	 *
	 * As an approximation, this timestep size is similar to the max. stable
	 * time step size over a few time steps for RK1 and RK2.
	 *
	 * A value of 0 deactivates the filter
	 */
	double ci_gaussian_filter_dt_norm = 0.0;

	/*
	 * Exponent of the Gaussian filter
	 */
	double ci_gaussian_filter_exp_N = 0.0;



	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << " + rexi_method: " << rexi_method << std::endl;
		std::cout << " + use_half_poles: " << use_half_poles << std::endl;
		std::cout << " + rexi_normalization: " << normalization << std::endl;
		std::cout << " + use_direct_solution: " << use_direct_solution << std::endl;
		std::cout << " + beta_cutoff: " << beta_cutoff << std::endl;
		std::cout << " + use_extended_modes: " << use_sphere_extended_modes << std::endl;
		std::cout << " + rexi_sphere_solver_preallocation: " << sphere_solver_preallocation << std::endl;
		std::cout << " [REXI Terry]" << std::endl;
		std::cout << " + h: " << terry_h << std::endl;
		std::cout << " + M: " << terry_M << std::endl;
		std::cout << " + L: " << terry_L << std::endl;
		std::cout << " [REXI File]" << std::endl;

#if 0
		std::cout << " + file_faf_dir: " << file_faf_dir << std::endl;
		std::cout << " + file_N: " << file_N << std::endl;
		std::cout << " + file_h: " << file_h << std::endl;
		std::cout << " + file_test_min: " << file_test_min << std::endl;
		std::cout << " + file_test_max: " << file_test_max << std::endl;
		std::cout << " + file_max_error_double_precision: " << file_max_error_double_precision << std::endl;
#endif

		std::cout << " + file_filename_alpha: " << file_filename_alpha << std::endl;
		std::cout << " + file_filename_beta: " << file_filename_beta << std::endl;
		std::cout << " [REXI CI]" << std::endl;
		std::cout << " + ci_n: " << ci_n << std::endl;
		std::cout << " + ci_primitive: " << ci_primitive << std::endl;
		std::cout << " + ci_max_real: " << ci_max_real << std::endl;
		std::cout << " + ci_max_imag: " << ci_max_imag << std::endl;
		std::cout << " + ci_s_real: " << ci_s_real << std::endl;
		std::cout << " + ci_s_imag: " << ci_s_imag << std::endl;
		std::cout << " + ci_mu: " << ci_mu << std::endl;
		std::cout << " + ci_gaussian_filter_scale_a: " << ci_gaussian_filter_scale_a << std::endl;
		std::cout << " + ci_gaussian_filter_dt_norm: " << ci_gaussian_filter_dt_norm << std::endl;
		std::cout << " + ci_gaussian_filter_exp_N: " << ci_gaussian_filter_exp_N << std::endl;
		std::cout << std::endl;
	}



	void outputProgParams()
	{
		std::cout << "" << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << "	--rexi-method [str]	Choose REXI method ('terry', 'file'), default:0" << std::endl;
		std::cout << "	--rexi-half [bool]			Use half REXI poles, default:1" << std::endl;
		std::cout << "	--rexi-normalization [bool]		Use REXI normalization around geostrophic balance, default:1" << std::endl;
		std::cout << "	--rexi-sphere-preallocation [bool]	Use preallocation of SPH-REXI solver coefficients, default:1" << std::endl;
		std::cout << "	--rexi-ext-modes [int]	Use this number of extended modes in spherical harmonics" << std::endl;
		std::cout << "	--rexi-use-direct-solution [bool]	Use direct solution (analytical) for REXI, default:0" << std::endl;
		std::cout << "	--rexi-beta-cutoff [float]	Cutoff of REXI coefficients if the absolute value of beta exceed this value, default:0" << std::endl;
		std::cout << "	--rexi-add-u0 [bool]	Add U0 to the sum of REXI terms, default:0" << std::endl;
		std::cout << "  REXI Terry:" << std::endl;
		std::cout << "	--rexi-terry-h [float]			REXI parameter h" << std::endl;
		std::cout << "	--rexi-terry-m [int]				REXI parameter M" << std::endl;
		std::cout << "	--rexi-terry-l [int]				REXI parameter L" << std::endl;
		std::cout << "  REXI File:" << std::endl;
		std::cout << "	--rexi-file-alpha [string]			File with alpha coefficients" << std::endl;
		std::cout << "	--rexi-file-beta [string]			File with beta coefficients" << std::endl;
#if 0
		std::cout << "	--rexi-file-faf-dir [string]			Directory with FAF coefficients" << std::endl;
		std::cout << "	--rexi-file-N [int]		Number of rational basis functions, default:auto" << std::endl;
		std::cout << "	--rexi-file-h [double]	Spacing of rational basis functions, default:auto" << std::endl;
		std::cout << "	--rexi-file-test-min [double]	Set minimum test interval, default:0" << std::endl;
		std::cout << "	--rexi-file-test-max [double]	Set maximum test interval, default:0" << std::endl;
		std::cout << "	--rexi-file-test-abs [double]	Set min/max test interval, default:0" << std::endl;
		std::cout << "	--rexi-file-max-error [double]	Maximum allowed error within test interval, default:0" << std::endl;
		std::cout << "	--rexi-file-filename [string]	Filename of REXI coefficients, default:''" << std::endl;
#endif
		std::cout << "  REXI CI:" << std::endl;
		std::cout << "	--rexi-ci-n [double]	Number of quadrature points, default: 64" << std::endl;
		std::cout << "	--rexi-ci-primitive [string]	Primitive ('circle', 'rectangle'), default: 'circle'" << std::endl;
		std::cout << "	--rexi-ci-max-real [double]	Maximum real Eigenvalue to shift circle, default: -999" << std::endl;
		std::cout << "	--rexi-ci-max-imag [double]	Maximum imag Eigenvalue to shift circle, default: -999" << std::endl;
		std::cout << "	--rexi-ci-sx [double]	Size of primitive in real, default: 1" << std::endl;
		std::cout << "	--rexi-ci-sy [double]	Size of primitive in imag, default: 1" << std::endl;
		std::cout << "	--rexi-ci-mu [double]	Shift, default: 0" << std::endl;
		std::cout << "	--rexi-ci-gaussian-filter-scale [double]	Gaussian filter scaling parameter, default: 0 (no filter)" << std::endl;
		std::cout << "	--rexi-ci-gaussian-filter-dt-norm [double]	Gaussian filter normalization parameter (max timestep size for explicit method), default: 0 (no filter)" << std::endl;
		std::cout << "	--rexi-ci-gaussian-filter-exp-N [double]	Gaussian filter exponent, larger values lead to sharper function form, default: 2 (Gaussian-shaped)" << std::endl;
		std::cout << "" << std::endl;
	}



	void setup_longOptionList(
			struct option io_long_options[],		///< string and meta information for long options
			int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
			int i_max_options					///< maximum number of options
	)
	{
		io_long_options[io_next_free_program_option] = {"rexi-terry-h", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-m", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-l", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-half", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-normalization", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-sphere-preallocation", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-use-direct-solution", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ext-modes", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-beta-cutoff", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-add-u0", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-method", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

#if 0
		io_long_options[io_next_free_program_option] = {"rexi-file-faf-dir", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-n", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-h", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-test-min", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-test-max", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-test-abs", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-max-error", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;
#endif

		io_long_options[io_next_free_program_option] = {"rexi-file-alpha", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-file-beta", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-n", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-primitive", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-max-real", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-max-imag", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-sx", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-sy", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-mu", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-gaussian-filter-scale", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-gaussian-filter-dt-norm", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ci-gaussian-filter-exp-N", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;
	}

	/**
	 * Callback method to setup the values for the option with given index.
	 *
	 * \return Number of processed options or 0 in case of processed arguments
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
			case 0:	terry_h = atof(optarg);	return 0;
			case 1:	terry_M = atoi(optarg);	return 0;
			case 2:	terry_L = atoi(optarg);	return 0;
			case 3:	use_half_poles = atoi(optarg);	return 0;
			case 4:	normalization = atoi(optarg);	return 0;
			case 5:	sphere_solver_preallocation = atoi(optarg);	return 0;
			case 6:	use_direct_solution = atoi(optarg);	return 0;
			case 7:	use_sphere_extended_modes = atoi(optarg);	return 0;
			case 8: beta_cutoff = atof(optarg);	return 0;
			case 9: rexi_terms_add_u0 = atoi(optarg);	return 0;

			case 10:		rexi_method = optarg;	return 0;
#if 0
			case 10:		file_faf_dir = optarg;	return 0;
			case 11:	file_N = atoi(optarg);	return 0;
			case 12:	file_h = atof(optarg);	return 0;
			case 13:	file_test_min = atof(optarg);	return 0;
			case 14:	file_test_max = atof(optarg);	return 0;
			case 15:	file_test_max = atof(optarg);	file_test_min = -file_test_max;	return 0;
			case 16:	file_max_error_double_precision = atof(optarg);	return 0;
#endif
			case 11:		file_filename_alpha = optarg;	return 0;
			case 12:	file_filename_beta = optarg;	return 0;

			case 13:	ci_n = atoi(optarg);	return 0;
			case 14:	ci_primitive = optarg;	return 0;
			case 15:	ci_max_real = atof(optarg);	return 0;
			case 16:	ci_max_imag = atof(optarg);	return 0;
			case 17:	ci_s_real = atof(optarg);	return 0;
			case 18:	ci_s_imag = atof(optarg);	return 0;
			case 19:	ci_mu = atof(optarg);	return 0;
			case 20:	ci_gaussian_filter_scale_a = atof(optarg);	return 0;
			case 21:	ci_gaussian_filter_dt_norm = atof(optarg);	return 0;
			case 22:	ci_gaussian_filter_exp_N = atof(optarg);	return 0;
		}


		if (rexi_method != "" && rexi_method == "terry" && rexi_method == "file")
			FatalError("Invalid argument for '--rexi-method='");

		return 23;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
