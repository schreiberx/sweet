/*
 * REXI_SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_REXI_SIMULATION_VARIABLES_HPP_
#define SRC_REXI_SIMULATION_VARIABLES_HPP_

#include <vector>
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
	 * Extend modes for certain operations
	 */
	int use_sphere_extended_modes = 2;

	/**
	 * Use REXI preallocation
	 */
	bool sphere_solver_preallocation = true;


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

	/**
	 * Reduction to half
	 */
	int terry_reduce_to_half = 0;


	/**
	 * Normalization for geostrophic modes
	 */
	int terry_normalization = 0;




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



	/***************************************************
	 * REXI FILES
	 */

	/*
	 * String for REXI coefficients
	 *
	 * Format:
	 * [[function_name:]filepath],[[function_name:]filepath]
	 *
	 * with function_name \in {"phi0", "phi1", ..., "varphi0", "varphi1", ...}
	 * default value for function_name is "phi0"
	 *
	 * E.g.
	 * phi0:rexi_phi0.data,phi1:rexi_phi1.data
	 */
	std::string rexi_files;


public:
	class REXIFile
	{
	public:
		std::string function_name;
		std::string filename;

		REXIFile(
				const std::string &i_function_name,
				const std::string &i_filename
		)	:
			function_name(i_function_name),
			filename(i_filename)
		{
		}
	};

public:
	std::vector<REXIFile> p_rexi_files_processed;


	void outputConfig()
	{
		std::cout << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << " + rexi_method: " << rexi_method << std::endl;
		std::cout << "REXI generic parameters:" << std::endl;
		std::cout << " + use_extended_modes: " << use_sphere_extended_modes << std::endl;
		std::cout << " + rexi_sphere_solver_preallocation: " << sphere_solver_preallocation << std::endl;

		std::cout << " [REXI Files]" << std::endl;
		std::cout << " + rexi_files: " << rexi_files << std::endl;

		for (std::vector<REXIFile>::iterator iter = p_rexi_files_processed.begin(); iter != p_rexi_files_processed.end(); iter++)
			std::cout << "     \"" << iter->function_name << "\" : \"" << iter->filename << "\"" << std::endl;

		std::cout << " [REXI Terry]" << std::endl;
		std::cout << " + h: " << terry_h << std::endl;
		std::cout << " + M: " << terry_M << std::endl;
		std::cout << " + L: " << terry_L << std::endl;
		std::cout << " [REXI CI]" << std::endl;
		std::cout << " + ci_n: " << ci_n << std::endl;
		std::cout << " + ci_primitive: " << ci_primitive << std::endl;
		std::cout << " + ci_max_real: " << ci_max_real << std::endl;
		std::cout << " + ci_max_imag: " << ci_max_imag << std::endl;
		std::cout << " + ci_s_real: " << ci_s_real << std::endl;
		std::cout << " + ci_s_imag: " << ci_s_imag << std::endl;
		std::cout << " + ci_mu: " << ci_mu << std::endl;
		std::cout << std::endl;
	}



	void outputProgParams()
	{
		std::cout << "" << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << "	--rexi-method [str]	Choose REXI method ('terry', 'file'), default:0" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-use-direct-solution [bool]	Use direct solution (analytical) for REXI, default:0" << std::endl;
		std::cout << "	--rexi-sphere-preallocation [bool]	Use preallocation of SPH-REXI solver coefficients, default:1" << std::endl;
		std::cout << "	--rexi-ext-modes [int]	Use this number of extended modes in spherical harmonics" << std::endl;
		std::cout << std::endl;
		std::cout << "  REXI file interface:" << std::endl;
		std::cout << "	--rexi-files [str]	REXI files: [function_name0:]filepath0,[function_name1:]filepath1,..." << std::endl;
		std::cout << std::endl;
		std::cout << "  T-REXI (Terry et al. method) [DEPRECATED interface]:" << std::endl;
		std::cout << "	--rexi-terry-h [float]			T-REXI parameter h" << std::endl;
		std::cout << "	--rexi-terry-m [int]			T-REXI parameter M" << std::endl;
		std::cout << "	--rexi-terry-l [int]			T-REXI parameter L" << std::endl;
		std::cout << "	--rexi-terry-reduce-to-half [int]	T-REXI reduction to half" << std::endl;
		std::cout << "	--rexi-terry-normalization [int]	T-REXI normalization" << std::endl;
		std::cout << std::endl;
		std::cout << "  CI-REXI (Cauchy Contour Integral-based method) [DEPRECATED interface]:" << std::endl;
		std::cout << "	--rexi-ci-n [double]	Number of quadrature points, default: 64" << std::endl;
		std::cout << "	--rexi-ci-primitive [string]	Primitive ('circle', 'rectangle'), default: 'circle'" << std::endl;
		std::cout << "	--rexi-ci-max-real [double]	Maximum real Eigenvalue to shift circle, default: -999" << std::endl;
		std::cout << "	--rexi-ci-max-imag [double]	Maximum imag Eigenvalue to shift circle, default: -999" << std::endl;
		std::cout << "	--rexi-ci-sx [double]	Size of primitive in real, default: 1" << std::endl;
		std::cout << "	--rexi-ci-sy [double]	Size of primitive in imag, default: 1" << std::endl;
		std::cout << "	--rexi-ci-mu [double]	Shift, default: 0" << std::endl;
		std::cout << "" << std::endl;
	}



	void setup_longOptionList(
			struct option io_long_options[],		///< string and meta information for long options
			int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
			int i_max_options					///< maximum number of options
	)
	{
		io_long_options[io_next_free_program_option] = {"rexi-method", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;


		// Generic REXI options
//		io_long_options[io_next_free_program_option] = {"rexi-use-direct-solution", required_argument, 0, 256+io_next_free_program_option};
//		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-sphere-preallocation", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-ext-modes", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;


		// Files
		io_long_options[io_next_free_program_option] = {"rexi-files", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;


		// Terry
		io_long_options[io_next_free_program_option] = {"rexi-terry-h", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-m", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-l", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-reduce-to-half", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"rexi-terry-normalization", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;


		// Cauchy
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
	}



	/*
	 * Source: https://www.techiedelight.com/split-string-cpp-using-delimiter/
	 */
	static
	void split_string(
			std::string const &i_string,
			const char needle,
			std::vector<std::string> &o_split_string
	)
	{
		std::size_t start;
		std::size_t end = 0;

		while ((start = i_string.find_first_not_of(needle, end)) != std::string::npos)
		{
			end = i_string.find(needle, start);
			o_split_string.push_back(i_string.substr(start, end - start));
		}
	}



	/**
	 * Callback method to setup the values for the option with given index.
	 *
	 * \return Number of processed options or 0 in case of no processed arguments
	 */
	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		bool rexi_files_given = false;
		switch(i_option_index)
		{
			case 0:		rexi_method = optarg;	return 0;

			case 1:		sphere_solver_preallocation = atoi(optarg);	return 0;
			case 2:		use_sphere_extended_modes = atoi(optarg);	return 0;

			case 3:		rexi_files = optarg;	rexi_files_given = true; break;

			case 4:		terry_h = atof(optarg);	return 0;
			case 5:		terry_M = atoi(optarg);	return 0;
			case 6:		terry_L = atoi(optarg);	return 0;
			case 7:		terry_reduce_to_half = atoi(optarg);	return 0;
			case 8:		terry_normalization = atoi(optarg);	return 0;

			case 9:	ci_n = atoi(optarg);	return 0;
			case 10:	ci_primitive = optarg;	return 0;
			case 11:	ci_max_real = atof(optarg);	return 0;
			case 12:	ci_max_imag = atof(optarg);	return 0;
			case 13:	ci_s_real = atof(optarg);	return 0;
			case 14:	ci_s_imag = atof(optarg);	return 0;
			case 15:	ci_mu = atof(optarg);	return 0;
		}

		if (rexi_files_given)
		{
			p_rexi_files_processed.empty();

			std::vector<std::string> split_rexi_files;
			split_string(
					rexi_files,
					',',
					split_rexi_files
				);

			for (std::size_t i = 0; i < split_rexi_files.size(); i++)
			{
				std::vector<std::string> split_split_rexi_files;
				split_string(
						split_rexi_files[i],
						':',
						split_split_rexi_files
					);

				if (split_split_rexi_files.size() == 1)
				{
					p_rexi_files_processed.push_back(REXIFile("", split_rexi_files[i]));
				}
				else if (split_split_rexi_files.size() == 2)
				{
					p_rexi_files_processed.push_back(REXIFile(split_split_rexi_files[0], split_split_rexi_files[1]));
				}
				else
				{
					std::cerr << "Unable to process REXI file parameters" << std::endl;
					std::cerr << "Current parameter: " << rexi_files << std::endl;
					std::cerr << "Current REXI coefficient file: " << split_rexi_files[i] << std::endl;
					FatalError("Invalid format");
				}
			}

			return 0;
		}

		if (rexi_method != "" && rexi_method == "terry" && rexi_method == "file")
			FatalError("Invalid argument for '--rexi-method='");

		return 16;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
