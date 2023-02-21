/*
 * REXI_SimulationVariables.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SRC_REXI_SIMULATION_VARIABLES_HPP_
#define SRC_REXI_SIMULATION_VARIABLES_HPP_

#include <vector>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <sweet/SWEETError.hpp>

#include <sweet/shacks/ShackInterface.hpp>


/**
 * Exponential integration
 */
struct EXP_SimulationVariables	:
	public sweet::ClassDictionaryInterface
{
	/**
	 * Choose EXP solver method
	 */
	std::string exp_method = "direct";

	/**
	 * Use REXI preallocation
	 */
	bool sphere_solver_preallocation = true;


	/***************************************************
	 * Taylor EXP
	 */

	/*
	 * Number of Taylor expansions
	 */
	int taylor_num_expansions = -1;


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

	/***************************************************
	 * Direct EXP
	 */
	int exp_direct_precompute_phin = 0;


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




	void printProgramArguments(const std::string& i_prefix = "")
	{

		std::cout << "" << std::endl;
		std::cout << "REXI:" << std::endl;
		std::cout << "	--rexi-method [str]	Choose REXI method ('terry', 'file', 'direct'), default:0" << std::endl;
		std::cout << std::endl;
		std::cout << "	--rexi-sphere-preallocation [bool]	Use preallocation of SPH-REXI solver coefficients, default:1" << std::endl;
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
		std::cout << "  EXP direct:" << std::endl;
		std::cout << "	--exp_direct-precompute-phin [int]	Precompute phin in direct exp. solver (only available for SWE on the plane), default:0" << std::endl;
		std::cout << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		bool rexi_files_given;

		i_pa.getArgumentValueByKey("--rexi-method", exp_method);	// depreated
		i_pa.getArgumentValueByKey("--exp-method", exp_method);

		i_pa.getArgumentValueByKey("--taylor_num_expansions", taylor_num_expansions);
		i_pa.getArgumentValueByKey("--rexi-sphere-preallocation", sphere_solver_preallocation);

		if (i_pa.getArgumentValueByKey("--rexi-files", rexi_files))
			rexi_files_given = true;

		/*
		 * T-REXI implemented in SWEET (deprecated!!!)
		 */
		i_pa.getArgumentValueByKey("--rexi-terry-h", terry_h);
		i_pa.getArgumentValueByKey("--rexi-terry-m", terry_M);
		i_pa.getArgumentValueByKey("--rexi-terry-l", terry_L);
		i_pa.getArgumentValueByKey("--rexi-terry-reduce-to-half", terry_reduce_to_half);
		i_pa.getArgumentValueByKey("--rexi-terry-normalization", terry_normalization);

		i_pa.getArgumentValueByKey("--rexi-ci-n", ci_n);
		i_pa.getArgumentValueByKey("--rexi-ci-primitive", ci_primitive);
		i_pa.getArgumentValueByKey("--rexi-ci-max-real", ci_max_real);
		i_pa.getArgumentValueByKey("--rexi-ci-max-imag", ci_max_imag);
		i_pa.getArgumentValueByKey("--rexi-ci-sx", ci_s_real);
		i_pa.getArgumentValueByKey("--rexi-ci-sy", ci_s_imag);
		i_pa.getArgumentValueByKey("--rexi-ci-mu", ci_mu);

		i_pa.getArgumentValueByKey("--exp-direct-precompute-phin", exp_direct_precompute_phin);

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
					SWEETError("Invalid format");
				}
			}

			return -1;
		}

		if (exp_method != "" && exp_method == "terry" && exp_method == "file")
			SWEETError("Invalid argument for '--rexi-method='");



		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "EXP:" << std::endl;
		std::cout << " + exp_method: " << exp_method << std::endl;
		std::cout << " [EXP Taylor]" << std::endl;
		std::cout << " + taylor_num_expansions: " << taylor_num_expansions << std::endl;
		std::cout << "REXI generic parameters:" << std::endl;
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
		std::cout << " [EXP Direct]" << std::endl;
		std::cout << " + exp_direct_precompute_phin: " << exp_direct_precompute_phin << std::endl;
		std::cout << std::endl;
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
