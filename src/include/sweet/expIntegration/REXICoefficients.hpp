/*
 *  Created on: Dec 9, 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_REXI_REXI_COEFFICIENTS_HPP__
#define SRC_INCLUDE_REXI_REXI_COEFFICIENTS_HPP__

#include <vector>
#include <complex>
#include <sweet/libmath/DQStuff.hpp>
#include <fstream>

namespace sweet
{

template <typename T = double>
class REXICoefficients
{
public:
	typedef std::complex<T> TComplex;

	std::vector<TComplex> alphas;
	std::vector<TComplex> betas;
	TComplex gamma;

	std::string filename;
	std::string function_name;


	/**
	 * Constructor
	 */
	REXICoefficients()	:
		gamma(0)
	{
	}

	/**
	 * Load binary double precision complex value
	 */
	static
	void load_binary_to_complex(
		std::ifstream &i_file,
		TComplex &o_output
	)
	{
		double val[2];
		i_file.read((char*)val, sizeof(double)*2);

		o_output.real(val[0]);
		o_output.imag(val[1]);
	}


	/*
	 * Load array of double precision complex values
	 */
	static
	void load_binary_to_vector(
		std::ifstream &i_file,
		std::size_t i_N,
		std::vector<TComplex> &o_output
	)
	{
		o_output.resize(i_N);
		for (std::size_t i = 0; i < i_N; i++)
		{
			double val[2];
			i_file.read((char*)val, sizeof(double)*2);

			o_output[i].real(val[0]);
			o_output[i].imag(val[1]);
		}
	}



public:
	/*
	 * File format
	 *
	 * # N 2				<- load 2 coefficients
	 * # function phi0		<- string describing function name
	 * # binary [0/1]		<- values are stored in binary format, all double precision and complex-valued
	 * # gamma				<- gamma coefficient, must be real-valued only
	 * 1.0 0.1
	 * # alphas				<- alphas start here
	 * 1 -3.2
	 * 2.3 3
	 * # betas				<- betas start here
	 * 5.1 -1
	 * 3 1.3
	 */
	bool load_from_file(
			const std::string &i_filename
	)
	{
		filename = i_filename;
		std::ifstream infile(i_filename, std::ios::in | std::ios::binary);

		if (!infile.is_open())
			SWEETError(std::string("Unable to open file ")+i_filename);

		bool binary = false;

		std::string line;
		int abg_mode = -1;	// 0: alpha, 1: beta, 2: gamma

		int N = -1;

		int line_nr = 0;
		while (std::getline(infile, line))
		{
			line_nr++;

			{
				std::string match = "# N ";
				if (match == line.substr(0, match.length()))
				{
					N = std::atoi(line.substr(match.length()).c_str());
					continue;
				}
			}

			{
				std::string match = "# binary ";
				if (match == line.substr(0, match.length()))
				{
					binary = std::atoi(line.substr(match.length()).c_str());
					continue;
				}
			}

			{
				std::string match = "# function_name ";
				if (match == line.substr(0, match.length()))
				{
					function_name = line.substr(match.length()).c_str();
					continue;
				}
			}

			{
				std::string match = "# alphas";
				if (match == line.substr(0, match.length()))
				{
					if (binary)
					{
						load_binary_to_vector(infile, N, alphas);
						continue;
					}

					abg_mode = 0;
					continue;
				}
			}

			{
				std::string match = "# betas";
				if (match == line.substr(0, match.length()))
				{
					if (binary)
					{
						load_binary_to_vector(infile, N, betas);
						continue;
					}

					abg_mode = 1;
					continue;
				}
			}

			{
				std::string match = "# gamma";
				if (match == line.substr(0, match.length()))
				{
					if (binary)
					{
						load_binary_to_complex(infile, gamma);
						continue;
					}

					abg_mode = 2;
					continue;
				}
			}

			bool real_valued_only = false;

			// search for tab
			std::size_t pos = line.find('\t');

			if (pos == std::string::npos)
			{
				// if tab wasn't found, search for whitespace
				pos = line.find(' ');

				if (pos == std::string::npos)
					real_valued_only = true;
			}

			T val_real, val_imag;
			if (real_valued_only)
			{
				val_real = DQStuff::fromString<T>(line);
				val_imag = 0;
			}
			else
			{
				val_real = DQStuff::fromString<T>(line.substr(0, pos));
				val_imag = DQStuff::fromString<T>(line.substr(pos+1));
			}

			std::complex<T> val(val_real, val_imag);

			if (abg_mode == 0)
			{
				alphas.push_back(val);
				//std::cout << "alphas: " << val << std::endl;
				continue;
			}
			if (abg_mode == 1)
			{
				betas.push_back(val);
				//std::cout << "betas: " << val << std::endl;
				continue;
			}

			if (abg_mode == 2)
			{
				if (val_imag != 0)
					SWEETError("Gamma must be real-only!");

				gamma = val_real;
				//std::cout << "gamma: " << val << std::endl;
				continue;
			}

			std::cout << "Error in line " << line_nr << std::endl;
			std::cout << "Line content: " << line << std::endl;
			SWEETError("Line contains bogus data");
		}

		if ((int)alphas.size() != N || (int)betas.size() != N)
		{
			std::cout << "alphas.size: " << alphas.size() << std::endl;
			std::cout << "betas.size: " << betas.size() << std::endl;
			std::cout << "N: " << N << std::endl;
			SWEETError("Size doesn't match!");
		}

		return true;
	}
};

}

#endif
