/*
 * RexiNGCoefficients.hpp
 *
 *  Created on: 4 Jul 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_REXI_REXIFILECOEFFICIENTS_HPP_
#define SRC_INCLUDE_REXI_REXIFILECOEFFICIENTS_HPP_

#include <string>
#include <fstream>
#include <vector>


template <typename T>
class REXI_File_Coefficients
{
public:
	std::string function_name;		///< "phi0"; "phi1"; "phi2"; etc.

	int N;							///< Number of approximation poles

	T max_error;
	T max_error_double_precision;

	T test_min;
	T test_max;

	T basis_function_scaling;
	T basis_function_spacing;
	T basis_function_rat_shift;

	std::vector< std::complex<T> > weights_cplx;

	bool reduce_to_half = true;

	std::string filename;

	static
	constexpr
	T None()
	{
		return std::numeric_limits<T>::quiet_NaN();
	}

	static
	constexpr
	bool isNone(T i_value)
	{
		return std::isnan(i_value);
	}


public:
	REXI_File_Coefficients()
	{
		reset();
	}



private:
	void reset()
	{
		function_name = "";
		N = 0;						///< Number of approximation poles
		max_error = std::numeric_limits<T>::quiet_NaN();
		max_error_double_precision = std::numeric_limits<T>::quiet_NaN();
		test_min = std::numeric_limits<T>::quiet_NaN();
		test_max = std::numeric_limits<T>::quiet_NaN();
		basis_function_scaling = std::numeric_limits<T>::quiet_NaN();
		basis_function_spacing = std::numeric_limits<T>::quiet_NaN();
		basis_function_rat_shift = std::numeric_limits<T>::quiet_NaN();

		weights_cplx.resize(0);

		filename = "";
	}


public:
	void output()
	{
		std::cout << "N: " << N << std::endl;
		std::cout << "function_name: " << function_name << std::endl;
		std::cout << "max_error: " << max_error << std::endl;
		std::cout << "max_error_double_precision: " << max_error_double_precision << std::endl;
		std::cout << "test_min: " << test_min << std::endl;
		std::cout << "test_max: " << test_max << std::endl;
		std::cout << "basis_function_scaling: " << basis_function_scaling << std::endl;
		std::cout << "basis_function_spacing: " << basis_function_spacing << std::endl;
		std::cout << "basis_function_rat_shift: " << basis_function_rat_shift << std::endl;
	}

	void outputWeights()
	{
//		std::cout << "weights:" << std::endl;
		for (int i = 0; i < N; i++)
			std::cout << "faf_weight[" << i << "] = " << weights_cplx[i] << std::endl;
	}


public:
	bool load_from_file(std::string &i_filename)
	{
		filename = i_filename;

		std::ifstream infile(i_filename);

		if (!infile.is_open())
			FatalError(std::string("Unable to open file ")+i_filename);

		std::string line;

		while (std::getline(infile, line))
		{
			{
				std::string match = "# max_error ";
				if (match == line.substr(0, match.length()))
				{
					max_error = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# max_error_double_precision ";
				if (match == line.substr(0, match.length()))
				{
					max_error_double_precision = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# N ";
				if (match == line.substr(0, match.length()))
				{
					N = std::atoi(line.substr(match.length()).c_str());
					continue;
				}
			}

			{
				std::string match = "# function_name ";
				if (match == line.substr(0, match.length()))
				{
					function_name = line.substr(match.length());
					continue;
				}
			}

			{
				std::string match = "# basis_function_spacing ";
				if (match == line.substr(0, match.length()))
				{
					basis_function_spacing = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# basis_function_scaling ";
				if (match == line.substr(0, match.length()))
				{
					basis_function_scaling = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# basis_function_rat_shift ";
				if (match == line.substr(0, match.length()))
				{
					basis_function_rat_shift = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# test_min ";
				if (match == line.substr(0, match.length()))
				{
					test_min = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# test_max ";
				if (match == line.substr(0, match.length()))
				{
					test_max = DQStuff::fromString<T>(line.substr(match.length()));
					continue;
				}
			}

			{
				std::string match = "# ";
				if (match == line.substr(0, match.length()))
					continue;
			}

			std::size_t pos = line.find('\t');

			if (pos == std::string::npos)
				FatalError(std::string("No tab found in file ")+i_filename);

			std::complex<T> val(
					DQStuff::fromString<T>(line.substr(0, pos)),
					DQStuff::fromString<T>(line.substr(pos+1))
			);

			weights_cplx.push_back(val);
		}

		return true;
	}
};



#endif /* SRC_INCLUDE_REXI_REXIFILECOEFFICIENTS_HPP_ */
