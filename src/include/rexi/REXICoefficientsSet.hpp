/*
 * REXICoefficientsSet.hpp
 *
 *  Created on: Dec 9, 2018
 *      Author: martin
 */

#ifndef SRC_INCLUDE_REXI_COEFFICIENTS_SET_HPP_
#define SRC_INCLUDE_REXI_COEFFICIENTS_SET_HPP_

#include <vector>
#include <complex>
#include "REXICoefficients.hpp"
#include <sweet/StringSplit.hpp>


template <typename T = double>
class REXICoefficientsSet
{
public:
	typedef std::complex<T> TComplex;

	std::vector< REXICoefficients<> > rexiCoefficientVector;


	/**
	 * Load REXI coefficients from filenames
	 */
	void setup_from_files(
			const std::string &i_rexi_filenames
	)
	{
		rexiCoefficientVector.clear();

		// do simply nothing if there's no file name
		if (i_rexi_filenames == "")
			return;

		std::vector<std::string> rexi_filenames = StringSplit::split(i_rexi_filenames, ",");

		for (auto iter = rexi_filenames.begin(); iter != rexi_filenames.end(); iter++)
		{
			std::string& rexi_filename = *iter;

			std::vector<std::string> split2 = StringSplit::split(rexi_filename, ":");

			if (split2.size() == 0)
			{
				SWEETError("Strange things....");
			}
			else if (split2.size() == 1)
			{
				REXICoefficients<T> rexiCoefficients;
				rexiCoefficients.load_from_file(split2[0]);
			}
			else if (split2.size() == 2)
			{
				REXICoefficients<T> rexiCoefficients;
				rexiCoefficients.load_from_file(split2[1]);

				if (rexiCoefficients.filename != "")
					if (rexiCoefficients.function_name != split2[0])
						SWEETError("Function name mismatch!");

				this->rexiCoefficientVector.push_back(rexiCoefficients);
			}
			else
			{
				SWEETError("Too many split variable names");
			}
		}
	}


	const REXICoefficients<>& find_by_function_name(
			const std::string &i_function_name
	)	const
	{
		for (auto iter = rexiCoefficientVector.begin(); iter != rexiCoefficientVector.end(); iter++)
		{
			if (iter->function_name == i_function_name)
				return *iter;
		}

		SWEETError("Not found");
	}
};


#endif /* SRC_INCLUDE_REXI_REXICOEFFICIENTS_HPP_ */
