/*
 * IOParameters.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: martin
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_VARIABLES_IODATAPARAMETERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_VARIABLES_IODATAPARAMETERS_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <limits>

#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include <sweet/class_dict/ClassDictionaryInterface.hpp>

/**
 * Parameters related to input and output data
 */
struct IODataParameters	:
	public ClassDictionaryInterface
{
public:
	sweet::ErrorBase error;

	sweet::ErrorBase& getError()
	{
		return error;
	}

	/**
	 * Filenames of input data for setup (this has to be setup by each application individually)
	 */
	std::vector<std::string> initial_condition_data_filenames;

	/*
	 * Use "BINARY;filename1;filename2" to specify that the binary files should be read in binary format
	 */
	bool initial_condition_input_data_binary = false;


	/*
	 * Prefix of filename for outputConfig of data
	 */
	std::string output_file_name = "X";

	/*
	 * Output mode of variables
	 */
	std::string output_file_mode = "default";

	/*
	 * Prefix of filename for outputConfig of data
	 */
	double output_each_sim_seconds = -1;

	/*
	 *
	 */
	/// Simulation seconds for next outputConfig
	double output_next_sim_seconds = 0;

	/*
	 * Time scaling for outputConfig
	 *
	 * E.g. use scaling by 1.0/(60*60) to output days instead of seconds
	 */
	double output_time_scale = 1.0;

	/*
	 * Precision for floating point outputConfig to std::cout and std::endl
	 */
	int output_floating_point_precision = std::numeric_limits<double>::digits10 + 1;


	void setup_initial_condition_filenames(
			const std::string i_string
	)
	{
		std::size_t last_pos = 0;
		for (std::size_t pos = 0; i_string[pos] != '\0'; pos++)
		{
			if (i_string[pos] != ';')
				continue;

			initial_condition_data_filenames.push_back(i_string.substr(last_pos, pos-last_pos));
			last_pos = pos+1;
		}

		initial_condition_data_filenames.push_back(i_string.substr(last_pos));

		if (initial_condition_data_filenames.size() > 0)
		{
			if (initial_condition_data_filenames[0] == "BINARY")
			{
				initial_condition_data_filenames.erase(initial_condition_data_filenames.begin());
				initial_condition_input_data_binary = true;
			}
		}
	}

	void setup_longOptionsList(
			struct option *long_options,
			int &next_free_program_option
	)
	{
	}

	int setup_longOptionValue(
			int i_option_index,		///< Index relative to the parameters setup in this class only, starts with 0
			const char *i_value		///< Value in string format
	)
	{
		switch(i_option_index)
		{
		case 0:
			output_file_name = i_value;
			return -1;

		case 1:
			output_file_mode = i_value;
			return -1;
		}

		return 2;
	}


	void outputVariables(std::string i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << i_prefix << "INPUT/OUTPUT:" << std::endl;
		for (std::size_t i = 0; i < initial_condition_data_filenames.size(); i++)
			std::cout << i_prefix << "   + filename " << i << " " << initial_condition_data_filenames[i] << std::endl;
		std::cout << i_prefix << " + input_data_binary: " << initial_condition_input_data_binary << std::endl;
		std::cout << i_prefix << " + output_file_name " << output_file_name << std::endl;
		std::cout << i_prefix << " + output_file_mode " << output_file_mode << std::endl;
		std::cout << i_prefix << " + output_each_sim_seconds: " << output_each_sim_seconds << std::endl;
		std::cout << i_prefix << " + output_next_sim_seconds: " << output_next_sim_seconds << std::endl;
		std::cout << i_prefix << " + output_time_scale: " << output_time_scale << std::endl;
		std::cout << i_prefix << " + output_floating_point_precision: " << output_floating_point_precision << std::endl;
		std::cout << i_prefix << std::endl;
	}


	void outputProgramArguments(std::string i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << i_prefix << "IOData:" << std::endl;
		std::cout << i_prefix << "	--output-file-name [string]		String specifying the name of the output file" << std::endl;
		std::cout << i_prefix << "	--output-file-mode [string]		Format of output file, default: default" << std::endl;
	}


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (!i_pa.getArgumentValueByKey("output-file-name", output_file_name))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		if (!i_pa.getArgumentValueByKey("output-file-mode", output_file_mode))
		{
			if (error.errorForwardFrom(i_pa.error))
				return false;
		}

		return true;
	}

};




#endif /* SRC_PROGRAMS_SWE_SPHERE_VARIABLES_IODATAPARAMETERS_HPP_ */
