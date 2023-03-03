/*
 * ShackIOData.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKIODATA_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKIODATA_HPP_


#include <getopt.h>
#include <iomanip>
#include <string>
#include <iostream>
#include <limits>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>


namespace sweet
{

/**
 * Input and output data
 */
class ShackIOData	:
		public ShackInterface
{
public:
	/// filenames of input data for setup (this has to be setup by each application individually)
	std::vector<std::string> initial_condition_data_filenames;

	/// use "BINARY;filename1;filename2" to specify that the binary files should be read in binary format
	bool initial_condition_input_data_binary = false;


	/// prefix of filename for outputConfig of data
	std::string output_file_name = "X";

	/// output mode of variables
	std::string output_file_mode = "default";

	/// prefix of filename for outputConfig of data
	double output_each_sim_seconds = -1;

	/// Simulation seconds for next outputConfig
	double output_next_sim_seconds = 0;

	/// time scaling for outputConfig
	/// e.g. use scaling by 1.0/(60*60) to output days instead of seconds
	double output_time_scale = 1.0;

	/// precision for floating point outputConfig to std::cout and std::endl
	int output_floating_point_precision = std::numeric_limits<double>::digits10 + 1;


	/// set verbosity of simulation
	int verbosity = 0;

	/// activate GUI mode?
	bool gui_enabled = (SWEET_GUI == 0 ? false : true);



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

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << "IOData:" << std::endl;
		std::cout << "	--output-file-name [string]		String specifying the name of the output file" << std::endl;
		std::cout << "	--output-file-mode [string]		Format of output file, default: default" << std::endl;
		std::cout << "	-v [int]			verbosity level" << std::endl;
		std::cout << "	-G [0/1]			graphical user interface" << std::endl;

		std::cout << "" << std::endl;
	}

	bool processProgramArguments(ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--output-file-mode", output_file_mode);
		i_pa.getArgumentValueBy2Keys("--output-file-name", "-O", output_file_name);
		i_pa.getArgumentValueByKey("-d", output_floating_point_precision);
		i_pa.getArgumentValueByKey("-o", output_each_sim_seconds);

		std::string tmp;
		if (i_pa.getArgumentValueByKey("-i", tmp))
			setup_initial_condition_filenames(optarg);

		if (i_pa.error.exists())
			return error.forwardWithPositiveReturn(i_pa.error);

		if (output_file_mode == "default")
		{
			output_file_mode = "bin";

			if (output_file_name == "X")
				output_file_name = "output_%s_t%020.8f.sweet";
		}
		else
		{
			if (output_file_name == "X")
			{
				if (output_file_mode == "csv")
					output_file_name = "output_%s_t%020.8f.csv";
				else if (output_file_mode == "bin")
					output_file_name = "output_%s_t%020.8f.sweet";
				else if (output_file_mode == "csv_spec_evol")
					output_file_name = "output_%s_t%020.8f.txt";
				else
					return error.set("Unknown filemode '"+output_file_mode+"'");
			}
		}

		if (output_file_name == "-")
			output_file_name = "";

		if (output_floating_point_precision >= 0)
		{
			std::cout << std::setprecision(output_floating_point_precision);
			std::cerr << std::setprecision(output_floating_point_precision);
		}

		i_pa.getArgumentValueByKey("-G", gui_enabled);
		i_pa.getArgumentValueByKey("-v", verbosity);

		ERROR_FORWARD_WITH_RETURN_BOOLEAN(i_pa);
	}

	virtual void printShack(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "INPUT/OUTPUT:" << std::endl;
		for (std::size_t i = 0; i < initial_condition_data_filenames.size(); i++)
			std::cout << "    - filename: " << i << " " << initial_condition_data_filenames[i] << std::endl;
		std::cout << " + input_data_binary: " << initial_condition_input_data_binary << std::endl;
		std::cout << " + output_file_name: " << output_file_name << std::endl;
		std::cout << " + output_file_mode: " << output_file_mode << std::endl;
		std::cout << " + output_each_sim_seconds: " << output_each_sim_seconds << std::endl;
		std::cout << " + output_next_sim_seconds: " << output_next_sim_seconds << std::endl;
		std::cout << " + output_time_scale: " << output_time_scale << std::endl;
		std::cout << " + output_floating_point_precision: " << output_floating_point_precision << std::endl;

		std::cout << " + verbosity: " << verbosity << std::endl;
		std::cout << " + gui_enabled: " << gui_enabled << std::endl;

		std::cout << std::endl;
	}
};

}

#endif
