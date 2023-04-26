/*
 * Author: JOAO STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_FILE_OUTPUT_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_FILE_OUTPUT_HPP_

#include <fstream>

// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include "ShackODEScalar.hpp"


class ODEScalar_FileOutput
{
public:
	sweet::ErrorBase error;

	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackODEScalar *shackODEScalar;

	void setup(
			sweet::ShackIOData *i_shackIOData,
			sweet::ShackTimestepControl *i_shackTimestepControl,
			ShackODEScalar *i_shackODEScalar
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		shackODEScalar = i_shackODEScalar;
	}

	void clear()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackODEScalar = nullptr;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const double u,
			const char* i_name		///< name of output variable
	)
	{
		char buffer[1024];

		const char* filename_template = shackIOData->output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << u;
		file.close();

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output(
			double u
	)
	{
		if (shackIOData->output_file_name.length() == 0)
			return;

		std::cout << "Writing output files at simulation time: " << shackTimestepControl->current_simulation_time << " secs" << std::endl;

		if (shackIOData->output_file_mode == "csv")
		{
			std::string output_filename;

			output_filename = write_file_csv(u, "u");
			output_reference_filenames += ";"+output_filename;
		}
		else if (shackIOData->output_file_mode == "bin")
		{
			SWEETError("Bin output not available for ODEScalar");
		}
		else
		{
			SWEETError("Unknown output file mode '"+shackIOData->output_file_mode+"'");
		}
	}
};

#endif
