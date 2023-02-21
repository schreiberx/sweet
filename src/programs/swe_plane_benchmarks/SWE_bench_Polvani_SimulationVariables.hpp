/*
 * SWEPolvani_SimulationVariables.hpp
 *
 *  Created on: 14 Oct 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_SWE_POLVANI_SIMULATION_VARIABLES_HPP_
#define SRC_SWE_POLVANI_SIMULATION_VARIABLES_HPP_


#include <getopt.h>
#include <sweet/SWEETError.hpp>

struct SWEPolvani_SimulationVariables	:
	public sweet::ClassDictionaryInterface
{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
	/**
	 * Rossby number
	 */
	double r = 0.01;

	/**
	 * Froude number
	 */
	double f = 0.04;
#endif


	void outputConfig()
	{
		printClass();
	}


	void outputProgParams()
	{
		printProgramArguments();
	}


	void setup_longOptionList(
			struct option io_long_options[],		///< string and meta information for long options
			int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
			int i_max_options					///< maximum number of options
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		io_long_options[io_next_free_program_option] = {"polvani-rossby", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;

		io_long_options[io_next_free_program_option] = {"polvani-froude", required_argument, 0, 256+io_next_free_program_option};
		io_next_free_program_option++;
#endif
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
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		switch(i_option_index)
		{
			case 0:	r = atof(optarg);	return -1;
			case 1:	f = atof(optarg);	return -1;
		}
#endif

		return 2;
	}



	void printProgramArguments(const std::string& i_prefix = "")
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		std::cout << "" << std::endl;
		std::cout << "Polvani benchmark settings (on the plane):" << std::endl;
		std::cout << "	--polvani-rossby [float]	Choose Rossby number, default:0" << std::endl;
		std::cout << "	--polvani-froude [float]	Choose Froude number, default:0" << std::endl;
		std::cout << "" << std::endl;
#endif
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		i_pa.getArgumentValueByKey("--polvani-rossby", r);
		i_pa.getArgumentValueByKey("--polvani-froude", f);
#endif

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		std::cout << std::endl;
		std::cout << "SWEPolvani:" << std::endl;
		std::cout << " + Rossby number: " << r << std::endl;
		std::cout << " + Froude number: " << f << std::endl;
		std::cout << std::endl;
#endif
	}
};



#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
