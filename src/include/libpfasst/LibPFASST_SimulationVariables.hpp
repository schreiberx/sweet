#ifndef SRC_LIBPFASST_SIMULATION_VARIABLES_HPP_
#define SRC_LIBPFASST_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>
#include <string>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "../sweet/shacks/ShackInterface.hpp"

struct LibPFASST_SimulationVariables	:
	public sweet::ClassDictionaryInterface
{

    /**
     * Number of space-time levels in the PFASST algorithms
     */
    int nlevels = 1;

    /**
     * Number of (ML)SDC iterations
     */
    int niters = 8;

    /** 
     * Number of sweeps per level
     */
    std::string nsweeps_str = "";
    std::vector<int> nsweeps;

    /**
     * Number of SDC temporal nodes on the finest level
     */
    int nnodes = 5;

    /**
     * Type of SDC nodes
     */
    std::string nodes_type = "SDC_GAUSS_LOBATTO";

    /**
     * Coarsening multiplier for the spatial coarsening
     */
    double coarsening_multiplier = 0.5;


    /**
     * Use the RK stepper instead of the spectral deferred correction sweeps
     */
    bool use_rk_stepper = false;

    /**
     * hyperviscosity values
     */
    std::string hyperviscosity_4_str = "";
    std::string hyperviscosity_6_str = "";
    std::string hyperviscosity_2_str = "";
    std::string hyperviscosity_8_str = "";
    std::vector<double> hyperviscosity_2;
    std::vector<double> hyperviscosity_4;
    std::vector<double> hyperviscosity_6;
    std::vector<double> hyperviscosity_8;
    std::unordered_map<std::string, bool> hyperviscosity_on_field = 
	std::unordered_map<std::string, bool>({{"phi_pert", true},
										   {"div", true},
										   {"vrt", true}});

private:
    template <typename T>
    std::string _print_vector(const std::vector<T> & input_vector)
    {
        std::stringstream result;
        std::copy(input_vector.begin(), input_vector.end(), std::ostream_iterator<T>(result, "\t"));
        return result.str();
    }

    std::string _print_hyperviscosity_fields()
    {
        std::string hv_fields{};
        for (const auto & it : hyperviscosity_on_field)
        {
            if (!it.second) continue;
            hv_fields += it.first + "  ";
        }
        return hv_fields;
    }

    void _set_hyperviscosity_fields(const std::string & hv_fields)
    {
        if (hv_fields == "all")
        {
            return;
        }
        for (auto & it : hyperviscosity_on_field)
        {
            it.second = false;
        }
        if (hv_fields == "none")
        {
            return;
        }
        std::vector<std::string> substrings = StringSplit::split(hv_fields, ",");
        // first, set all entries to false
        // now, set entries true that are contained in the partial strings.
        for (const auto & substring : substrings)
        {
            if (hyperviscosity_on_field.find(substring) == hyperviscosity_on_field.end())
            {
                SWEETError("'" + substring + "' not a known field for hyperviscosity!");
            }
            hyperviscosity_on_field.at(substring) = true;
        }
    }
        
public:
    void outputConfig()
    {
    	printClass();
    }

    void printOptions()
    {
    	printProgramArguments();
    }

    void setup_longOptionList(
            struct option io_long_options[],		///< string and meta information for long options
            int &io_next_free_program_option,	///< number of free options, has to be increased for each new option
            int i_max_options					///< maximum number of options
    )
    {
        io_long_options[io_next_free_program_option] = {"libpfasst-nlevels", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-niters", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-nsweeps", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-nnodes", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-nodes-type", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-coarsening-multiplier", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-use-rk-stepper", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-u2", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-u4", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-u6", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-u8", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;

        io_long_options[io_next_free_program_option] = {"libpfasst-u-fields", required_argument, 0, 256+io_next_free_program_option};
        io_next_free_program_option++;
    }


    /**
     * @brief finalize hyperviscosity values using loaded LibPFASST parameters
     * 
     */
    void postprocess_hyperviscosity()
    {
        hyperviscosity_2 = std::vector<double>(nlevels, 0);
        hyperviscosity_4 = std::vector<double>(nlevels, 0);
        hyperviscosity_6 = std::vector<double>(nlevels, 0);
        hyperviscosity_8 = std::vector<double>(nlevels, 0);
        const std::string delimiter = ",";
        if (hyperviscosity_2_str != "")
        {
            StringSplit::split_n_doubles(hyperviscosity_2_str, hyperviscosity_2, delimiter);
        }
        if (hyperviscosity_4_str != "")
        {
            StringSplit::split_n_doubles(hyperviscosity_4_str, hyperviscosity_4, delimiter);
        }
        if (hyperviscosity_6_str != "")
        {
            StringSplit::split_n_doubles(hyperviscosity_6_str, hyperviscosity_6, delimiter);
        }
        if (hyperviscosity_8_str != "")
        {
            StringSplit::split_n_doubles(hyperviscosity_8_str, hyperviscosity_8, delimiter);
        }
    }

    /**
     * @brief finalize nsweeps values using loaded LibPFASST parameters
     * 
     */
    void postprocess_nsweeps()
    {
        nsweeps = std::vector<int>(nlevels, 1);
        const std::string delimiter = ",";
        if (nsweeps_str != "")
        {
            StringSplit::split_n_ints(nsweeps_str, nsweeps, delimiter);
        }
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
            case 0:	 nlevels                 = atoi(optarg); return -1;
            case 1:	 niters                  = atoi(optarg); return -1;
            case 2:  nsweeps_str             = optarg; return -1;
            case 3:	 nnodes                  = atoi(optarg); return -1;
            case 4:	 nodes_type              = optarg; 	     return -1;
            case 5:	 coarsening_multiplier   = atof(optarg); return -1;
            case 6:  use_rk_stepper          = atoi(optarg); return -1;
            case 7:  hyperviscosity_2_str    = optarg; return -1;
            case 8:  hyperviscosity_4_str    = optarg; return -1;
            case 9:  hyperviscosity_6_str    = optarg; return -1;
            case 10: hyperviscosity_8_str    = optarg; return -1;
            case 11: _set_hyperviscosity_fields(optarg); return -1;
        }
        return 12;
    }


	void printProgramArguments(const std::string& i_prefix = "")
	{
        std::cout << ""                                                                                                                        << std::endl;
        std::cout << "LibPFASST:"                                                                                                              << std::endl;
        std::cout << "\t--libpfasst-nlevels [int]                       LibPFASST parameter nlevels, default: 1"                               << std::endl;
        std::cout << "\t--libpfasst-niters [int]                        LibPFASST parameter niters, default: 8"                                << std::endl;
        std::cout << "\t--libpfasst-nsweeps [ints]                      LibPFASST parameter nsweeps, default: 1 on all levels"                 << std::endl;
        std::cout << "\t--libpfasst-nnodes [int]                        LibPFASST parameter nnodes, default: 5"                                << std::endl;
        std::cout << "\t--libpfasst-nodes-type [string]                 LibPFASST parameter nodes-type, default: SDC_GAUSS_LOBATTO"                << std::endl;
        std::cout << "\t--libpfasst-coarsening-multiplier [float]       LibPFASST parameter coarsening-multiplier, default: 0.5"               << std::endl;
        std::cout << "\t--libpfasst-use-rk-stepper [bool]               LibPFASST parameter use the Runge-Kutta stepper, default: false" << std::endl;
        std::cout << "\t--libpfasst-u2 [floats]                         Hyperviscosity of order 2, default: 0 on all levels"             << std::endl;
        std::cout << "\t--libpfasst-u4 [floats]                         Hyperviscosity of order 4, default: 0 on all levels"             << std::endl;
        std::cout << "\t--libpfasst-u6 [floats]                         Hyperviscosity of order 6, default: 0 on all levels"             << std::endl;
        std::cout << "\t--libpfasst-u8 [floats]                         Hyperviscosity of order 8, default: 0 on all levels"             << std::endl;
        std::cout << "\t--libpfasst-u-fields [string]                   Set fields for hyperviscosity, default: all"                     << std::endl;
        std::cout << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--libpfasst-nlevels", nlevels);
		i_pa.getArgumentValueByKey("--libpfasst-niters", niters);
		i_pa.getArgumentValueByKey("--libpfasst-nsweeps", nsweeps_str);
		i_pa.getArgumentValueByKey("--libpfasst-nnodes", nnodes);
		i_pa.getArgumentValueByKey("--libpfasst-nodes-type", nodes_type);
		i_pa.getArgumentValueByKey("--libpfasst-coarsening-multiplier", coarsening_multiplier);
		i_pa.getArgumentValueByKey("--libpfasst-use-rk-stepper", use_rk_stepper);
		i_pa.getArgumentValueByKey("--libpfasst-u2", hyperviscosity_2_str);
		i_pa.getArgumentValueByKey("--libpfasst-u4", hyperviscosity_4_str);
		i_pa.getArgumentValueByKey("--libpfasst-u6", hyperviscosity_6_str);
		i_pa.getArgumentValueByKey("--libpfasst-u8", hyperviscosity_8_str);

		std::string tmp;
		i_pa.getArgumentValueByKey("--libpfasst-u-fields", tmp);
		_set_hyperviscosity_fields(tmp);

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
        std::cout << std::endl;
        std::cout << "LibPFASST:" << std::endl;
        std::cout << " + nlevels: "                 << nlevels                 << std::endl;
        std::cout << " + niters: "                  << niters                  << std::endl;
        std::cout << " + nsweeps: "                 << _print_vector<int>(nsweeps)       << std::endl;
        std::cout << " + nnodes: "                  << nnodes                  << std::endl;
        std::cout << " + nodes_type: "              << nodes_type              << std::endl;
        std::cout << " + coarsening_multiplier: "   << coarsening_multiplier   << std::endl;
        std::cout << " + use_rk_stepper: "          << use_rk_stepper          << std::endl;
        std::cout << " + hyperviscosity order 2 [from coarse to fine]:  " << _print_vector<double>(hyperviscosity_2) << std::endl;
        std::cout << " + hyperviscosity order 4 [from coarse to fine]:  " << _print_vector<double>(hyperviscosity_4) << std::endl;
        std::cout << " + hyperviscosity order 6 [from coarse to fine]:  " << _print_vector<double>(hyperviscosity_6) << std::endl;
        std::cout << " + hyperviscosity order 8 [from coarse to fine]:  " << _print_vector<double>(hyperviscosity_8) << std::endl;
        std::cout << " + apply hyperviscosity on field(s):              " << _print_hyperviscosity_fields() << std::endl;
        std::cout << std::endl;
	}
};

#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
