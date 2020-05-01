#ifndef SRC_LIBPFASST_SIMULATION_VARIABLES_HPP_
#define SRC_LIBPFASST_SIMULATION_VARIABLES_HPP_

#include <unistd.h>
#include <getopt.h>
#include <string>

struct LibPFASST_SimulationVariables
{
  
  /**
   * Number of space-time levels in the PFASST algorithms
   */
  int nlevels = 2;

  /**
   * Number of (ML)SDC iterations
   */
  int niters = 8;

  /** 
   * Number of sweeps on the coarse level
   */
  int nsweeps_coarse = 1;

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
   * Use rexi as a linear solve in f2comp
   */ 
  bool use_rexi = false;

  /**
   * Treat the Coriolis force implicitly
   */
  bool implicit_coriolis_force = true;

  /**
   * Use the RK stepper instead of the spectral deferred correction sweeps
   */
  bool use_rk_stepper = false;

  void outputConfig()
  {
    std::cout << std::endl;
    std::cout << "LibPFASST:" << std::endl;
    std::cout << " + nlevels: "                 << nlevels                 << std::endl;
    std::cout << " + niters: "                  << niters                  << std::endl;
    std::cout << " + nsweeps_coarse: "          << nsweeps_coarse          << std::endl;
    std::cout << " + nnodes: "                  << nnodes                  << std::endl;
    std::cout << " + nodes_type: "              << nodes_type              << std::endl;
    std::cout << " + coarsening_multiplier: "   << coarsening_multiplier   << std::endl;
    std::cout << " + use_rexi: "                << use_rexi                << std::endl;
    std::cout << " + implicit_coriolis_force: " << implicit_coriolis_force << std::endl;
    std::cout << " + use_rk_stepper: "          << use_rk_stepper          << std::endl;
    std::cout                                                              << std::endl;
  }

  void printOptions()
  {
    std::cout << ""                                                                                                                      << std::endl;
    std::cout << "LibPFASST:"                                                                                                            << std::endl;
    std::cout << "	--libpfasst-nlevels [int]			LibPFASST parameter nlevels, default: 2"                         << std::endl;
    std::cout << "	--libpfasst-niters [int]                        LibPFASST parameter niters, default: 8"                          << std::endl;
    std::cout << "	--libpfasst-nsweeps_coarse [int]                LibPFASST parameter nsweeps_coarse, default: 1"                  << std::endl;
    std::cout << "	--libpfasst-nnodes [int]			LibPFASST parameter nnodes, default: 5"                          << std::endl;
    std::cout << "	--libpfasst-nodes_type [string]			LibPFASST parameter nodes_type, default: SDC_GAUSS_LOBATTO"      << std::endl;
    std::cout << "	--libpfasst-coarsening_multiplier [float]	LibPFASST parameter coarsening_multiplier, default: 0.5"         << std::endl;
    std::cout << "	--libpfasst-use_rexi [bool]	                LibPFASST parameter use_rexi, default: false"                    << std::endl;
    std::cout << "	--libpfasst-implicit_coriolis_force [bool]      LibPFASST parameter implicit_coriolis_force, default: true"      << std::endl;
    std::cout << "      --libpfasst-use_rk_stepper [bool]               LibPFASST parameter use the Runge-Kutta stepper, default: false" << std::endl;
    std::cout << ""                                                                                                                      << std::endl;
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

    io_long_options[io_next_free_program_option] = {"libpfasst-nsweeps_coarse", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;
    
    io_long_options[io_next_free_program_option] = {"libpfasst-nnodes", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;
    
    io_long_options[io_next_free_program_option] = {"libpfasst-nodes_type", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;
    
    io_long_options[io_next_free_program_option] = {"libpfasst-coarsening_multiplier", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;

    io_long_options[io_next_free_program_option] = {"libpfasst-use_rexi", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;

    io_long_options[io_next_free_program_option] = {"libpfasst-implicit_coriolis_force", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;

    io_long_options[io_next_free_program_option] = {"libpfasst-use_rk_stepper", required_argument, 0, 256+io_next_free_program_option};
    io_next_free_program_option++;

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
      case 0:	nlevels                 = atoi(optarg);	return -1;
      case 1:	niters                  = atoi(optarg);	return -1;
      case 2:   nsweeps_coarse          = atoi(optarg); return -1;
      case 3:	nnodes                  = atoi(optarg);	return -1;
      case 4:	nodes_type              = optarg; 	return -1;
      case 5:	coarsening_multiplier   = atof(optarg);	return -1;
      case 6:   use_rexi                = atoi(optarg); return -1;
      case 7:   implicit_coriolis_force = atoi(optarg); return -1;
      case 8:   use_rk_stepper          = atoi(optarg); return -1;
      }
    
    return 9;
  }
  

};

#endif /* SRC_SIMULATION_VARIABLES_HPP_ */
