/*
 * Author: Thibaut LUNET <thibaut.lunet@tuhh.de>
 */
#ifndef SRC_INCLUDE_SWEET_SHACKS_PARALLELIZATION_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_PARALLELIZATION_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>

#if SWEET_MPI
#   include <mpi.h>
#endif

namespace sweet
{

class ShackParallelization : 
    public ShackInterface
{
public:
    /**
     * Wether MPI is used or not
     */
    bool useMPI = false;

    /**
     * MPI size of the main communicator
     */
    int mpiSize = 1;

    /**
     * MPI rank of the current process
     */
    int mpiRank = 0;

    /**
     * Wether or not this process is root
     */
    bool isMPIRoot = true;

public:
    
    ShackParallelization() {
#if SWEET_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        useMPI = true;
        isMPIRoot = (mpiRank == 0);
#endif
    }

    virtual void printProgramArguments(const std::string &i_prefix = "") {}

    virtual bool processProgramArguments(ProgramArguments &i_pa) {
        return true;
    }

    virtual void printShack(
        const std::string& i_prefix = ""
    )
    {
        std::cout << i_prefix << std::endl;
        std::cout << i_prefix << "PARALLELIZATION:" << std::endl;
        std::cout << i_prefix << " + useMPI: " << useMPI << std::endl;
        std::cout << i_prefix << " + mpiSize: " << mpiSize << std::endl;
        std::cout << i_prefix << " + mpiRank: " << mpiRank << std::endl;

        std::cout << i_prefix << std::endl;
    }

};

} // namespace sweet

#endif
