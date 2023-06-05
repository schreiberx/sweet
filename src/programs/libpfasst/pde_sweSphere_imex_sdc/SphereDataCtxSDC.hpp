#ifndef _SPHERE_DATA_CTX_SDC_HPP_
#define _SPHERE_DATA_CTX_SDC_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include "../../pde_sweSphere/PDESWESphere_Diagnostics.hpp"
#include <sweet/core/sphere/SphereOperators.hpp>
#include <vector>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "../interface/LevelSingleton.hpp"

#include "../../pde_sweSphere/time/PDESWESphereTS_lg_erk_lc_n_erk.hpp"
#include "../../pde_sweSphere/time/PDESWESphereTS_lg_irk.hpp"

#include "../ShackLibPFASST.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the LevelSingleton and the sweet::ShackDictionary object

class SphereDataCtxSDC
{

public:
    sweet::ShackDictionary *shackDict;

    sweet::ShackSphereDataOps *shackSphereDataOps;
    sweet::ShackIOData *shackIOData;
    sweet::ShackTimestepControl *shackTimestepControl;
    ShackPDESWESphere *shackPDESWESphere;
    ShackLibPFASST *shackLibPFASST;
    ShackPDESWESphereTimeDiscretization *shackTimeDisc;

    PDESWESphere_Diagnostics diagnostics;

    bool shackRegistration(
        sweet::ShackDictionary *io_shackDict)
    {
        shackDict = io_shackDict;

        shackSphereDataOps = shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        shackIOData = shackDict->getAutoRegistration<sweet::ShackIOData>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        shackPDESWESphere = shackDict->getAutoRegistration<ShackPDESWESphere>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        shackTimeDisc = shackDict->getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        shackTimestepControl = shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        shackLibPFASST = shackDict->getAutoRegistration<ShackLibPFASST>();
        ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(*shackDict);

        return true;
    }

    SphereDataCtxSDC()
    {
    }

    bool setup(sweet::ShackDictionary *i_shackDict, LevelSingleton *i_singleton, int *i_nnodes)
    {
        shackDict = i_shackDict;
        levelSingleton = i_singleton;

        int rank = 0;
        int nprocs = 0;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        if (!shackDict)
        {
            SWEETError("SphereDataCtx: shackDict pointer is NULL!");
        }

        if (!levelSingleton)
        {
            SWEETError("SphereDataCtx: levelSingleton pointer is NULL!");
        }

        // initialize the timesteppers
        timestepper_lg_erk_lc_n_erk = new PDESWESphereTS_lg_erk_lc_n_erk;
        timestepper_lg_irk = new PDESWESphereTS_lg_irk;

        timestepper_lg_erk_lc_n_erk->shackRegistration(shackDict);
        timestepper_lg_irk->shackRegistration(shackDict);

        int timestepping_order = 1;
        // Strang splitting version that should be used:
        int version_id = 1;

        // setup lg_irk timestepper
        timestepper_lg_erk_lc_n_erk->setup(&(*levelSingleton).ops, timestepping_order, version_id);
        timestepper_lg_irk->setup(&(*levelSingleton).ops, timestepping_order, shackTimestepControl->current_timestepSize);

        // initialize the residuals
        residuals.resize(nprocs, std::vector<double>(0, 0.));

        // initialize the diagnostics object
        diagnostics.setup(
            &(*i_singleton).ops,
            shackPDESWESphere,
            0);

        return true;
    }

    // Destructor
    ~SphereDataCtxSDC()
    {
        delete timestepper_lg_erk_lc_n_erk;
        delete timestepper_lg_irk;
    }

    // Getter for the sphere data configuration
    sweet::SphereData_Config *get_sphere_data_config() const
    {
        return &(levelSingleton->sphereDataConfig);
    }

    // Getter for the sphere data configuration
    PDESWESphere_BenchmarksCombined *get_swe_benchmark() const
    {
        return &(levelSingleton->benchmarks);
    }

    // Getter for the sphere data operators
    sweet::SphereOperators *get_sphere_operators() const
    {
        return &(levelSingleton->ops);
    }

#if 0
    // Getter for the sphere data operators with no dealiasing
    sweet::SphereOperators *get_sphere_operators_nodealiasing() const
    {
        return &(levelSingleton->opNoDealiasing);
    }
#endif

    // Getter for the sphere diagnostics
    PDESWESphere_Diagnostics *get_sphere_diagnostics()
    {
        return &diagnostics;
    }

    // Getter for the explicit timestepper
    PDESWESphereTS_lg_erk_lc_n_erk *get_lg_erk_lc_n_erk_timestepper() const
    {
        return timestepper_lg_erk_lc_n_erk;
    }

    // Getter for the implicit timestepper
    PDESWESphereTS_lg_irk *get_lg_irk_timestepper() const
    {
        return timestepper_lg_irk;
    }

    // Getter for the simulationVariables object
    sweet::ShackDictionary *get_simulation_variables() const
    {
        return shackDict;
    }

    // Getter for the number of levels
    int get_number_of_levels() const
    {
        return 1;
    }

    // Save the physical invariants
    void save_physical_invariants(
        int i_niter)
    {
        time.push_back(shackTimestepControl->current_timestepSize * i_niter);
        mass.push_back(diagnostics.total_mass);
        energy.push_back(diagnostics.total_energy);
        potentialEnstrophy.push_back(diagnostics.total_potential_enstrophy);
    }

    // Getters for the time and invariants vectors
    const std::vector<double> &get_time() const { return time; }
    const std::vector<double> &get_mass() const { return mass; }
    const std::vector<double> &get_energy() const { return energy; }
    const std::vector<double> &get_potential_enstrophy() const { return potentialEnstrophy; }

    // Getters for the residuals
    const std::vector<std::vector<double>> &get_residuals() const { return residuals; }
    std::vector<std::vector<double>> &get_residuals() { return residuals; }

protected:
    // Pointer to the LevelSingleton object
    LevelSingleton *levelSingleton;

    // Pointer to the lg_erk_lc_n_erk timestepper (used for ceval_f1, ceval_f2)
    PDESWESphereTS_lg_erk_lc_n_erk *timestepper_lg_erk_lc_n_erk;

    // Pointer to the lg_irk timestepper (used for ccomp_f2)
    PDESWESphereTS_lg_irk *timestepper_lg_irk;

    // Saved Residuals for each processor
    std::vector<std::vector<double>> residuals;

    // Some constructors and operator= are disabled
    SphereDataCtxSDC(const SphereDataCtxSDC &);
    SphereDataCtxSDC &operator=(const SphereDataCtxSDC &);

    // Vectors used for plotting
    std::vector<double> time;
    std::vector<double> mass;
    std::vector<double> energy;
    std::vector<double> potentialEnstrophy;
};

#endif // _SPHERE_DATA_CTX_SDC_HPP_
