#ifndef _SPHERE_DATA_CTX_SDC_HPP_
#define _SPHERE_DATA_CTX_SDC_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <vector>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "../libpfasst_interface/LevelSingleton.hpp"

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_irk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the LevelSingleton and the sweet::ShackDictionary object

class SphereDataCtxSDC {

public:

    // Constructor
    SphereDataCtxSDC(
        sweet::ShackDictionary *i_shackDict,
        LevelSingleton *i_singleton,
        int* i_nnodes
        ) 
        : shackDict(i_shackDict), levelSingleton(i_singleton)
    {
        int rank   = 0;
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

        // use first order integration in time for all pieces (only order supported)
        if (shackDict->disc.timestepping_order != -1)
        {
            std::cout << "WARNING: Supplying the timestepping order manually is not supported!" << std::endl;
        }
        int timestepping_order  = 1; 
        
        // initialize the timesteppers
        timestepper_lg_erk_lc_n_erk = new SWE_Sphere_TS_lg_erk_lc_n_erk(*shackDict, levelSingleton->op);
        timestepper_lg_irk = new SWE_Sphere_TS_lg_irk(*shackDict, levelSingleton->op);

        // setup lg_erk_lc_n_erk timestepper
        int num_off_diagonals = 1; // TODO why?
        timestepper_lg_erk_lc_n_erk->setup(timestepping_order, num_off_diagonals);
        // TODO what to do about the topography?

        // setup lg_irk timestepper
        timestepper_lg_irk->setup(timestepping_order, shackDict->timecontrol.current_timestep_size);	

        // initialize the residuals
        residuals.resize(nprocs,std::vector<double>(0,0.));

        // initialize the diagnostics object
        sphereDiagnostics = new SphereHelpers_Diagnostics(
                            &(levelSingleton->dataConfig),
                            *shackDict,
                            0
                            );
    }

    // Destructor
    ~SphereDataCtxSDC() 
    {
        delete timestepper_lg_erk_lc_n_erk;
        delete timestepper_lg_irk;
        delete sphereDiagnostics;
    }

    // Getter for the sphere data configuration
    SphereData_Config* get_sphere_data_config() const 
    {
        return &(levelSingleton->dataConfig);
    }

    // Getter for the sphere data configuration
    BenchmarksSphereSWE* get_swe_benchmark() const 
    {
        return &(levelSingleton->benchmarks);
    }

    // Getter for the sphere data operators
    SphereOperators* get_sphere_operators() const
    {
        return &(levelSingleton->op);
    }

    // Getter for the sphere data operators with no dealiasing
    SphereOperators* get_sphere_operators_nodealiasing() const
    {
        return &(levelSingleton->opNoDealiasing);
    }

    // Getter for the sphere diagnostics
    SphereHelpers_Diagnostics* get_sphere_diagnostics() 
    {
        return sphereDiagnostics;
    }

    // Getter for the explicit timestepper
    SWE_Sphere_TS_lg_erk_lc_n_erk* get_lg_erk_lc_n_erk_timestepper() const
    {
        return timestepper_lg_erk_lc_n_erk;
    }

    // Getter for the implicit timestepper
    SWE_Sphere_TS_lg_irk* get_lg_irk_timestepper() const
    {
        return timestepper_lg_irk;
    }

    // Getter for the simulationVariables object
    sweet::ShackDictionary* get_simulation_variables() const 
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
                    int i_niter
                    ) 
    {
        time.push_back(shackDict->timecontrol.current_timestep_size * i_niter);
        mass.push_back(shackDict->diag.total_mass);
        energy.push_back(shackDict->diag.total_energy);
        potentialEnstrophy.push_back(shackDict->diag.total_potential_enstrophy);
    }

    // Getters for the time and invariants vectors
    const std::vector<double>& get_time()                const { return time; }
    const std::vector<double>& get_mass()                const { return mass; }
    const std::vector<double>& get_energy()              const { return energy; }
    const std::vector<double>& get_potential_enstrophy() const { return potentialEnstrophy; }

    // Getters for the residuals
    const std::vector<std::vector<double>>& get_residuals() const { return residuals; }
    std::vector<std::vector<double>>&       get_residuals()       { return residuals; }
    
protected:

    // Pointer to the sweet::ShackDictionary object
    sweet::ShackDictionary *shackDict;

    // Pointer to the LevelSingleton object
    LevelSingleton *levelSingleton;

    // Pointer to the lg_erk_lc_n_erk timestepper (used for ceval_f1, ceval_f2)
    SWE_Sphere_TS_lg_erk_lc_n_erk* timestepper_lg_erk_lc_n_erk;

    // Pointer to the lg_irk timestepper (used for ccomp_f2)
    SWE_Sphere_TS_lg_irk* timestepper_lg_irk;

    // Saved Residuals for each processor
    std::vector<std::vector<double>> residuals;

    // Diagnostics (mass, energy, enstrophy)
    SphereHelpers_Diagnostics* sphereDiagnostics;

    // Some constructors and operator= are disabled
    SphereDataCtxSDC() {};
    SphereDataCtxSDC(const SphereDataCtxSDC&);
    SphereDataCtxSDC& operator=(const SphereDataCtxSDC&);
    
    // Vectors used for plotting
    std::vector<double> time;
    std::vector<double> mass;
    std::vector<double> energy;
    std::vector<double> potentialEnstrophy;
};

#endif // _SPHERE_DATA_CTX_SDC_HPP_
