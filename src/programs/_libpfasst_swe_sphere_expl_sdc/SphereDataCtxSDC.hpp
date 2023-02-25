#ifndef _SPHERE_DATA_CTX_SDC_HPP_
#define _SPHERE_DATA_CTX_SDC_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <vector>
#include <sweet/core/SimulationVariables.hpp>
#include "../libpfasst_interface/LevelSingleton.hpp"

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_erk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the LevelSingleton and the SimulationVariables object

class SphereDataCtxSDC {

public:

    // Constructor
    SphereDataCtxSDC(
        SimulationVariables *i_simVars,
        LevelSingleton *i_singleton,
        int* i_nnodes
        ) 
        : simVars(i_simVars), levelSingleton(i_singleton)
    {
        int rank   = 0;
        int nprocs = 0;
        
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        
        if (!simVars)
        {
            SWEETError("SphereDataCtx: simVars pointer is NULL!");
        }
        
        if (!levelSingleton)
        {
            SWEETError("SphereDataCtx: levelSingleton pointer is NULL!");
        }
        
        // initialize the ln_erk time stepper
        timestepper_ln_erk = new SWE_Sphere_TS_ln_erk(*simVars, levelSingleton->op);

        // this is never used but this makes clear that with niters=1,
        // we're actually just calling ERK1
        timestepper_ln_erk->setup(1);
    
        // initialize the residuals
        residuals.resize(nprocs,std::vector<double>(0,0.));

        // initialize the diagnostics object
        sphereDiagnostics = new SphereHelpers_Diagnostics(
                            &(levelSingleton->dataConfig),
                            *simVars,
                            0
                            );
    }

    // Destructor
    ~SphereDataCtxSDC() 
    {
        delete timestepper_ln_erk;
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

    // Getter for the sphere data configuration with no dealiasing
    SphereData_Config* get_sphere_data_config_nodealiasing() const 
    {
        return &(levelSingleton->dataConfigNoDealiasing);
    }

    // Getter for the sphere data operators
    SphereOperators_SphereData* get_sphere_operators() const
    {
        return &(levelSingleton->op);
    }

    // Getter for the sphere data operators with no dealiasing
    SphereOperators_SphereData* get_sphere_operators_nodealiasing() const
    {
        return &(levelSingleton->opNoDealiasing);
    }

    // Getter for the sphere diagnostics
    SphereHelpers_Diagnostics* get_sphere_diagnostics() 
    {
        return sphereDiagnostics;
    }

    // Getter for the explicit timestepper
    SWE_Sphere_TS_ln_erk* get_ln_erk_timestepper() const
    {
        return timestepper_ln_erk;
    }

    // Getter for the simulationVariables object
    SimulationVariables* get_simulation_variables() const 
    { 
        return simVars;
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
        time.push_back(simVars->timecontrol.current_timestep_size * i_niter);
        mass.push_back(simVars->diag.total_mass);
        energy.push_back(simVars->diag.total_energy);
        potentialEnstrophy.push_back(simVars->diag.total_potential_enstrophy);
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

    // Pointer to the SimulationVariables object
    SimulationVariables *simVars;

    // Pointer to the LevelSingleton object
    LevelSingleton *levelSingleton;

    // Pointer to the ln_erk timestepper
    SWE_Sphere_TS_ln_erk* timestepper_ln_erk;

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
