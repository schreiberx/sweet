#ifndef _SPHERE_DATA_CTX_HPP_
#define _SPHERE_DATA_CTX_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <vector>
#include <sweet/core/SimulationVariables.hpp>
#include "../libpfasst_interface/LevelSingleton.hpp"

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_irk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the level singletons and the SimulationVariables object

class SphereDataCtx {

public:

    // Constructor
    SphereDataCtx(SimulationVariables *i_simVars,
                  std::vector<LevelSingleton> *i_singletons,
                  int* i_nnodes) : 
                  simVars(i_simVars),
                  levelSingletons(i_singletons)
    {
        int rank   = 0;
        int nprocs = 0;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

        if (!simVars) 
        {
            SWEETError("SphereDataCtx: simVars pointer is NULL!");
        }

        if (!levelSingletons) 
        {
            SWEETError("SphereDataCtx: levelSingletons pointer is NULL!");
        }
        
        // use first order integration in time for all pieces (only order supported)
        if (simVars->disc.timestepping_order != -1)
        {
            std::cout << "WARNING: Supplying the timestepping order manually is not supported!" << std::endl;
        }
        int timestepping_order  = 1; 

        // resize vectors
        timestepper_lg_erk_lc_n_erk.resize(levelSingletons->size());
        timestepper_lg_irk.resize(levelSingletons->size());

        // Strang splitting version that should be used:
        int version_id = 1;

        // initialize timesteppers for each level
        for (unsigned int level = 0; level < levelSingletons->size(); ++level) 
        {
            timestepper_lg_erk_lc_n_erk[level] = new SWE_Sphere_TS_lg_erk_lc_n_erk(*simVars, levelSingletons->at(level).op);
            timestepper_lg_erk_lc_n_erk[level]->setup(timestepping_order, version_id);
                
            timestepper_lg_irk[level] = new SWE_Sphere_TS_lg_irk(*simVars, levelSingletons->at(level).op);
            timestepper_lg_irk[level]->setup(timestepping_order, simVars->timecontrol.current_timestep_size);
        }
        
        // initialize the residuals
        residuals.resize(nprocs, std::vector<double>(0,0.));

        // initialize the diagnostics object
        sphereDiagnostics = new SphereHelpers_Diagnostics(
                                &(levelSingletons->back().dataConfig),
                                *simVars,
                                0);
    }

    // Destructor
    ~SphereDataCtx() 
    {
        for (auto & p : timestepper_lg_erk_lc_n_erk)
        {
            delete p;
        }
        for (auto & p : timestepper_lg_irk)
        {
            delete p;
        }
        delete sphereDiagnostics;
    }

    // Getter for the sphere data configuration at level i_level
    SphereData_Config* get_sphere_data_config(int i_level) const 
    {
        return &(levelSingletons->at(i_level).dataConfig);
    }

    // Getter for the sphere data configuration at level i_level
    BenchmarksSphereSWE* get_swe_benchmark(int i_level) const 
    {
        return &(levelSingletons->at(i_level).benchmarks);
    }

    // Getter for the sphere data configuration with no dealiasing at the fine level
    SphereData_Config* get_sphere_data_config_nodealiasing() const 
    {
        return &(levelSingletons->back().dataConfigNoDealiasing);
    }

    // Getter for the sphere data operators at level i_level
    SphereOperators_SphereData* get_sphere_operators(int i_level) const
    {
        return &(levelSingletons->at(i_level).op);
    }

    // Getter for the sphere data operators with no dealiasing at the fine level
    SphereOperators_SphereData* get_sphere_operators_nodealiasing() const
    {
        return &(levelSingletons->back().opNoDealiasing);
    }

    // Getter for the sphere diagnostics at the fine level
    SphereHelpers_Diagnostics* get_sphere_diagnostics() 
    {
        return sphereDiagnostics;
    }

    // Getter for the explicit timestepper
    SWE_Sphere_TS_lg_erk_lc_n_erk* get_lg_erk_lc_n_erk_timestepper(int i_level) const
    {
        return timestepper_lg_erk_lc_n_erk.at(i_level);
    }

    // Getter for the implicit timestepper
    SWE_Sphere_TS_lg_irk* get_lg_irk_timestepper(int i_level) const
    {
        return timestepper_lg_irk.at(i_level);
    }

    // Getter for the simulationVariables object
    SimulationVariables* get_simulation_variables() const 
    { 
        return simVars;
    }

    // Getter for the number of levels
    int get_number_of_levels() const 
    {
        return levelSingletons->size();
    }
	      
    // Save the physical invariants
    void save_physical_invariants(int i_niter) 
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
    const std::vector<std::vector<double> >& get_residuals() const { return residuals; }
    std::vector<std::vector<double> >&       get_residuals()       { return residuals; }
    
protected:

    // Pointer to the SimulationVariables object
    SimulationVariables *simVars;

    // Pointer to the LevelSingletons vector
    std::vector<LevelSingleton> *levelSingletons;

    // Pointer to the lg_erk_lc_n_erk timestepper (used for ceval_f1, ceval_f2)
    std::vector<SWE_Sphere_TS_lg_erk_lc_n_erk*> timestepper_lg_erk_lc_n_erk;

    // Pointer to the lg_irk timestepper (used for ccomp_f2)
    std::vector<SWE_Sphere_TS_lg_irk*> timestepper_lg_irk;

    // Saved Residuals for each processor
    std::vector<std::vector<double> > residuals;

    // Diagnostics (mass, energy, enstrophy)
    SphereHelpers_Diagnostics* sphereDiagnostics;

    // Some contructors and operator= are disabled
    SphereDataCtx() {};
    SphereDataCtx(const SphereDataCtx&);
    SphereDataCtx& operator=(const SphereDataCtx&);
    
    // Vectors used for plotting
    std::vector<double> time;
    std::vector<double> mass;
    std::vector<double> energy;
    std::vector<double> potentialEnstrophy;

};

#endif
