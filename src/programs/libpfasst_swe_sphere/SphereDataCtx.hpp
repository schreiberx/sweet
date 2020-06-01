#ifndef _SPHERE_DATA_CTX_HPP_
#define _SPHERE_DATA_CTX_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_irk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_exp.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk_n_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_erk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the level singletons and the SimulationVariables object

class SphereDataCtx {

public:
  
  // Contructor
  SphereDataCtx(
               SimulationVariables *i_simVars,
	       std::vector<LevelSingleton> *i_singletons,
	       int* i_nnodes
	       ) 
    : simVars(i_simVars),
      levelSingletons(i_singletons)
  {
    int rank   = 0;
    int nprocs = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (!simVars) 
      SWEETError("SphereDataCtx: simVars pointer is NULL!");

    if (!levelSingletons) 
      SWEETError("SphereDataCtx: levelSingletons pointer is NULL!");

    // initialize the time steppers from SWEET

    timestepper_lg_irk_lc_n_erk = new SWE_Sphere_TS_lg_irk_lc_n_erk(
								    *simVars,
								    ((*levelSingletons)[levelSingletons->size()-1].op)
								    );
    timestepper_lg_irk_lc_n_erk->setup(2,2,1);
    timestepper_l_irk_n_erk = new SWE_Sphere_TS_l_irk_n_erk(
							    *simVars,
							    ((*levelSingletons)[levelSingletons->size()-1].op)
							    );
    timestepper_l_irk_n_erk->setup(2,2,1);

    timestepper_ln_erk = new SWE_Sphere_TS_ln_erk(
						  *simVars,
						  ((*levelSingletons)[levelSingletons->size()-1].op)
						  );
    timestepper_ln_erk->setup(4);



    if (simVars->libpfasst.use_exp)
      {
	if (simVars->libpfasst.implicit_coriolis_force) 
	  {
	    timestepper_l_exp.resize(levelSingletons->size());
	    timestepper_l_erk_n_erk.resize(levelSingletons->size());
	  }
	else 
	  SWEETError("exp-based libPFASST with explicit coriolis force not implemented yet");
      }
    else
      {
	if (simVars->libpfasst.implicit_coriolis_force)
	  {
	    timestepper_l_irk.resize(levelSingletons->size());
	    timestepper_l_erk_n_erk.resize(levelSingletons->size());
	  }
	else 
	  {
	    timestepper_lg_irk.resize(levelSingletons->size());
	    if (!simVars->benchmark.use_topography)
	      timestepper_lg_erk_lc_n_erk.resize(levelSingletons->size());
	    else 
	      timestepper_lg_erk_lc_n_t_erk.resize(levelSingletons->size());
	  }
      }

    for (unsigned int level = 0; level < levelSingletons->size(); ++level) 
      {
	// select first order integration in time for explicit
	// and first order integration for implicit (only order currently supported)
	simVars->disc.timestepping_order  = 2; 
	simVars->disc.timestepping_order2 = 2; 
		
	// these timesteppers contain the functions called by LibPFASST 
	if (simVars->libpfasst.implicit_coriolis_force) 
	  {
	    timestepper_l_erk_n_erk[level] = 
	      new SWE_Sphere_TS_l_erk_n_erk(
					    *simVars,
					    ((*levelSingletons)[level].op)
					    );
	    timestepper_l_erk_n_erk[level]->setup(simVars->disc.timestepping_order,
						  simVars->disc.timestepping_order2);
	  }
	else
	  {
  	      if (simVars->benchmark.use_topography)
	        {
  	          timestepper_lg_erk_lc_n_t_erk[level] = 
	          new SWE_Sphere_TS_lg_erk_lc_n_erk(
						*simVars,
						((*levelSingletons)[level].op)
						);
	          timestepper_lg_erk_lc_n_t_erk[level]->setup(simVars->disc.timestepping_order, 0);
                }
	      else
		{
  	          timestepper_lg_erk_lc_n_erk[level] = 
	          new SWE_Sphere_TS_lg_erk_lc_n_erk(
						*simVars,
						((*levelSingletons)[level].op)
						);
	          timestepper_lg_erk_lc_n_erk[level]->setup(simVars->disc.timestepping_order,1);
                }
	  }

        if (simVars->libpfasst.use_exp)
	  {
            timestepper_l_exp[level] =
	      new SWE_Sphere_TS_l_exp(
	                               *simVars,
				       ((*levelSingletons)[level].op)
				       );

            timestepper_l_exp[level]->setup(
            			simVars->rexi,
					     "phi0",
					     simVars->timecontrol.current_timestep_size,
					     simVars->sim.sphere_use_fsphere,
					     false

					     );
 
	  }
	else
	  {
	    simVars->disc.timestepping_order  = 1; 
	    simVars->disc.timestepping_order2 = 1; 


	    if (simVars->libpfasst.implicit_coriolis_force)
	      {
		timestepper_l_irk[level] = 
		  new SWE_Sphere_TS_l_irk(
					  *simVars,
					  ((*levelSingletons)[level].op)
					  );
		
		timestepper_l_irk[level]->setup(
						simVars->disc.timestepping_order,
						simVars->timecontrol.current_timestep_size
					);
	      }
	    else
	      {
		timestepper_lg_irk[level] = 
		  new SWE_Sphere_TS_lg_irk(
					  *simVars,
					  ((*levelSingletons)[level].op)
					  );
 		
		timestepper_lg_irk[level]->setup(simVars->disc.timestepping_order,
						 simVars->timecontrol.current_timestep_size);	
	      }
	  }
      }

    // initialize the residuals
    residuals.resize(nprocs,std::vector<double>(0,0.));

    // initialize the diagnostics object
    sphereDiagnostics = new SphereHelpers_Diagnostics(
					      &((*levelSingletons)[levelSingletons->size()-1].dataConfig),
					      *simVars,
					      0
					      );
	
  }

  // Destructor
  ~SphereDataCtx() 
  {
    int m = 0;
    if (!simVars->benchmark.use_topography)
      m = (timestepper_l_erk_n_erk.size() > timestepper_lg_erk_lc_n_erk.size()) 
	? timestepper_l_erk_n_erk.size()
	: timestepper_lg_erk_lc_n_erk.size();
    else 
      m = (timestepper_l_erk_n_erk.size() > timestepper_lg_erk_lc_n_t_erk.size()) 
	? timestepper_l_erk_n_erk.size()
	: timestepper_lg_erk_lc_n_t_erk.size();


    for (int level = 0; level < m; ++level) 
    {
      if (simVars->libpfasst.use_exp)
	{
	  if (simVars->libpfasst.implicit_coriolis_force)
	    delete timestepper_l_exp[level];
	}
      else 
	{
	  if (simVars->libpfasst.implicit_coriolis_force)
	    delete timestepper_l_irk[level];
	  else
	    delete timestepper_lg_irk[level];
	}
      
      if (simVars->libpfasst.implicit_coriolis_force)
	delete timestepper_l_erk_n_erk[level];
      else 
	{
	  if (!simVars->benchmark.use_topography)
	    delete timestepper_lg_erk_lc_n_erk[level];
	  else 
	    delete timestepper_lg_erk_lc_n_t_erk[level];
	}
    }

    delete sphereDiagnostics;
  }

  // Getter for the sphere data configuration at level i_level
  SphereData_Config* get_sphere_data_config(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].dataConfig);
  }

  // Getter for the sphere data configuration at level i_level
  BenchmarksSphereSWE* get_swe_benchmark(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].benchmarks);
  }

  // Getter for the sphere data configuration with no dealiasing at the fine level
  SphereData_Config* get_sphere_data_config_nodealiasing() const 
  {
    return &((*levelSingletons)[levelSingletons->size()-1].dataConfigNoDealiasing);
  }

  // Getter for the sphere data operators at level i_level
  SphereOperators_SphereData* get_sphere_operators(
					int i_level
					) const
  {
    return &((*levelSingletons)[i_level].op);
  }

  // Getter for the sphere data operators with no dealiasing at the fine level
  SphereOperators_SphereData* get_sphere_operators_nodealiasing() const
  {
    return &((*levelSingletons)[levelSingletons->size()-1].opNoDealiasing);
  }


  // Getter for the sphere diagnostics at the fine level
  SphereHelpers_Diagnostics* get_sphere_diagnostics() 
  {
    return sphereDiagnostics;
  }

  // Getter for the linear implicit nonlinear explicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_erk_n_erk* get_l_erk_n_erk_timestepper(
							 int i_level
							 ) const
  {
    if (!simVars->libpfasst.implicit_coriolis_force)
      return NULL;
    else
      return timestepper_l_erk_n_erk[i_level];
  }

// Getter for the linear (gravitational) implicit linear (coriolis) and nonlinear explicit SWEET time stepper at level i_level
  SWE_Sphere_TS_lg_erk_lc_n_erk* get_lg_erk_lc_n_erk_timestepper(
								 int i_level
								 ) const
  {
    if (simVars->libpfasst.implicit_coriolis_force || simVars->benchmark.use_topography)
      return NULL;
    else
      return timestepper_lg_erk_lc_n_erk[i_level];
  }

  // Getter for the linear implicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_irk* get_l_irk_timestepper(
					     int i_level
					     ) const
  {
    if (simVars->libpfasst.use_exp || !simVars->libpfasst.implicit_coriolis_force)
      return NULL;
    else
      return timestepper_l_irk[i_level];
  }

  // Getter for the linear implicit SWEET time stepper at level i_level
  SWE_Sphere_TS_lg_irk* get_lg_irk_timestepper(
					       int i_level
					       ) const
  {
    if (simVars->libpfasst.use_exp || simVars->libpfasst.implicit_coriolis_force)
      return NULL;
    else
      return timestepper_lg_irk[i_level];
  }

  // Getter for the exp linear implicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_exp* get_l_exp_timestepper(
					       int i_level
					       ) const
  {
    if (!simVars->libpfasst.use_exp || !simVars->libpfasst.implicit_coriolis_force)
      return NULL;
    else
      return timestepper_l_exp[i_level];
  }
 
  // Getter for linear gravitational implicit linear Coriolis nonlinear explicit SWEET time stepper at the fine level
  SWE_Sphere_TS_lg_irk_lc_n_erk* get_lg_irk_lc_n_erk_timestepper() const
  {
    return timestepper_lg_irk_lc_n_erk;
  }

  // Getter for linear implicit nonlinear explicit SWEET time stepper at the fine level
  SWE_Sphere_TS_l_irk_n_erk* get_l_irk_n_erk_timestepper() const
  {
    return timestepper_l_irk_n_erk;
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
    return levelSingletons->size();
  }

  // Setter for the timestepping params
  void setup_time_steps(
			double i_t, 
			double i_dt
			)
  {
    simVars->timecontrol.current_timestep_size   = i_dt;
    simVars->timecontrol.current_simulation_time = i_t;
    simVars->timecontrol.max_simulation_time     = i_dt;
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
  const std::vector<std::vector<double> >& get_residuals() const { return residuals; }
  std::vector<std::vector<double> >&       get_residuals()       { return residuals; }
    
protected:

  // Pointer to the SimulationVariables object
  SimulationVariables *simVars;

  // Pointer to the LevelSingletons vector
  std::vector<LevelSingleton> *levelSingletons;

  // Pointer to the SWE_Sphere time integrator (implicit linear part, explicit nonlinear part)
  std::vector<SWE_Sphere_TS_l_erk_n_erk*>       timestepper_l_erk_n_erk;
  std::vector<SWE_Sphere_TS_lg_erk_lc_n_erk*>   timestepper_lg_erk_lc_n_erk;
  std::vector<SWE_Sphere_TS_lg_erk_lc_n_erk*> timestepper_lg_erk_lc_n_t_erk;
  std::vector<SWE_Sphere_TS_l_irk*>             timestepper_l_irk;
  std::vector<SWE_Sphere_TS_lg_irk*>            timestepper_lg_irk;
  std::vector<SWE_Sphere_TS_l_exp*>            timestepper_l_exp;

  SWE_Sphere_TS_l_irk_n_erk*     timestepper_l_irk_n_erk;
  SWE_Sphere_TS_lg_irk_lc_n_erk* timestepper_lg_irk_lc_n_erk;
  SWE_Sphere_TS_ln_erk*          timestepper_ln_erk;

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
