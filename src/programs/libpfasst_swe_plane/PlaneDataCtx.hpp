#ifndef _PLANE_DATA_CTX_HPP_
#define _PLANE_DATA_CTX_HPP_

#include <vector>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include "LevelSingleton.hpp"
#include "mpi.h"

#include "SWE_Plane_TS_l_erk_n_erk.hpp"
#include "SWE_Plane_TS_l_irk_n_erk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_ln_erk.hpp"
//#include "SWE_Plane_TS_l_rexi.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the level singletons and the SimulationVariables object

class PlaneDataCtx {

public:
  
  // Contructor
  PlaneDataCtx(
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
      SWEETError("PlaneDataCtx: simVars pointer is NULL!");

    if (!levelSingletons) 
      SWEETError("PlaneDataCtx: levelSingletons pointer is NULL!");

    // initialize the time steppers from SWEET

    timestepper_l_irk_n_erk = new SWE_Plane_TS_l_irk_n_erk(
							    *simVars,
							    ((*levelSingletons)[levelSingletons->size()-1].op)
							    );
    timestepper_l_irk_n_erk->setup(2,2);

    timestepper_ln_erk = new SWE_Plane_TS_ln_erk(
						  *simVars,
						  ((*levelSingletons)[levelSingletons->size()-1].op)
						  );
    timestepper_ln_erk->setup(4);

    if (simVars->libpfasst.use_rexi) 
      {
        //timestepper_l_rexi.resize(levelSingletons->size());
	//timestepper_l_erk_n_erk.resize(levelSingletons->size());
      }
    else
      {
	timestepper_l_irk.resize(levelSingletons->size());
	timestepper_l_erk_n_erk.resize(levelSingletons->size());
      }

    for (int level = 0; level < levelSingletons->size(); ++level) 
      {
	// select first order integration in time for explicit
	// and first order integration for implicit (only order currently supported)
	simVars->disc.timestepping_order  = 2; 
	simVars->disc.timestepping_order2 = 2; 
		
	// these timesteppers contain the functions called by LibPFASST 
	timestepper_l_erk_n_erk[level] = 
	  new SWE_Plane_TS_l_erk_n_erk(
				       *simVars,
				       ((*levelSingletons)[level].op)
				       );
	timestepper_l_erk_n_erk[level]->setup(simVars->disc.timestepping_order,
					      simVars->disc.timestepping_order2);
	
        if (simVars->libpfasst.use_rexi) 
	  {
            // timestepper_l_rexi[level] = 
	    //   new SWE_Plane_TS_l_rexi(
	    // 			      *simVars,
	    // 			      ((*levelSingletons)[level].op)
	    // 			      );

            // timestepper_l_rexi[level]->setup(simVars->rexi,
	    // 				     "phi0",
	    // 				     simVars->timecontrol.current_timestep_size
	    // 				     );
	    
	  }
	else
	  {
	    simVars->disc.timestepping_order  = 1; 
	    simVars->disc.timestepping_order2 = 1; 
	    
	    timestepper_l_irk[level] = 
	      new SWE_Plane_TS_l_irk(
				     *simVars,
				     ((*levelSingletons)[level].op)
				     );
		
	    timestepper_l_irk[level]->setup(simVars->disc.timestepping_order);
	  }
      }

    // initialize the residuals
    residuals.resize(nprocs,std::vector<double>(0,0.));

    // initialize the diagnostics object
    planeDiagnostics = new PlaneDiagnostics();
	
  }

  // Destructor
  ~PlaneDataCtx() 
  {
    int m = timestepper_l_erk_n_erk.size();

    for (int level = 0; level < m; ++level) 
    {
      //if (simVars->libpfasst.use_rexi)
	//delete timestepper_l_rexi[level];
      //else 
      delete timestepper_l_irk[level];
      
      delete timestepper_l_erk_n_erk[level];
    }
        delete planeDiagnostics;
  }

  // Getter for the plane data configuration at level i_level
  PlaneDataConfig* get_plane_data_config(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].dataConfig);
  }

  // Getter for the plane data configuration with no dealiasing at the fine level
  PlaneDataConfig* get_plane_data_config_nodealiasing() const 
  {
    return &((*levelSingletons)[levelSingletons->size()-1].dataConfigNoDealiasing);
  }

  // Getter for the plane data operators at level i_level
  PlaneOperators* get_plane_operators(
					int i_level
					) const
  {
    return &((*levelSingletons)[i_level].op);
  }

  // Getter for the plane data operators with no dealiasing at the fine level
  PlaneOperators* get_plane_operators_nodealiasing() const
  {
    return &((*levelSingletons)[levelSingletons->size()-1].opNoDealiasing);
  }


  // Getter for the plane diagnostics at the fine level
  PlaneDiagnostics* get_plane_diagnostics() 
  {
    return planeDiagnostics;
  }

  // Getter for the linear implicit nonlinear explicit SWEET time stepper at level i_level
  SWE_Plane_TS_l_erk_n_erk* get_l_erk_n_erk_timestepper(
							 int i_level
							 ) const
  {
    return timestepper_l_erk_n_erk[i_level];
  }


  // Getter for the linear implicit SWEET time stepper at level i_level
  SWE_Plane_TS_l_irk* get_l_irk_timestepper(
					     int i_level
					     ) const
  {
    if (simVars->libpfasst.use_rexi)
      return NULL;
    else
      return timestepper_l_irk[i_level];
  }
  
  // Getter for the REXI linear implicit SWEET time stepper at level i_level
  // SWE_Plane_TS_l_rexi* get_l_rexi_timestepper( 
  // 					       int i_level
  // 					       ) const
  // {
  //   if (!simVars->libpfasst.use_rexi || !simVars->libpfasst.implicit_coriolis_force)
  //     return NULL;
  //   else
  //     return timestepper_l_rexi[i_level];
  // }
 

  // Getter for linear implicit nonlinear explicit SWEET time stepper at the fine level
  SWE_Plane_TS_l_irk_n_erk* get_l_irk_n_erk_timestepper() const 
  {
    return timestepper_l_irk_n_erk;
  }


  // Getter for the explicit timestepper
  SWE_Plane_TS_ln_erk* get_ln_erk_timestepper() const 
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

  // Pointer to the SWE_Plane time integrator (implicit linear part, explicit nonlinear part)
  std::vector<SWE_Plane_TS_l_erk_n_erk*>       timestepper_l_erk_n_erk;
  std::vector<SWE_Plane_TS_l_irk*>             timestepper_l_irk;
  //std::vector<SWE_Plane_TS_l_rexi*>            timestepper_l_rexi; 

  SWE_Plane_TS_l_irk_n_erk*     timestepper_l_irk_n_erk;
  SWE_Plane_TS_ln_erk*          timestepper_ln_erk;

  // Saved Residuals for each processor
  std::vector<std::vector<double> > residuals;

  // Diagnostics (mass, energy, enstrophy)
  PlaneDiagnostics* planeDiagnostics;

  // Some contructors and operator= are disabled
  PlaneDataCtx() {};
  PlaneDataCtx(const PlaneDataCtx&);
  PlaneDataCtx& operator=(const PlaneDataCtx&);
 
  // Vectors used for plotting
  std::vector<double> time;
  std::vector<double> mass;
  std::vector<double> energy;
  std::vector<double> potentialEnstrophy;

};

#endif
