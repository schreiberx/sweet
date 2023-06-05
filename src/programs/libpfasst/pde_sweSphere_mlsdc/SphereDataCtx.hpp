#ifndef _SPHERE_DATA_CTX_HPP_
#define _SPHERE_DATA_CTX_HPP_

#include <vector>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "../interface/LevelSingleton.hpp"

#include "../../pde_sweSphere/PDESWESphere_BenchmarksCombined.hpp"
#include "../../pde_sweSphere/PDESWESphere_Diagnostics.hpp"

#include "../../pde_sweSphere/time/PDESWESphereTS_lg_erk_lc_n_erk.hpp"
#include "../../pde_sweSphere/time/PDESWESphereTS_lg_irk.hpp"


#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "../../pde_sweSphere/ShackPDESWESphere.hpp"
#include "../../pde_sweSphere/time/ShackPDESWESphereTimeDiscretization.hpp"
#include "../ShackLibPFASST.hpp"

/**
 * Class containing the context necessary to evaluate the right-hand sides
 *
 * Currently only contains a pointer to the level singletons and the sweet::ShackDictionary object

 */
class SphereDataCtx
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

public:

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
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

	SphereDataCtx()
	{
	}

	bool setup(
			sweet::ShackDictionary *i_shackDict,
			std::vector<LevelSingleton> *i_singletons,
			std::vector<int> *i_nnodes
	)
	{
		shackDict = i_shackDict;
		levelSingletons = i_singletons;

		int rank   = 0;
		int nprocs = 0;

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

		if (!shackDict)
		{
			SWEETError("SphereDataCtx: shackDict pointer is NULL!");
		}

		if (!levelSingletons)
		{
			SWEETError("SphereDataCtx: levelSingletons pointer is NULL!");
		}

		// use first order integration in time for all pieces (only order supported)
		if (shackTimeDisc->timestepping_order != -1)
		{
			std::cout << "WARNING: Supplying the timestepping order manually is not supported!" << std::endl;
		}

		int timestepping_order = 1;

		// resize vectors
		timestepper_lg_erk_lc_n_erk.resize(levelSingletons->size());
		timestepper_lg_irk.resize(levelSingletons->size());

		// Strang splitting version that should be used:
		int version_id = 1;

		// initialize timesteppers for each level
		for (unsigned int level = 0; level < levelSingletons->size(); level++)
		{
			timestepper_lg_erk_lc_n_erk[level] = new PDESWESphereTS_lg_erk_lc_n_erk;
			timestepper_lg_erk_lc_n_erk[level]->shackRegistration(shackDict);
			timestepper_lg_erk_lc_n_erk[level]->setup(&(levelSingletons->at(level).ops), timestepping_order, version_id);

			timestepper_lg_irk[level] = new PDESWESphereTS_lg_irk;
			timestepper_lg_irk[level]->shackRegistration(shackDict);
			timestepper_lg_irk[level]->setup(&(levelSingletons->at(level).ops), timestepping_order, shackTimestepControl->current_timestepSize);
		}

		// initialize the residuals
		residuals.resize(nprocs, std::vector<double>(0,0.));

		diagnostics.setup(
				&((*i_singletons)[shackLibPFASST->nlevels-1].ops),
				shackPDESWESphere,
				0
		);

		return true;
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
	}

	// Getter for the sphere data configuration at level i_level
	sweet::SphereData_Config* get_sphere_data_config(int i_level) const
	{
		return &(levelSingletons->at(i_level).sphereDataConfig);
	}

	// Getter for the sphere data configuration at level i_level
	PDESWESphere_BenchmarksCombined* get_swe_benchmark(int i_level) const
	{
		return &(levelSingletons->at(i_level).benchmarks);
	}

	// Getter for the sphere data operators at level i_level
	sweet::SphereOperators* get_sphere_operators(int i_level) const
	{
		return &(levelSingletons->at(i_level).ops);
	}

#if 0
	// Getter for the sphere data operators with no dealiasing at the fine level
	sweet::SphereOperators* get_sphere_operators_nodealiasing() const
	{
		return &(levelSingletons->back().opNoDealiasing);
	}
#endif

	// Getter for the sphere diagnostics at the fine level
	PDESWESphere_Diagnostics* get_sphere_diagnostics()
	{
		return &diagnostics;
	}

	// Getter for the explicit timestepper
	PDESWESphereTS_lg_erk_lc_n_erk* get_lg_erk_lc_n_erk_timestepper(int i_level) const
	{
		return timestepper_lg_erk_lc_n_erk.at(i_level);
	}

	// Getter for the implicit timestepper
	PDESWESphereTS_lg_irk* get_lg_irk_timestepper(int i_level) const
	{
		return timestepper_lg_irk.at(i_level);
	}

	// Getter for the simulationVariables object
	sweet::ShackDictionary* get_simulation_variables() const
	{
		return shackDict;
	}

	// Getter for the number of levels
	int get_number_of_levels() const
	{
		return levelSingletons->size();
	}

	// Save the physical invariants
	void save_physical_invariants(int i_niter)
	{
		time.push_back(shackTimestepControl->current_timestepSize * i_niter);
		mass.push_back(diagnostics.total_mass);
		energy.push_back(diagnostics.total_energy);
		potentialEnstrophy.push_back(diagnostics.total_potential_enstrophy);
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

	// Pointer to the LevelSingletons vector
	std::vector<LevelSingleton> *levelSingletons;

	// Pointer to the lg_erk_lc_n_erk timestepper (used for ceval_f1, ceval_f2)
	std::vector<PDESWESphereTS_lg_erk_lc_n_erk*> timestepper_lg_erk_lc_n_erk;

	// Pointer to the lg_irk timestepper (used for ccomp_f2)
	std::vector<PDESWESphereTS_lg_irk*> timestepper_lg_irk;

	// Saved Residuals for each processor
	std::vector<std::vector<double> > residuals;


	SphereDataCtx(const SphereDataCtx&);
	SphereDataCtx& operator=(const SphereDataCtx&);

	// Vectors used for plotting
	std::vector<double> time;
	std::vector<double> mass;
	std::vector<double> energy;
	std::vector<double> potentialEnstrophy;

};

#endif
