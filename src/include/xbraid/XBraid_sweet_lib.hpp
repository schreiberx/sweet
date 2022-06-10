#include <braid.hpp>
#include <parareal/Parareal_GenericData.hpp>

#if SWEET_XBRAID_SCALAR
	#include <parareal/Parareal_GenericData_Scalar.hpp>
	#include "src/programs/ode_scalar_timeintegrators/ODE_Scalar_TimeSteppers.hpp"
	typedef ODE_Scalar_TimeSteppers t_tsmType;
	#define N 1
#elif SWEET_XBRAID_PLANE
	#include <parareal/Parareal_GenericData_PlaneData_Spectral.hpp>
	#if SWEET_XBRAID_PLANE_BURGERS
		typedef Burgers_Plane_TimeSteppers t_tsmType;
		#define N 2
	#elif SWEET_XBRAID_PLANE_SWE
		typedef SWE_Plane_TimeSteppers t_tsmType;
		#define N 3
	#endif
#elif SWEET_XBRAID_SPHERE
	#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>
	typedef SWE_Sphere_TimeSteppers t_tsmType;
	#define N 3
#endif

#if SWEET_XBRAID_PLANE
	// Grid Mapping (staggered grid)
	PlaneDataGridMapping gridMapping;
#endif



class sweet_BraidVector;
class sweet_BraidApp;



/* --------------------------------------------------------------------
 * XBraid vector 
 * Stores the state of the simulation for a given time step
 * Define BraidVector, can contain anything, and be named anything
 * --> Put all time-dependent information here
 * -------------------------------------------------------------------- */
class sweet_BraidVector
{
public:
	Parareal_GenericData*	data = nullptr;

#if SWEET_XBRAID_PLANE
	PlaneDataConfig* planeDataConfig;
#elif SWEET_XBRAID_SPHERE
	SphereDataConfig* sphereDataConfig;
#endif


	sweet_BraidVector(
#if SWEET_XBRAID_PLANE
				PlaneDataConfig* i_planeDataConfig
#elif SWEET_XBRAID_SPHERE
				SphereDataConfig* i_sphereDataConfig
#endif
	)
#if SWEET_XBRAID_PLANE
		: planeDataConfig(i_planeDataConfig)
#elif SWEET_XBRAID_SPHERE
		: sphereDataConfig(i_sphereDataConfig)
#endif
	{
		this->allocate_data();
	}

	~sweet_BraidVector()
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
	}

	sweet_BraidVector(const sweet_BraidVector &i_vector)
	{
#if SWEET_XBRAID_PLANE
		this->planeDataConfig = i_vector.planeDataConfig;
#elif SWEET_XBRAID_SPHERE
		this->sphereDataConfig = i_vector.sphereDataConfig;
#endif
		*this->data = *i_vector.data;
	};

	sweet_BraidVector& operator=(const sweet_BraidVector &i_vector)
	{
#if SWEET_XBRAID_PLANE
		this->planeDataConfig = i_vector.planeDataConfig;
#elif SWEET_XBRAID_SPHERE
		this->sphereDataConfig = i_vector.sphereDataConfig;
#endif
		*this->data = *i_vector.data;
		return *this;
	};

	sweet_BraidVector operator+(
			const sweet_BraidVector &i_vector
	)	const
	{
#if SWEET_XBRAID_SCALAR
		sweet_BraidVector out;
#elif SWEET_XBRAID_PLANE
		sweet_BraidVector out(this->planeDataConfig);
#elif SWEET_XBRAID_SPHERE
		sweet_BraidVector out(this->sphereDataConfig);
#endif
		*out.data = *this->data;
		*out.data += *i_vector.data;
		return out;
	}

	sweet_BraidVector operator*(
			const double i_value
	)	const
	{
#if SWEET_XBRAID_SCALAR
		sweet_BraidVector out;
#elif SWEET_XBRAID_PLANE
		sweet_BraidVector out(this->planeDataConfig);
#elif SWEET_XBRAID_SPHERE
		sweet_BraidVector out(this->sphereDataConfig);
#endif
		*out.data = *this->data;
		*out.data *= i_value;
		return out;
	}


	void allocate_data()
	{
#if SWEET_XBRAID_SCALAR
		{
			data = new Parareal_GenericData_Scalar<N>;
			data->allocate_data();
			//this->set_time(this->timeframe_end);
		}

#elif SWEET_XBRAID_PLANE
		{
			data = new Parareal_GenericData_PlaneData_Spectral<N>;
			data->setup_data_config(this->planeDataConfig);
			data->allocate_data();
			//this->setup_data_config(this->planeDataConfig);
			//this->set_time(this->timeframe_end);
		}

#elif SWEET_XBRAID_SPHERE
		{
			data = new Parareal_GenericData_SphereData_Spectral<N>;
			data->setup_data_config(this->sphereDataConfig);
			data->allocate_data();
			//this->setup_data_config(this->sphereDataConfig);
			//this->set_time(this->timeframe_end);
		}
#endif
	}



};


/**
 * operator to support operations such as:
 *
 * 1.5 * arrayData;
 *
 * Otherwise, we'd have to write it as arrayData*1.5
 *
 */
inline
static
sweet_BraidVector operator*(
		const double i_value,
		const sweet_BraidVector &i_vector
)
{
	return i_vector * i_value;
}



// Wrapper for BRAID's App object·
// --> Put all time INDEPENDENT information here
class sweet_BraidApp
			: public BraidApp
{

public:

	SimulationVariables*		simVars;
	double				dt;
	std::vector<t_tsmType*>		timeSteppers;

#if SWEET_XBRAID_PLANE
	PlaneOperators*		op_plane;
	PlaneDataConfig*	planeDataConfig;
#elif SWEET_XBRAID_SPHERE
	SphereOperators*	op_sphere;
	SphereDataConfig*	sphereDataConfig;
#endif


	int			size_buffer;

	// We will need the MPI Rank
	int rank;
public:

	// Constructor·
	sweet_BraidApp( MPI_Comm		i_comm_t,
			int			i_rank,
			double			i_tstart,
			double			i_tstop,
			int 			i_ntime,
			SimulationVariables*	i_simVars
#if SWEET_XBRAID_PLANE
			,
			PlaneDataConfig* i_planeDataConfig,
			PlaneOperators* i_op_plane
#elif SWEET_XBRAID_SPHERE
			,
			SphereDataConfig* i_sphereDataConfig,
			SphereOperators* i_op_sphere
#endif
			)
		:
			BraidApp(i_comm_t, i_tstart, i_tstop, i_ntime),
			rank(i_rank),
			simVars(i_simVars)
#if SWEET_XBRAID_PLANE
			,
			planeDataConfig(i_planeDataConfig),
			op_plane(i_op_plane)
#elif SWEET_XBRAID_SPHERE
			,
			sphereDataConfig(i_sphereDataConfig)
			op_sphere(i_op_sphere)
#endif
	{
	}

	virtual ~sweet_BraidApp()
	{
	}

public:
	void setup(BraidCore& i_core)
	{


		/////////////////////////////////////////////////
		// get parameters from simVars and set to Core //
		/////////////////////////////////////////////////

		i_core.SetMaxLevels(this->simVars->xbraid.xbraid_max_levels);
		/////i_core.SetIncrMaxLevels();

		i_core.SetSkip(this->simVars->xbraid.xbraid_skip);

		i_core.SetMinCoarse(this->simVars->xbraid.xbraid_min_coarse);

		///i_core.SetRelaxOnlyCG(this->simVars->xbraid.xbraid_relax_only_cg);

		i_core.SetNRelax(-1, this->simVars->xbraid.xbraid_nrelax);
		if (this->simVars->xbraid.xbraid_nrelax0 > -1)
			i_core.SetNRelax(0, this->simVars->xbraid.xbraid_nrelax0);

		i_core.SetAbsTol(this->simVars->xbraid.xbraid_tol);
		i_core.SetRelTol(this->simVars->xbraid.xbraid_tol);

		i_core.SetTemporalNorm(this->simVars->xbraid.xbraid_tnorm);

		i_core.SetCFactor(-1, this->simVars->xbraid.xbraid_cfactor);
		if (this->simVars->xbraid.xbraid_cfactor0 > -1)
			i_core.SetCFactor(0, this->simVars->xbraid.xbraid_cfactor0);

		///i_core.SetiPeriodic(this->simVars->xbraid.xbraid_periodic);

		////i_core.SetResidual();

		i_core.SetMaxIter(this->simVars->xbraid.xbraid_max_iter);

		i_core.SetPrintLevel(this->simVars->xbraid.xbraid_print_level);

		i_core.SetSeqSoln(this->simVars->xbraid.xbraid_use_seq_soln);

		i_core.SetAccessLevel(this->simVars->xbraid.xbraid_access_level);

		i_core.SetNFMG(this->simVars->xbraid.xbraid_fmg);

		//i_core.SetiNFMGVcyc(this->simVars->xbraid.xbraid_fmgvcyc);

		i_core.SetStorage(this->simVars->xbraid.xbraid_storage);

		//i_core.SetRevertedRanks(this->simVars->xbraid.xbraid_reverted_ranks);

		////i_core.SetRefine(this->simVars->xbraid.xbraid_refine);
		///i_core.SetMaxRefinements(this->simVars->xbraid.xbraid_max_Refinements);


		this->setup();
	}

public:
	void setup()
	{
		////////////////////////////////////
		// get tsm and tso for each level //
		////////////////////////////////////
		std::vector<std::string> tsms = this->getTimeSteppingMethodFromParameters();
		std::vector<int> tsos = this->getTimeSteppingOrderFromParameters(1);
		std::vector<int> tsos2 = this->getTimeSteppingOrderFromParameters(2);

		// create a timeSteppers instance for each level
		for (int level = 0; level < this->simVars->xbraid.xbraid_max_levels; level++)
		{

#if SWEET_XBRAID_SCALAR
			ODE_Scalar_TimeSteppers* tsm = new ODE_Scalar_TimeSteppers;
			tsm->setup(
					//tsms[level],
					//tsos[level],
					*this->simVars
				);
#elif SWEET_XBRAID_PLANE
	#if SWEET_XBRAID_PLANE_SWE
			SWE_Plane_TimeSteppers* tsm = new SWE_Plane_TimeSteppers;
			tsm->setup(
					tsms[level],
					tsos[level],
					tsos2[level],
					*this->op_plane,
					*this->simVars
				);
	#elif SWEET_XBRAID_PLANE_BURGERS
			Burgers_Plane_TimeSteppers* tsm = new Burgers_Plane_TimeSteppers;
			tsm->setup(
					tsms[level],
					tsos[level],
					tsos2[level],
					*this->op_plane,
					*this->simVars
				);
	#endif
#elif SWEET_XBRAID_SPHERE
			SWE_Sphere_TimeSteppers* tsm = new SWE_Sphere_TimeSteppers;
			tsm->setup(
						tsms[level],
						*this->op_sphere,
						*this->simVars
					);
#endif

			this->timeSteppers.push_back(tsm);

		}


		// get buffer size
#if SWEET_XBRAID_SCALAR
		this->size_buffer = N * sizeof(double);
#elif SWEET_XBRAID_PLANE
		this->size_buffer = N * planeDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
#elif SWEET_XBRAID_SPHERE
		this->size_buffer = N * sphereDataConfig->spectral_array_data_number_of_elements * sizeof(std::complex<double>);
#endif

	}


private:
	/*
	 * Get timestepping methods for each of the N levels
	 * Input string must contain 1, 2 or N tsm names separated by comma:
	 *  - 1 tsm name: same tsm for all levels;
	 *  - 2 tsm names: first one is for level 0 (finest); all coarse levels use the second one
	 *  - N tsm names: i-th level uses the i-th tsm.
	 */
	std::vector<std::string> getTimeSteppingMethodFromParameters()
	{
		std::vector<std::string> tsm = {};
		std::stringstream all_tsm = std::stringstream(this->simVars->xbraid.xbraid_timestepping_method);

		while (all_tsm.good())
		{
			std::string str;
			getline(all_tsm, str, ',');
			tsm.push_back(str);
		}

		if ( ! (tsm.size() == 1 || tsm.size() == 2 || tsm.size() == this->simVars->xbraid.xbraid_max_levels ) )
			SWEETError("xbraid_timestepping_method must contain 1, 2 or N timestepping names.");

		// all levels use same tsm
		if (tsm.size() == 1)
			for (int level = 1; level < this->simVars->xbraid.xbraid_max_levels; level++)
				tsm.push_back(tsm[0]);

		// all coarse levels use same tsm
		if (tsm.size() == 2)
			for (int level = 2; level < this->simVars->xbraid.xbraid_max_levels; level++)
				tsm.push_back(tsm[1]);


		return tsm;
	}

	/*
	 * Get timestepping order for each of the N levels
	 * Input string must contain 1, 2 or N orders separated by comma:
	 *  - 1 tso: same tso for all levels;
	 *  - 2 tso: first one is for level 0 (finest); all coarse levels use the second one
	 *  - N tso: i-th level uses the i-th tso.
	 */
	std::vector<int> getTimeSteppingOrderFromParameters(int o)
	{
		std::vector<int> tso = {};
		std::stringstream all_tso;
		if (o == 1)
			all_tso = std::stringstream(this->simVars->xbraid.xbraid_timestepping_order);
		else if (o == 2)
			all_tso = std::stringstream(this->simVars->xbraid.xbraid_timestepping_order2);
		else
			SWEETError("Wrong parameter for getting timestepping order.");

		while (all_tso.good())
		{
			std::string str;
			getline(all_tso, str, ',');
			tso.push_back(stoi(str));
		}

		if ( ! (tso.size() == 1 || tso.size() == 2 || tso.size() == this->simVars->xbraid.xbraid_max_levels ) )
			SWEETError("xbraid_timestepping_order must contain 1, 2 or N timestepping orders.");

		// all levels use same tso
		if (tso.size() == 1)
			for (int level = 1; level < this->simVars->xbraid.xbraid_max_levels; level++)
				tso.push_back(tso[0]);

		// all coarse levels use same tso
		if (tso.size() == 2)
			for (int level = 2; level < this->simVars->xbraid.xbraid_max_levels; level++)
				tso.push_back(tso[1]);

		return tso;
	}



public:
	sweet_BraidVector* create_new_vector()
	{
#if SWEET_XBRAID_SCALAR
		sweet_BraidVector* U = new sweet_BraidVector;
#elif SWEET_XBRAID_PLANE
		sweet_BraidVector* U = new sweet_BraidVector(this->planeDataConfig);
#elif SWEET_XBRAID_SPHERE
		sweet_BraidVector* U = new sweet_BraidVector(this->sphereDataConfig);
#endif
		return U;
	}


public:
	/* --------------------------------------------------------------------
	 * Time integrator routine that performs the update
	 *   u_i = Phi_i(u_{i-1}) + g_i 
	 * 
	 * When Phi is called, u is u_{i-1}.
	 * The return value is that u is set to u_i upon completion
	 * -------------------------------------------------------------------- */
	braid_Int
	Step(
			braid_Vector		io_U,
			braid_Vector		i_ustop,
			braid_Vector		i_fstop,
			BraidStepStatus&	io_status
			)
	{

		sweet_BraidVector* U = (sweet_BraidVector*) io_U;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;
		int nlevels;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
		io_status.GetLevel(&level);
		io_status.GetNLevels(&nlevels);

		/////std::cout << "NLevels " << nlevels << std::endl;
		/////std::cout << "tstart tstop " << tstart << " " << tstop << std::endl;
		/////std::cout << "level " << level << std::endl;

		this->simVars->timecontrol.current_simulation_time = tstart;
		this->simVars->timecontrol.current_timestep_size = tstop - tstart;

		this->timeSteppers[level]->master->run_timestep(
								U->data,
								this->simVars->timecontrol.current_timestep_size,
								this->simVars->timecontrol.current_simulation_time
		);

		/* Tell XBraid no refinement */
		io_status.SetRFactor(1);

		return 0;
	}
		
		/* --------------------------------------------------------------------
		 * -------------------------------------------------------------------- */

	virtual braid_Int
	Residual(
				braid_Vector			i_ustop,
				braid_Vector			o_r,
				BraidStepStatus&		io_status
		)
	{

		//braid_StepStatus& status = (braid_StepStatus&) io_status;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;

		/* Grab status of current time step */
		io_status.GetTstartTstop(&tstart, &tstop);
	
		/* Grab level */
		io_status.GetLevel(&level);
	
		/* Set the new dt in the user's manager*/
		this->dt = tstop - tstart;


		sweet_BraidVector* u = (sweet_BraidVector*) i_ustop;
		sweet_BraidVector* r = (sweet_BraidVector*) o_r;
		std::cout << "SOLUTION " << tstop << " " << u->data->get_pointer_to_data_Scalar()->simfields[0] << std::endl;
		std::cout << "RESIDUAL " << tstop << " " << r->data->get_pointer_to_data_Scalar()->simfields[0] << std::endl;

		/////////////////* Now, set up the discretization matrix.  Use the XBraid level to index
		//////////////// * into the matrix lookup table */
		////////////////if( app->dt_A[level] == -1.0 ){
		////////////////	app->nA++;
		////////////////	app->dt_A[level] = tstop-tstart;
	
		////////////////	setUpImplicitMatrix( app->man );
		////////////////	app->A[level] = app->man->A;
	
		////////////////	/* Set up the PFMG solver using r->x as dummy vectors. */
		////////////////	setUpStructSolver( app->man, r->x, r->x );
		////////////////	app->solver[level] = app->man->solver;
		////////////////}
	
		/////////////////* Set the correct solver for the current level */
		////////////////app->man->timeSteppers = app->timeSteppers[solver];
	
		/////////////////* Compute residual Ax */
		//////////////////////app->man->A = app->A[level];
		////////////////comp_res(app->man, ustop->x, r->x, tstart, tstop);
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a vector object for a given time point.
	 * This function is only called on the finest level.
	 * -------------------------------------------------------------------- */
	braid_Int
	Init(
			double		i_t,
			braid_Vector*	o_U
			)
	{

		sweet_BraidVector* U = create_new_vector();

		if( i_t == this->tstart )
		{
			std::cout << "Setting initial solution " << std::endl;
	#if SWEET_XBRAID_SCALAR
			double u0 = atof(simVars->bogus.var[1].c_str());
			U->data->dataArrays_to_GenericData_Scalar(u0);
	
	#elif SWEET_XBRAID_PLANE
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			PlaneData_Spectral t0_prog_v(planeDataConfig);
	
		#if SWEET_XBRAID_PLANE_SWE
			SWEPlaneBenchmarksCombined swePlaneBenchmarks;
			swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, *simVars, *op_plane);
		#elif SWEET_XBRAID_PLANE_BURGERS
			PlaneData_Physical t0_prog_u_phys(t0_prog_u.planeDataConfig);
			PlaneData_Physical t0_prog_v_phys(t0_prog_v.planeDataConfig);
			if (simVars->disc.space_grid_use_c_staggering)
			{
				t0_prog_u_phys.physical_update_lambda_array_indices(
							[&](int i, int j, double &io_data)
					{
						double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
						double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
						io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
					}
				);
				t0_prog_v_phys.physical_update_lambda_array_indices(
							[&](int i, int j, double &io_data)
					{
					io_data = 0.0;
					}
				);
			}
			else
			{
				t0_prog_u_phys.physical_update_lambda_array_indices(
							[&](int i, int j, double &io_data)
					{
						double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
						double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
						io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
					}
				);
		
				t0_prog_v_phys.physical_update_lambda_array_indices(
							[&](int i, int j, double &io_data)
					{
					io_data = 0.0;
					}
				);
			}
			t0_prog_u.loadPlaneDataPhysical(t0_prog_u_phys);
			t0_prog_v.loadPlaneDataPhysical(t0_prog_v_phys);
		#endif
	
	
			U->data->dataArrays_to_GenericData_PlaneData_Spectral(
		#if SWEET_XBRAID_PLANE_SWE
										t0_prog_h_pert,
		#endif
										t0_prog_u,
										t0_prog_v);
	
	#elif SWEET_XBRAID_SPHERE
			SphereData_Spectral t0_prog_phi_pert(sphereDataConfig);
			SphereData_Spectral t0_prog_vrt(sphereDataConfig);
			SphereData_Spectral t0_prog_div(sphereDataConfig);
	
			BenchmarksSphereSWE sphereBenchmarks;
			sphereBenchmarks.setup(*simVars, *op_sphere);
			sphereBenchmarks.master->get_initial_state(t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	
			Parareal_SimulationInstace::dataArrays_to_GenericData_SphereData_Spectral(U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
	
	
		}
		else if (this->simVars->xbraid.xbraid_use_rand)
		{
			std::cout << "Setting random solution " << std::endl;
	#if SWEET_XBRAID_SCALAR
			double u0 = ((double)braid_Rand())/braid_RAND_MAX;
			U->data->dataArrays_to_GenericData_Scalar(u0);
	#elif SWEET_XBRAID_PLANE
	
		#if SWEET_XBRAID_PLANE_SWE
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			PlaneData_Physical t0_prog_h_phys(planeDataConfig);
		#endif
	
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			PlaneData_Spectral t0_prog_v(planeDataConfig);
	
			PlaneData_Physical t0_prog_u_phys(planeDataConfig);
			PlaneData_Physical t0_prog_v_phys(planeDataConfig);
	
		#if SWEET_XBRAID_PLANE_SWE
			t0_prog_h_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = this->simVars->sim.h0 + ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
		#endif
			t0_prog_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
			t0_prog_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
	
		#if SWEET_XBRAID_PLANE_SWE
			t0_prog_h_pert.loadPlaneDataPhysical(t0_prog_h_phys);
		#endif
			t0_prog_u.loadPlaneDataPhysical(t0_prog_u_phys);
			t0_prog_v.loadPlaneDataPhysical(t0_prog_v_phys);
	
			U->data->dataArrays_to_GenericData_PlaneData_Spectral(
		#if SWEET_XBRAID_PLANE_SWE
										t0_prog_h_pert,
		#endif
										t0_prog_u,
										t0_prog_v);
	
	
	#elif SWEET_XBRAID_SPHERE
			SphereData_Spectral t0_prog_phi_pert(sphereDataConfig);
			SphereData_Spectral t0_prog_vrt(sphereDataConfig);
			SphereData_Spectral t0_prog_div(sphereDataConfig);
	
			SphereData_Physical t0_prog_phi_pert_phys(sphereDataConfig);
			SphereData_Physical t0_prog_vrt_phys(sphereDataConfig);
			SphereData_Physical t0_prog_div_phys(sphereDataConfig);
	
			t0_prog_phi_pert_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = this->simVars->sim.h0 + ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
			t0_prog_vrt_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
			t0_prog_div_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			);
	
			t0_prog_phi_pert.loadPlaneDataPhysical(t0_prog_phi_pert_phys);
			t0_prog_vrt.loadPlaneDataPhysical(t0_prog_vrt_phys);
			t0_prog_div.loadPlaneDataPhysical(t0_prog_div_phys);
	
			this->dataArrays_to_GenericData_SphereData_Spectral(U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
		}
		else
		{
			std::cout << "Setting zero solution " << std::endl;
			/* Sets U as an all zero vector*/
	#if SWEET_XBRAID_SCALAR
			double zero = 0;
			U->data->dataArrays_to_GenericData_Scalar(zero);
	#elif SWEET_XBRAID_PLANE
		#if SWEET_XBRAID_PLANE_SWE
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			t0_prog_h_pert.spectral_set_zero();
		#endif
	
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			t0_prog_u.spectral_set_zero();
	
			PlaneData_Spectral t0_prog_v(planeDataConfig);
			t0_prog_v.spectral_set_zero();

			U->data->dataArrays_to_GenericData_PlaneData_Spectral(
		#if SWEET_XBRAID_PLANE_SWE
										t0_prog_h_pert,
		#endif
										t0_prog_u,
										t0_prog_v);
	
	#elif SWEET_XBRAID_SPHERE
			SphereData_Spectral t0_prog_phi_pert(sphereDataConfig);
			SphereData_Spectral t0_prog_vrt(sphereDataConfig);
			SphereData_Spectral t0_prog_div(sphereDataConfig);
	
			t0_prog_phi_pert.spectral_set_zero();
			t0_prog_vrt.spectral_set_zero();
			t0_prog_div_pert.spectral_set_zero();
	
			this->dataArrays_to_GenericData_SphereData_Spectral(U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
		}

		*o_U = (braid_Vector) U;

		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a copy of a vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Clone(
			braid_Vector	i_U,
			braid_Vector*	o_V
		)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;
		sweet_BraidVector* V = create_new_vector();
		*V = *U;
		*o_V = (braid_Vector) V;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Destroy vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Free(
			braid_Vector	i_U)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;
		delete U;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * Compute vector sum y = alpha*x + beta*y.
	 * -------------------------------------------------------------------- */
	braid_Int
	Sum(
			double			i_alpha,
			braid_Vector		i_X,
			double			i_beta,
			braid_Vector		io_Y
		)
	{

		sweet_BraidVector* X = (sweet_BraidVector*) i_X;
		sweet_BraidVector* Y = (sweet_BraidVector*) io_Y;

		*Y = *X * i_alpha + *Y * i_beta;

		return 0;
	}


	/* --------------------------------------------------------------------
	 * User access routine to spatial solution vectors and allows for user
	 * output.  The default XBraid parameter of access_level=1, calls 
	 * my_Access only after convergence and at every time point.
	 * -------------------------------------------------------------------- */
	braid_Int
	Access(
				braid_Vector		i_U,
				BraidAccessStatus&	io_astatus
			)
	{
		double     tstart         = (this->tstart);
		double     tstop          = (this->tstop);
		////int        nt             = (this->nt);
	
		double     rnorm, disc_err, t;
		int        iter, level, done, index, myid, it;
		char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];


		sweet_BraidVector *U = (sweet_BraidVector*) i_U;

		/* Retrieve current time from Status Object */
		////braid_AccessStatusGetT(astatus, &t);
		io_astatus.GetT(&t);
		io_astatus.GetTIndex(&it);
		io_astatus.GetIter(&iter);
		io_astatus.GetLevel(&level);
	
		/* Retrieve XBraid State Information from Status Object */
		///////////MPI_Comm_rank(app->comm_x, &myid);
		///////////braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
		///////////braid_AccessStatusGetResidual(astatus, &rnorm);


		if(level == 0)
		{

			// Output physical solution to file
			if (simVars->xbraid.xbraid_store_iterations)
				this->output_data_file(
							U->data,
							iter,
							it,
							t,
							false
				);
			// Compute and store errors w.r.t. ref solution
			if (simVars->xbraid.xbraid_load_ref_csv_files)
				this->store_parareal_error(
								k + 1,
								i,
								pVars->path_ref_csv_files,
								"ref"
				);
			// Compute and store errors w.r.t. fine (serial) solution
			if (simVars->xbraid.xbraid_load_fine_csv_files)
				this->store_parareal_error(
								k + 1,
								i,
								pVars->path_fine_csv_files,
								"fine"
				);


			// Store residual
			if (it == 0) {
				double res;
				io_astatus.GetResidual(&res);
				this->output_residual_file(res,
								iter);
			}

		return 0;
	}



	/* --------------------------------------------------------------------
	 * Compute norm of a spatial vector 
	 * -------------------------------------------------------------------- */
	braid_Int
	SpatialNorm(
			braid_Vector  i_U,
			double*       o_norm)
	{
		sweet_BraidVector* U = (sweet_BraidVector*) i_U;
		*o_norm = U->data->reduce_norm2();
		return 0;
	}


	/* --------------------------------------------------------------------
	 * Return buffer size needed to pack one spatial braid_Vector.  Here the
	 * vector contains one double at every grid point and thus, the buffer 
	 * size is the number of grid points.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufSize(
			int*			o_size,
			BraidBufferStatus&	o_status)
	{
		*o_size = this->size_buffer;
		return 0;
	}


	/* --------------------------------------------------------------------
	 * Pack a braid_Vector into a buffer.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufPack(
			braid_Vector		i_U,
			void*			o_buffer,
			BraidBufferStatus&	o_status
		)
	{

		sweet_BraidVector* U = (sweet_BraidVector*) i_U;

#if SWEET_XBRAID_SCALAR
		double* dbuffer = (double*) o_buffer;
#else
		std::complex<double>* dbuffer = (std::complex<double>*) o_buffer;
#endif

		U->data->serialize(dbuffer);

		o_status.SetSize( this->size_buffer );
		return 0;
	}

	/* --------------------------------------------------------------------
	 * Unpack a buffer and place into a braid_Vector
	 * -------------------------------------------------------------------- */
	braid_Int
	BufUnpack(
			void*			i_buffer,
			braid_Vector*		o_U,
			BraidBufferStatus&	status
		)
	{

		sweet_BraidVector* U = create_new_vector();

#if SWEET_XBRAID_SCALAR
		double* dbuffer = (double*) i_buffer;
#else
		std::complex<double>* dbuffer = (std::complex<double>*) i_buffer;
#endif

		U->data->deserialize(dbuffer);

		*o_U = (braid_Vector) U;

		return 0;
	}


///////////////////////////
// FILE OUTPUT FUNCTIONS //
///////////////////////////
// Copied from Parareal; find an unified formulation!

	void output_residual_file(
			double res,
			int iteration_id
	)
	{

		char buffer[1024];

		const char* filename_template = "residual_iter%03d.csv";
		sprintf(buffer, filename_template, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(16);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << res;

		file.close();

	}


	void output_data_file(
			Parareal_GenericData* i_data,
			int iteration_id,
			int time_slice_id,
			double t,
			bool output_initial_data = false
	)
	{
#if SWEET_XBRAID_SCALAR
		double u_out;
		////if (output_initial_data)
		////	this->GenericData_Scalar_to_dataArrays(this->parareal_data_start, u_out);
		////else
		////	this->GenericData_Scalar_to_dataArrays(this->parareal_data_output, u_out);
		i_data->GenericData_Scalar_to_dataArrays(u_out);

		// Dump  data in csv, if output filename is not empty
		if (simVars->iodata.output_file_name.size() > 0)
		{
			std::string output_filenames = "";
			output_filenames = write_file_xbraid_scalar(u_out, "prog_u", iteration_id, t, output_initial_data);
		}

#elif SWEET_XBRAID_PLANE
	#if SWEET_XBRAID_PLANE_BURGERS

		PlaneData_Spectral dummy(this->planeDataConfig);
		PlaneData_Spectral u_out(this->planeDataConfig);
		PlaneData_Spectral v_out(this->planeDataConfig);
		/////if (output_initial_data)
		/////	this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_start, u_out, v_out);
		/////else
		/////	this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_output, u_out, v_out);
		i_data->GenericData_PlaneData_Spectral_to_dataArrays(u_out, v_out);

		PlaneData_Physical u_out_phys = u_out.toPhys();
		PlaneData_Physical v_out_phys = v_out.toPhys();

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		PlaneData_Physical t_u(planeDataConfig);
		PlaneData_Physical t_v(planeDataConfig);

		if (simVars->disc.space_grid_use_c_staggering) // Remap in case of C-grid
		{
			gridMapping.mapCtoA_u(u_out_phys, t_u);
			gridMapping.mapCtoA_v(v_out_phys, t_v);
		}
		else
		{
			t_u = u_out_phys;
			t_v = v_out_phys;
		}

		// Dump  data in csv, if output filename is not empty
		if (simVars->iodata.output_file_name.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_xbraid_plane(t_u, "prog_u", iteration_id, output_initial_data);
			output_filenames += ";" + write_file_xbraid_plane(t_v, "prog_v", iteration_id, output_initial_data);

		}

		write_file_spec_amp_phase_xbraid_plane(u_out, "prog_u", iteration_id, output_initial_data);

		if (simVars->misc.compute_errors)
		{
			PlaneData_Spectral ana = compute_errors2(u_out, v_out);

			write_file_xbraid_plane(ana.toPhys(),"analytical",iteration_id,time_slice_id);
			write_file_spec_amp_phase_xbraid_plane(ana.toPhys(), "analytical", iteration_id, output_initial_data);
		}

	#elif SWEET_XBRAID_PLANE_SWE

		PlaneData_Spectral h_out(this->planeDataConfig);
		PlaneData_Spectral u_out(this->planeDataConfig);
		PlaneData_Spectral v_out(this->planeDataConfig);
		///if (output_initial_data)
		///	this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_start, h_out, u_out, v_out);
		///else
		///	this->GenericData_PlaneData_Spectral_to_dataArrays(this->parareal_data_output, h_out, u_out, v_out);
		i_data->GenericData_PlaneData_Spectral_to_dataArrays(h_out, u_out, v_out);

		PlaneData_Physical h_out_phys = h_out.toPhys();
		PlaneData_Physical u_out_phys = u_out.toPhys();
		PlaneData_Physical v_out_phys = v_out.toPhys();

		// Save .vtk files for visualizing in paraview
		std::ostringstream ss2;
		if (output_initial_data)
			ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
		else
			ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		std::string filename2 = ss2.str();
		h_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		PlaneData_Physical t_h(planeDataConfig);
		PlaneData_Physical t_u(planeDataConfig);
		PlaneData_Physical t_v(planeDataConfig);

		if (simVars->disc.space_grid_use_c_staggering) // Remap in case of C-grid
		{
			t_h = h_out_phys;
			gridMapping.mapCtoA_u(u_out_phys, t_u);
			gridMapping.mapCtoA_v(v_out_phys, t_v);
		}
		else
		{
			t_h = h_out_phys;
			t_u = u_out_phys;
			t_v = v_out_phys;
		}

		// Dump  data in csv, if output filename is not empty
		if (simVars->iodata.output_file_name.size() > 0)
		{
			std::string output_filenames = "";

			output_filenames = write_file_xbraid_plane(t_h, "prog_h_pert", iteration_id, t, output_initial_data);
			output_filenames += ";" + write_file_xbraid_plane(t_u, "prog_u", iteration_id, t, output_initial_data);
			output_filenames += ";" + write_file_xbraid_plane(t_v, "prog_v", iteration_id, t, output_initial_data);

			output_filenames += ";" + write_file_xbraid_plane(op_plane->ke(t_u,t_v).toPhys(),"diag_ke", iteration_id, t, output_initial_data);

			output_filenames += ";" + write_file_spec_xbraid_plane(op_plane->ke(t_u,t_v).toPhys(),"diag_ke_spec", iteration_id, t, output_initial_data);

			output_filenames += ";" + write_file_xbraid_plane(op_plane->vort(t_u, t_v).toPhys(), "diag_vort", iteration_id, t, output_initial_data);
			output_filenames += ";" + write_file_xbraid_plane(op_plane->div(t_u, t_v).toPhys(), "diag_div", iteration_id, t, output_initial_data);

			/////////if(this->compute_normal_modes){
			/////////	SWEETError("TODO");
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.geo, "nm_geo", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igwest, "nm_igwest", iteration_id, output_initial_data);
			/////////	///output_filenames += ";" + write_file_spec_parareal(normalmodes.igeast, "nm_igeast", iteration_id, output_initial_data);
			/////////}
			
		}
	#endif


#elif SWEET_XBRAID_SPHERE

		SphereData_Spectral phi_out(this->sphereDataConfig);
		SphereData_Spectral vrt_out(this->sphereDataConfig);
		SphereData_Spectral div_out(this->sphereDataConfig);
		if (output_initial_data)
			this->GenericData_SphereData_Spectral_to_dataArrays(this->parareal_data_start, phi_out, vrt_out, div_out);
		else
			this->GenericData_SphereData_Spectral_to_dataArrays(this->parareal_data_output, phi_out, vrt_out, div_out);
	
		SphereData_Physical phi_out_phys = phi_out.toPhys();
		SphereData_Physical vrt_out_phys = vrt_out.toPhys();
		SphereData_Physical div_out_phys = div_out.toPhys();
	
		///////////////// Save .vtk files for visualizing in paraview
		///////////////std::ostringstream ss2;
		///////////////if (output_initial_data)
		///////////////	ss2 << "output_slice" << time_slice_id - 1 << "_iter" << iteration_id << ".vtk";
		///////////////else
		///////////////	ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";
		///////////////std::string filename2 = ss2.str();
		///////////////phi_out_phys.file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// Dump  data in csv, if output filename is not empty
		if (simVars->iodata.output_file_name.size() > 0)
		{
			if (simVars->iodata.output_file_mode == "csv")
			{
				std::string output_filename;
	
				SphereData_Spectral h = phi_out_phys*(1.0/simVars->sim.gravitation);
				h += simVars->sim.h0;
	
				output_filename = write_file_csv_parareal_sphere(h, "prog_h", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;
	
				output_filename = write_file_csv_parareal_sphere(phi_out, "prog_phi_pert", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << " (min: " << phi_out_phys.physical_reduce_min() << ", max: " << phi_out_phys.physical_reduce_max() << ")" << std::endl;
	
				SphereData_Physical u(sphereDataConfig);
				SphereData_Physical v(sphereDataConfig);
	
				op_sphere->vrtdiv_to_uv(vrt_out_phys, div_out_phys, u, v);
	
				output_filename = write_file_csv_parareal_sphere(u, "prog_u", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_parareal_sphere(v, "prog_v", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_parareal_sphere(vrt_out, "prog_vrt", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << std::endl;
	
				output_filename = write_file_csv_parareal_sphere(div_out, "prog_div", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << std::endl;
	
				SphereData_Spectral potvrt = (phi_out/simVars->sim.gravitation)*vrt_out;
	
				output_filename = write_file_csv_parareal_sphere(potvrt, "prog_potvrt", iteration_id, output_initial_data);
				std::cout << " + " << output_filename << std::endl;
			}
			else if (simVars->iodata.output_file_mode == "bin")
			{
				std::string output_filename;
	
				{
					output_filename = write_file_bin_parareal_sphere(phi_out, "prog_phi_pert", iteration_id, output_initial_data);
					SphereData_Physical prog_phys = phi_out.toPhys();
	
					std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_parareal_sphere(vrt_out, "prog_vrt", iteration_id, output_initial_data);
					SphereData_Physical prog_phys = vrt_out.toPhys();
	
					std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
	
				{
					output_filename = write_file_bin_parareal_sphere(div_out, "prog_div", iteration_id, output_initial_data);
					SphereData_Physical prog_phys = div_out.toPhys();
	
					std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
				}
			}
			
		}


#endif
	};


#if SWEET_XBRAID_PLANE_BURGERS
	// For Burgers
	PlaneData_Spectral compute_errors2(
         const PlaneData_Spectral &i_planeData_u,
         const PlaneData_Spectral &i_planeData_v
	)
	{

		int analytic_solution;
		if (simVars->misc.compute_errors)
		{
			bool foundl = (simVars->disc.timestepping_method.find("l_")==0) || (simVars->disc.timestepping_method.find("_l_")!=std::string::npos);
			bool foundn = (simVars->disc.timestepping_method.find("n_")==0) || (simVars->disc.timestepping_method.find("_n_")!=std::string::npos);
			bool foundnl = (simVars->disc.timestepping_method.find("ln_")==0) || (foundl && foundn);
		
			if (foundnl)
				analytic_solution = 1;
			else if (foundl)
				analytic_solution = 2;
			else
				SWEETError("Computing errors for this timestepping-method is not possible");
		}



		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData_Physical u = i_planeData_u.toPhys();
		PlaneData_Physical v = i_planeData_v.toPhys();

		///// Analytical solution at current time on original grid
		///PlaneData_Spectral ts_u = t0_prog_u;
		///PlaneData_Spectral ts_v = t0_prog_v;

		PlaneData_Spectral ts_u(planeDataConfig);
		PlaneData_Spectral ts_v(planeDataConfig);
		PlaneData_Physical ts_u_phys(planeDataConfig);
		PlaneData_Physical ts_v_phys(planeDataConfig);

		if (simVars->misc.compute_errors)
		{
			//if (simVars.setup.benchmark_id > 51 && simVars.setup.benchmark_id < 65)
			if (simVars->disc.timestepping_method.find("forcing")!=std::string::npos)
			{
				if (simVars->disc.space_grid_use_c_staggering)
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				else
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars->disc.space_res_physical[0])*simVars->sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars->disc.space_res_physical[1])*simVars->sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_u(*simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				ts_u.loadPlaneDataPhysical(ts_u_phys);
				ts_v.loadPlaneDataPhysical(ts_v_phys);
			}
			else //if (simVars.setup.benchmark_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppersFine->ln_cole_hopf->run_timestep(
						 ts_u, ts_v,
						 //ts_u, ts_v,
						 simVars->timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppersFine->l_direct->run_timestep(
						 ts_u, ts_v,
						 //ts_u, ts_v,
						 simVars->timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).toPhys().physical_reduce_rms();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).toPhys().physical_reduce_rms();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).toPhys().physical_reduce_max_abs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).toPhys().physical_reduce_max_abs();

			return ts_u;
		}
		return nullptr;
	}

#endif


#if SWEET_XBRAID_SCALAR
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_xbraid_scalar(
			const double &i_u,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t,
			bool output_initial_data = false
		)
	{
		char buffer[1024];


		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(16);

		file << "#SWEET" << std::endl;
		file << "#FORMAT ASCII" << std::endl;
		file << "#PRIMITIVE SCALAR" << std::endl;

		file << i_u;

		file.close();

		return buffer;
	}
#endif

#if SWEET_XBRAID_SPHERE
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_csv_xbraid_sphere(
			const SphereData_Spectral &i_sphereData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			bool output_initial_data = false,
			bool i_phi_shifted = false
		)
	{
		char buffer[1024];

		// create copy
		SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;

	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin_xbraid_sphere(
			const SphereData_Spectral &i_sphereData,
			const char* i_name,
			int iteration_id,
			bool output_initial_data = false
	)
	{
		char buffer[1024];

		SphereData_Spectral sphereData(i_sphereData);
		//const char* filename_template = simVars.iodata.output_file_name.c_str();
		const char* filename_template = "output_%s_t%020.8f_iter%03d.sweet";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t * simVars->iodata.output_time_scale, iteration_id);
		//sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}
#endif


#if SWEET_XBRAID_PLANE
	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_xbraid_plane(
			const PlaneData_Physical &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t,
			bool output_initial_data = false
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_xbraid_plane(
			const PlaneData_Spectral &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t,
			bool output_initial_data = false
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		i_planeData.file_spectral_abs_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */
	std::string write_file_spec_amp_phase_xbraid_plane(
			const PlaneData_Spectral &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id,
			double t,
			bool output_initial_data = false
		)
	{

		char buffer[1024];

		const char* filename_template = "output_%s_amp_phase_t%020.8f_iter%03d.csv";
		if (output_initial_data)
			sprintf(buffer, filename_template, i_name, t, iteration_id);
		else
			sprintf(buffer, filename_template, i_name, t, iteration_id);

		std::ofstream file(buffer, std::ios_base::trunc);
		file << std::setprecision(12);

		for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
		{
			file << x << ", " << i_planeData.spectral_return_amplitude(0,x) << ", " << i_planeData.spectral_return_phase(0,x) << std::endl;
		}
		file.close();
		file.clear();

		return buffer;
	}

#endif





};








