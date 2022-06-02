#include <braid.hpp>

#include <parareal/Parareal_GenericData.hpp>
#if SWEET_XBRAID_SCALAR
	#include <parareal/Parareal_GenericData_Scalar.hpp>
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


typedef struct _braid_Vector_struct sweet_BraidVector;
class sweet_BraidApp;



/* --------------------------------------------------------------------
 * XBraid vector 
 * Stores the state of the simulation for a given time step
 * -------------------------------------------------------------------- */
typedef struct _braid_Vector_struct
{
	Parareal_GenericData*	data = nullptr;

	sweet_BraidVector()
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

	sweet_BraidVector(sweet_BraidVector &i_vector)
	{
		this->data = i_vector.data;
	};

	sweet_BraidVector& operator=(const sweet_BraidVector &i_vector)
	{
		this->data = i_vector.data;
		return *this;
	};

	sweet_BraidVector operator+(
			const sweet_BraidVector &i_vector
	)	const
	{
		sweet_BraidVector out;
		out->data = this->data + i_vector.data
		return out;
	}

	sweet_BraidVector operator*(
			const double i_value
	)	const
	{
		sweet_BraidVector out;
		out->data = this->data * i_value;
		return out;
	}


	void allocate_data()
	{
#if SWEET_PARAREAL_SCALAR
		{
			data = new Parareal_GenericData_Scalar<N>;
			out->allocate_data();
			out->set_time(this->timeframe_end);
			return out;
		}

#elif SWEET_PARAREAL_PLANE
		{
			data = new Parareal_GenericData_PlaneData_Spectral<N>;
			out->setup_data_config(this->planeDataConfig);
			out->allocate_data();
			return out;
		}

#elif SWEET_PARAREAL_SPHERE
		{
			data = new Parareal_GenericData_SphereData_Spectral<N>;
			out->setup_data_config(this->sphereDataConfig);
			out->allocate_data();
			return out;
		}
#endif
	}



} sweet_BraidVector;


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
	return i_vector*i_value;
}


class sweet_BraidApp
			: public BraidApp
{

	/* --------------------------------------------------------------------
	 * Time integrator routine that performs the update
	 *   u_i = Phi_i(u_{i-1}) + g_i 
	 * 
	 * When Phi is called, u is u_{i-1}.
	 * The return value is that u is set to u_i upon completion
	 * -------------------------------------------------------------------- */
	braid_Int
	Step(
			sweet_BraidVector		i_ustop,
			sweet_BraidVector		i_fstop,
			sweet_BraidVector		io_U,
			braid_StepStatus	io_status
			)
	{
		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;
	
		/* Grab status of current time step */
		braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
		braid_StepStatusGetLevel(status, &level);
	
		/* Set the new dt in the user's manager*/
		app->man->dt = tstop - tstart;
	
		/* Set the correct solver for the current level */
		app->man->timeSteppers = app->timeSteppers[solver];
	
		timeSteppers.master->run_timestep(
							io_U,
							simVars.timecontrol.current_timestep_size,
							simVars.timecontrol.current_simulation_time
		);
	
		/* Tell XBraid no refinement */
		braid_StepStatusSetRFactor(status, 1);
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * -------------------------------------------------------------------- */
	braid_Int
	Residual(
				sweet_BraidVector		i_ustop,
				sweet_BraidVector		o_r,
				braid_StepStatus	io_status
			)
	{
		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;
	
		/* Grab status of current time step */
		braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
	
		/* Grab level */
		braid_StepStatusGetLevel(status, &level);
	
		/* Set the new dt in the user's manager*/
		app->man->dt = tstop - tstart;
	
		/* Now, set up the discretization matrix.  Use the XBraid level to index
		 * into the matrix lookup table */
		if( app->dt_A[level] == -1.0 ){
			app->nA++;
			app->dt_A[level] = tstop-tstart;
	
			setUpImplicitMatrix( app->man );
			app->A[level] = app->man->A;
	
			/* Set up the PFMG solver using r->x as dummy vectors. */
			setUpStructSolver( app->man, r->x, r->x );
			app->solver[level] = app->man->solver;
		}
	
		/* Set the correct solver for the current level */
		app->man->timeSteppers = app->timeSteppers[solver];
	
		/* Compute residual Ax */
		//////app->man->A = app->A[level];
		comp_res(app->man, ustop->x, r->x, tstart, tstop);
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a vector object for a given time point.
	 * This function is only called on the finest level.
	 * -------------------------------------------------------------------- */
	braid_Int
	Init(
			double		i_t,
			sweet_BraidVector*	o_U
			)
	{
	
		///sweet_BraidVector * U = (sweet_BraidVector *) malloc( sizeof(sweet_BraidVector) );
	
		if( t == app->man->tstart )
		{
	#if SWEET_XBRAID_SCALAR
			double u0 = atof(simVars->bogus.var[1].c_str());
			this->dataArrays_to_GenericData_Scalar(o_U->data, u0);
	
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
	
	
			this->dataArrays_to_GenericData_PlaneData_Spectral(	U->data,
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
	
			this->dataArrays_to_GenericData_SphereData_Spectral(o_U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
	
	
		}
		else if (app->use_rand)
		{
	#if SWEET_XBRAID_SCALAR
			double u0 = ((double)braid_Rand())/braid_RAND_MAX;
			this->dataArrays_to_GenericData_Scalar(o_U->data, u0);
	#elif SWEET_XBRAID_PLANE
	
		#if SWEET_XBRAID_PLANE_SWE
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			PlaneData_Physical t0_prog_h_phys(t0_prog_u.planeDataConfig);
		#endif
	
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			PlaneData_Spectral t0_prog_v(planeDataConfig);
	
			PlaneData_Physical t0_prog_u_phys(t0_prog_u.planeDataConfig);
			PlaneData_Physical t0_prog_v_phys(t0_prog_v.planeDataConfig);
	
		#if SWEET_XBRAID_PLANE_SWE
			t0_prog_h_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = i_app->man->simVars->sim.h0 + ((double)braid_Rand())/braid_RAND_MAX;
				}
		#endif
			t0_prog_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			t0_prog_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
	
		#if SWEET_XBRAID_PLANE_SWE
			t0_prog_h.loadPlaneDataPhysical(t0_prog_u_phys);
		#endif
			t0_prog_u.loadPlaneDataPhysical(t0_prog_u_phys);
			t0_prog_v.loadPlaneDataPhysical(t0_prog_v_phys);
	
			this->dataArrays_to_GenericData_PlaneData_Spectral(	o_U->data,
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
					io_data = i_app->man->simVars->sim.h0 + ((double)braid_Rand())/braid_RAND_MAX;
				}
			t0_prog_vrt_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
			t0_prog_div_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
				{
					io_data = ((double)braid_Rand())/braid_RAND_MAX;
				}
	
			t0_prog_phi_pert.loadPlaneDataPhysical(t0_prog_phi_pert_phys);
			t0_prog_vrt.loadPlaneDataPhysical(t0_prog_vrt_phys);
			t0_prog_div.loadPlaneDataPhysical(t0_prog_div_phys);
	
			this->dataArrays_to_GenericData_SphereData_Spectral(o_U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
		}
		else
		{
			/* Sets U as an all zero vector*/
	#if SWEET_XBRAID_SCALAAR
			this->dataArrays_to_GenericData_Scalar(U->data, 0.);
	#elif SWEET_XBRAID_PLANE
		#if SWEET_XBRAID_PLANE_SWE
			PlaneData_Spectral t0_prog_h_pert(planeDataConfig);
			t0_prog_h_pert.spectral_set_zero();
		#endif
	
			PlaneData_Spectral t0_prog_u(planeDataConfig);
			t0_prog_u.spectral_set_zero();
	
			PlaneData_Physical t0_prog_u_phys(t0_prog_u.planeDataConfig);
			t0_prog_v.spectral_set_zero();
	
			this->dataArrays_to_GenericData_PlaneData_Spectral(	U->data,
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
	
			this->dataArrays_to_GenericData_SphereData_Spectral(o_U->data, t0_prog_phi_pert, t0_prog_vrt, t0_prog_div);
	#endif
		}
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Create a copy of a vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Clone(
			sweet_BraidVector	i_U,
			sweet_BraidVector*	o_V
		)
	{
		o_V->data = U->data;
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Destroy vector object.
	 * -------------------------------------------------------------------- */
	braid_Int
	Free(
			sweet_BraidVector	i_U)
	{
		delete i_U->data;
		///free(i_U);
	        delete i_U;
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Compute vector sum y = alpha*x + beta*y.
	 * -------------------------------------------------------------------- */
	braid_Int
	Sum(
			double		i_alpha,
			sweet_BraidVector	i_X,
			double		i_beta,
			sweet_BraidVector&	io_Y
		)
	{
	
		io_Y = i_X * i_alpha + io_Y * i_beta;
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * User access routine to spatial solution vectors and allows for user
	 * output.  The default XBraid parameter of access_level=1, calls 
	 * my_Access only after convergence and at every time point.
	 * -------------------------------------------------------------------- */
	braid_Int
	Access(
				sweet_BraidVector		i_U,
				braid_AccessStatus	io_astatus
			)
	{
		double     tstart         = (app->man->tstart);
		double     tstop          = (app->man->tstop);
		int        nt             = (app->man->nt);
	
		double     rnorm, disc_err, t;
		int        iter, level, done, index, myid;
		char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];
	
		/* Retrieve current time from Status Object */
		braid_AccessStatusGetT(astatus, &t);
	
		/* Retrieve XBraid State Information from Status Object */
		MPI_Comm_rank(app->comm_x, &myid);
		braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
		braid_AccessStatusGetResidual(astatus, &rnorm);
	
		if(level == 0)
		{
	/////   /* Print discretization error to screen for only final time */
	/////   index = ((t - tstart) / ((tstop - tstart)/nt) + 0.1);
	/////   compute_disc_err(app->man, u->x, t, app->e, &disc_err);
	/////   if( (t == app->man->tstop) && myid == 0 ) {
	/////      printf("\n  Discr. error         = %1.5e\n", disc_err);
	/////      printf("\n  my_Access():  Braid iter %d,  discr. error at final time:  %1.4e\n", iter, disc_err);
	
		}
	
		/* Write the norm of the discretization error to a separate file for each time step */
		if( app->man->output_files )
		{
	/////         sprintf(filename, "%s.iter%03d.time%07d", "ex-03.error_norm", iter, index);
	/////         output_error_file(app->man, t, disc_err, filename); 
		}
	
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Compute norm of a spatial vector 
	 * -------------------------------------------------------------------- */
	braid_Int
	SpatialNorm(
			sweet_BraidVector  i_U,
			double*       o_norm)
	{
		*o_norm = i_U->data->reduce_norm2();
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Return buffer size needed to pack one spatial sweet_BraidVector.  Here the
	 * vector contains one double at every grid point and thus, the buffer 
	 * size is the number of grid points.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufSize(
			int*			o_size,
			braid_BufferStatus	o_status)
	{
		*o_size = i_app->size_buffer * sizeof(double);
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Pack a sweet_BraidVector into a buffer.
	 * -------------------------------------------------------------------- */
	braid_Int
	BufPack(
			sweet_BraidVector		i_u,
			void*			o_buffer,
			braid_BufferStatus	o_status
		)
	{
	
		double* dbuffer = buffer;
	
		U->data->serialize(dbuffer);
	
		/* Return the number of bytes actually packed */
		braid_BufferStatusSetSize(o_status, i_app->size_buffer*sizeof(double) );
		return 0;
	}
	
	/* --------------------------------------------------------------------
	 * Unpack a buffer and place into a sweet_BraidVector
	 * -------------------------------------------------------------------- */
	braid_Int
	BufUnpack(
			void*			i_buffer,
			sweet_BraidVector*		o_U,
			braid_BufferStatus	status
			)
	{
		double* dbuffer = buffer;
		///sweet_BraidVector* U = (sweet_BraidVector *) malloc( sizeof(sweet_BraidVector) );
		sweet_BraidVector* U = new sweet_BraidVector;
	
		U->data->deserialize(dbuffer);
	
		o_U = U;
	
		return 0;
	}

};




/* --------------------------------------------------------------------
 * Simulation manager structure.  
 * Holds the needed simulation data structures, e.g., discretizaion 'e' 
 * or 'i', spatial distribution, and solver used at each time point
 */
///template <class t_tsmType, int N>
typedef struct _simulation_manager_struct
{
	MPI_Comm		comm;
	SimulationVariables*	simVars;
	t_tsmType*		timeSteppers;
} simulation_manager;

//////////////////////* --------------------------------------------------------------------
///////////////////// * XBraid app struct 
///////////////////// * -------------------------------------------------------------------- */
////////////////////////template <class t_tsmType, int N>
/////////////////////typedef struct _braid_App_struct {
/////////////////////	MPI_Comm			mpi_comm;		/* global communicator */
/////////////////////	MPI_Comm			mpi_comm_t;		/* communicator for parallelizing in time  */
/////////////////////	MPI_Comm			mpi_comm_x;		/* communicator for parallelizing in space  */
/////////////////////	int				mpi_pt;			/* number of processors in time  */
/////////////////////	simulation_manager*		man;			/* user's simulation manager structure --> defined in ex03-lib */
/////////////////////	t_tsmType**			timeSteppers;		/* nA sized array of solvers (one per time level) */
/////////////////////	int				use_rand;		/* binary, use random initial guess (1) or zero initial guess (0) */
/////////////////////
/////////////////////#if SWEET_XBRAID_PLANE
/////////////////////	PlaneDataConfig*		planeDataConfig;
/////////////////////#elif SWEET_XBRAID_SPHERE
/////////////////////	PlaneDataConfig*		sphereDataConfig;
/////////////////////#endif
/////////////////////
/////////////////////} sweet_App;

int print_app(sweet_BraidApp* app)
{
	int myid,i;
	MPI_Comm_rank( app->comm, &myid );
	printf("\n\nmyid:  %d,  App contents:\n", myid);
	printf("myid:  %d,  pt:            %d\n", myid, app->pt);
	printf("myid:  %d,  use_rand:      %d\n", myid, app->use_rand);
	printf("\nmyid:  %d,  Note that some object members like comm, comm_t, comm_x, man, A and solver cannot be printed\n\n", myid);
	return 0;
}


