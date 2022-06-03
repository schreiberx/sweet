// TODO:
// in all methods in App:
// sweet_BraidVector *u = (sweet_BraidVector*) u_; where u_ is braid_Vector


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


class sweet_BraidVector;
class sweet_BraidApp;


///////////////////////* --------------------------------------------------------------------
////////////////////// * Simulation manager structure.  
////////////////////// * Holds the needed simulation data structures, e.g., discretizaion 'e' 
////////////////////// * or 'i', spatial distribution, and solver used at each time point
////////////////////// */
/////////////////////////template <class t_tsmType, int N>
//////////////////////typedef struct _simulation_manager_struct
//////////////////////{
//////////////////////	MPI_Comm		comm;
//////////////////////	SimulationVariables*	simVars;
//////////////////////	t_tsmType*		timeSteppers;
//////////////////////	double			dt;
//////////////////////} simulation_manager;


////////////int print_app(sweet_BraidApp* app)
////////////{
////////////	int myid,i;
////////////	MPI_Comm_rank( app->comm, &myid );
////////////	printf("\n\nmyid:  %d,  App contents:\n", myid);
////////////	printf("myid:  %d,  pt:            %d\n", myid, app->pt);
////////////	printf("myid:  %d,  use_rand:      %d\n", myid, app->use_rand);
////////////	printf("\nmyid:  %d,  Note that some object members like comm, comm_t, comm_x, man, A and solver cannot be printed\n\n", myid);
////////////	return 0;
////////////}




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

	sweet_BraidVector(const sweet_BraidVector &i_vector)
	{
		this->data = i_vector.data;
	};

	sweet_BraidVector& operator=(const sweet_BraidVector &i_vector)
	{
		*this->data = *i_vector.data;
		return *this;
	};

	sweet_BraidVector operator+(
			const sweet_BraidVector &i_vector
	)	const
	{
		sweet_BraidVector out;
		*out.data = *this->data;
		*out.data += *i_vector.data;
		return out;
	}

	sweet_BraidVector operator*(
			const double i_value
	)	const
	{
		sweet_BraidVector out;
		*out.data = *this->data;
		*out.data *= i_value;
		return out;
	}


	void allocate_data()
	{
#if SWEET_PARAREAL_SCALAR
		{
			data = new Parareal_GenericData_Scalar<N>;
			//this->set_time(this->timeframe_end);
		}

#elif SWEET_PARAREAL_PLANE
		{
			data = new Parareal_GenericData_PlaneData_Spectral<N>;
			//this->setup_data_config(this->planeDataConfig);
			//this->set_time(this->timeframe_end);
		}

#elif SWEET_PARAREAL_SPHERE
		{
			data = new Parareal_GenericData_SphereData_Spectral<N>;
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

////////////	simulation_manager man;

	MPI_Comm		comm;
	SimulationVariables*	simVars;
	double			dt;
	t_tsmType**		timeSteppers;

#if SWEET_XBRAID_PLANE
	PlaneOperators*		op_plane;
#elif SWEET_XBRAID_SPHERE
	SphereOperators*	op_sphere;
#endif


	double			tstart;
	double			tstop;

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
			)
		: BraidApp(i_comm_t, i_tstart, i_tstop, i_ntime)
	{
		this->rank = i_rank;
		this->simVars = i_simVars;
	}

	~sweet_BraidApp()
	{
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
		braid_StepStatus& status = (braid_StepStatus&) io_status;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;

		/* Grab status of current time step */
		braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
		braid_StepStatusGetLevel(status, &level);

		/* Set the new dt in the user's manager*/
		this->dt = tstop - tstart;

		this->timeSteppers[level]->master->run_timestep(
								U->data,
								this->simVars->timecontrol.current_timestep_size,
								this->simVars->timecontrol.current_simulation_time
		);
	
		/* Tell XBraid no refinement */
		braid_StepStatusSetRFactor(status, 1);

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

		braid_StepStatus& status = (braid_StepStatus&) io_status;

		double tstart;             /* current time */
		double tstop;              /* evolve u to this time*/
		int level;

		/* Grab status of current time step */
		braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
	
		/* Grab level */
		braid_StepStatusGetLevel(status, &level);
	
		/* Set the new dt in the user's manager*/
		this->dt = tstop - tstart;
	
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


		sweet_BraidVector* U = new sweet_BraidVector;
	
		if( i_t == this->tstart )
		{
	#if SWEET_XBRAID_SCALAR
			double u0 = atof(simVars->bogus.var[1].c_str());
			this->dataArrays_to_GenericData_Scalar(U->data, u0);
	
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
	#if SWEET_XBRAID_SCALAR
			double u0 = ((double)braid_Rand())/braid_RAND_MAX;
			this->dataArrays_to_GenericData_Scalar(U->data, u0);
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
			/* Sets U as an all zero vector*/
	#if SWEET_XBRAID_SCALAR
			this->dataArrays_to_GenericData_Scalar(U->data, 0.);
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
		sweet_BraidVector* V = new sweet_BraidVector;
		V->data = U->data;
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
		delete U->data;
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
		int        iter, level, done, index, myid;
		char       filename[255], filename_mesh[255], filename_err[255], filename_sol[255];

		braid_AccessStatus& astatus = (braid_AccessStatus&) io_astatus;

		/* Retrieve current time from Status Object */
		braid_AccessStatusGetT(astatus, &t);
	
		/* Retrieve XBraid State Information from Status Object */
		///////////MPI_Comm_rank(app->comm_x, &myid);
		///////////braid_AccessStatusGetTILD(astatus, &t, &iter, &level, &done);
		///////////braid_AccessStatusGetResidual(astatus, &rnorm);
	
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
	//////////	if( this->output_files )
	//////////	{
	//////////         sprintf(filename, "%s.iter%03d.time%07d", "ex-03.error_norm", iter, index);
	//////////         output_error_file(app->man, t, disc_err, filename); 
	//////////	}
	
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
		*o_size = this->size_buffer * sizeof(double);
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

		/* Return the number of bytes actually packed */
		o_status.SetSize( this->size_buffer*sizeof(double) );
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

		sweet_BraidVector* U = new sweet_BraidVector;

#if SWEET_XBRAID_SCALAR
		double* dbuffer = (double*) i_buffer;
#else
		std::complex<double>* dbuffer = (std::complex<double>*) i_buffer;
#endif

		U->data->deserialize(dbuffer);

		*o_U = (braid_Vector) U;

		return 0;
	}

};





