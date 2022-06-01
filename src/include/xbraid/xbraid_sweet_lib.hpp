/* --------------------------------------------------------------------
 * Simulation manager structure.  
 * Holds the needed simulation data structures, e.g., discretizaion 'e' 
 * or 'i', spatial distribution, and solver used at each time point
 *
 *   comm                communicator 
 *   dim_x               spatial dimension
 *   K                   diffusion coefficient
 *   nlx                 local problem size in x-dimension
 *   nly                 local problem size in y-dimension
 *   nx                  global problem size in x-dimension
 *   ny                  global problem size in y-dimension
 *   tstart              initial time
 *   tstop               integration stop time
 *   nt                  number of time steps
 *   dx                  spatial step size in x-direction
 *   dy                  spatial step size in y-direction
 *   dt                  time step on finest grid
 *   forcing             consider non-zero forcing term
 *   vartype             variable type of of the structured spatial grid 
 *   grid_x              spatial grid
 *   stencil             discretization stencil object
 *   graph               graph object that determine the non-zero structure
 *                       of the discretization matrices
 *   A                   discretization matrix
 *   px                  number of processors in x-dimension
 *   py                  number of processors in y-dimension
 *   pi                  x-coordinate of position in processor grid
 *   pj                  y-coordinate of position in processor grid
 *   ilower              (dim_x)-dimensional array with integer indices of 
 *                       local space interval lower bounds
 *   iupper              (dim_x)-dimensional array with integer indices of
 *                       local space interval upper bounds
 *   object_type         object type of vector to access different hypre solvers 
 *   solver              hypre solver (for implicit time stepping)
 *   max_iter            maximum number of spatial MG iterations
 *   tol                 stopping tolerance for spatial MG
 *   explicit            use explicit discretization (1) or not (0)
 *   output_vis          save the error for GLVis visualization
 *   output_files        save the solution/error/error norm to files
 */
typedef struct _simulation_manager_struct
{
	MPI_Comm		comm;
	SimulationVariables*	simVars;
} simulation_manager;

/* --------------------------------------------------------------------
 * XBraid app struct 
 * -------------------------------------------------------------------- */
typedef struct _braid_App_struct {
	MPI_Comm			mpi_comm;		/* global communicator */
	MPI_Comm			mpi_comm_t;		/* communicator for parallelizing in time  */
	MPI_Comm			mpi_comm_x;		/* communicator for parallelizing in space  */
	int				mpi_pt;			/* number of processors in time  */
	simulation_manager*		man;			/* user's simulation manager structure --> defined in ex03-lib */
	SWEPlaneTimeSteppers**		timeSteppers;		/* nA sized array of solvers (one per time level) */
	int				use_rand;		/* binary, use random initial guess (1) or zero initial guess (0) */
} sweet_App;

int print_app(my_App * app)
{
	int myid,i;
	MPI_Comm_rank( app->comm, &myid );
	printf("\n\nmyid:  %d,  App contents:\n", myid);
	printf("myid:  %d,  pt:            %d\n", myid, app->pt);
	printf("myid:  %d,  use_rand:      %d\n", myid, app->use_rand);
	printf("\nmyid:  %d,  Note that some object members like comm, comm_t, comm_x, man, A and solver cannot be printed\n\n", myid);
	return 0;
}

/* --------------------------------------------------------------------
 * XBraid vector 
 * Stores the state of the simulation for a given time step
 * -------------------------------------------------------------------- */
typedef struct _braid_Vector_struct
{
	Parareal_GenericData	data;
} sweet_Vector;

/* --------------------------------------------------------------------
 * Time integrator routine that performs the update
 *   u_i = Phi_i(u_{i-1}) + g_i 
 * 
 * When Phi is called, u is u_{i-1}.
 * The return value is that u is set to u_i upon completion
 * -------------------------------------------------------------------- */
int
sweet_Step(
		sweet_App		i_app,
		sweet_Vector		i_ustop,
		sweet_Vector		i_fstop,
		sweet_Vector		io_U,
		braid_StepStatus	io_status
		)
{
	double tstart;             /* current time */
	double tstop;              /* evolve u to this time*/
	HYPRE_SStructVector  bstop;
	int level;
	int iters_taken = -1;

	/* Grab status of current time step */
	braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
	braid_StepStatusGetLevel(status, &level);

	/* Now, set up the discretization matrix.  Use the XBraid level to index
	* into the matrix lookup table */

	/* We need to "trick" the user's manager with the new dt */
	app->man->dt = tstop - tstart;

///////	/* Set up a new matrix */
///////	if( app->dt_A[level] == -1.0 )
///////	{
///////		app->nA++;
///////		app->dt_A[level] = tstop-tstart;
///////
///////		setUpImplicitMatrix( app->man );
///////		app->A[level] = app->man->A;
///////
///////		/* Set up the PFMG solver using u->x as dummy vectors. */
///////		setUpStructSolver( app->man, u->x, u->x );
///////		app->solver[level] = app->man->solver;
///////	}


	/* Time integration to next time point: Solve the system Ax = b.
	 * First, "trick" the user's manager with the right matrix and solver */ 
	///app->man->A = app->A[level];
	///app->man->solver = app->solver[level];


/////	/* Use level specific max_iter */
/////	if( level == 0 )
/////	   app->man->max_iter = app->max_iter_x[0];
/////	else
/////	   app->man->max_iter = app->max_iter_x[1];
/////
/////	/* Take step */
/////	if (fstop == NULL)
/////	{
/////	   bstop = NULL;
/////	}
/////	else
/////	{
/////	   bstop = fstop->x;
/////	}
/////	take_step(app->man, ustop->x, bstop, u->x, tstart, tstop, &iters_taken);
	timeSteppers.master->run_timestep(
						io_U,
						simVars.timecontrol.current_timestep_size,
						simVars.timecontrol.current_simulation_time
	);


//////	/* Store iterations taken */
//////	app->runtime_max_iter[level] = max_i( (app->runtime_max_iter[level]), iters_taken);

	/* Tell XBraid no refinement */
	braid_StepStatusSetRFactor(status, 1);

	return 0;
}

/* --------------------------------------------------------------------
 * -------------------------------------------------------------------- */
int
sweet_Residual(
			sweet_App		i_app,
			sweet_Vector		i_ustop,
			sweet_Vector		o_r,
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

	/* We need to "trick" the user's manager with the new dt */
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

	/* Compute residual Ax */
	app->man->A = app->A[level];
	comp_res(app->man, ustop->x, r->x, tstart, tstop);

	return 0;
}

/* --------------------------------------------------------------------
 * Create a vector object for a given time point.
 * This function is only called on the finest level.
 * -------------------------------------------------------------------- */
int
sweet_Init(
		sweet_App	i_app,
		double		i_t,
		sweet_Vector*	o_U
		)
{

	sweet_Vector * U = (sweet_Vector *) malloc( sizeof(sweet_Vector) );

	if( t == app->man->tstart )
	{
		/* Sets u_ptr as the initial condition */
		swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, *simVars, *op_plane);
		U->data = func(t0_prog_h_pert, t0_prog_u, t0_prog_v);
	}
	else if (app->use_rand)
	{
		/* This t-value will tell set_initial_condition() below to make u_ptr uniformly random */
		U->data = func_rand();
	}
	else
	{
		/* Sets u_ptr as an all zero vector*/
		U->data = func_zero();
	}

	*o_U = U;
	return 0;
}

/* --------------------------------------------------------------------
 * Create a copy of a vector object.
 * -------------------------------------------------------------------- */
int
sweet_Clone(
		sweet_App	i_app,
		sweet_Vector	i_U,
		sweet_Vector*	o_V
	)
{
	sweet_Vector *V = (sweet_Vector *) malloc(sizeof(sweet_Vector));

	V->data = U->data;

	*o_V = V;
	return 0;
}

/* --------------------------------------------------------------------
 * Destroy vector object.
 * -------------------------------------------------------------------- */
int
sweet_Free(
		sweet_App	i_app,
		sweet_Vector	i_U)
{
	delete i_U->data;
	free(i_U);

	return 0;
}

/* --------------------------------------------------------------------
 * Compute vector sum y = alpha*x + beta*y.
 * -------------------------------------------------------------------- */
int
sweet_Sum(
		sweet_App	i_app,
		double		i_alpha,
		sweet_Vector	i_X,
		double		i_beta,
		sweet_Vectori&	io_y
	)
{

	io_Y->data = i_alpha * i_X->data + i_beta * io_Y_data;

	return 0;
}

/* --------------------------------------------------------------------
 * User access routine to spatial solution vectors and allows for user
 * output.  The default XBraid parameter of access_level=1, calls 
 * my_Access only after convergence and at every time point.
 * -------------------------------------------------------------------- */
int
sweet_Access(
			sweet_App		i_app,
			sweet_Vector		i_U,
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
int
sweet_SpatialNorm(
		sweet_App     i_app,
		sweet_Vector  i_U,
		double*       o_norm)
{
	*o_norm = i_U->data->reduceNorm();
	return 0;
}

/* --------------------------------------------------------------------
 * Return buffer size needed to pack one spatial sweet_Vector.  Here the
 * vector contains one double at every grid point and thus, the buffer 
 * size is the number of grid points.
 * -------------------------------------------------------------------- */
int
sweet_BufSize(
		sweet_App		i_app,
		int*			o_size,
		braid_BufferStatus	o_status)
{
	*o_size = i_app->size_buffer * sizeof(double);
	return 0;
}

/* --------------------------------------------------------------------
 * Pack a sweet_Vector into a buffer.
 * -------------------------------------------------------------------- */
int
sweet_BufPack(
		sweet_App		i_app,
		sweet_Vector		i_u,
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
 * Unpack a buffer and place into a sweet_Vector
 * -------------------------------------------------------------------- */
int
sweet_BufUnpack(
		sweet_App		i_app,
		void*			i_buffer,
		sweet_Vector*		o_U,
		braid_BufferStatus	status
		)
{
	double* dbuffer = buffer;
	sweet_Vector* U = (sweet_Vector *) malloc( sizeof(sweet_Vector) );

	U->data->deserialize(dbuffer);

	o_U = U;

	return 0;
}
