#include "_hypre_utilities.h"
#include "HYPRE_sstruct_ls.h"
#include "_hypre_sstruct_mv.h"

#include "braid.h"
#include "braid_test.h"
#include "ex-03-lib.c"


/* --------------------------------------------------------------------
 * XBraid app struct 
 * -------------------------------------------------------------------- */
typedef struct _braid_App_struct {
	MPI_Comm                comm;             /* global communicator */
	MPI_Comm                comm_t;           /* communicator for parallelizing in time  */
	MPI_Comm                comm_x;           /* communicator for parallelizing in space  */
	int                     pt;               /* number of processors in time  */
	simulation_manager     *man;              /* user's simulation manager structure */
	HYPRE_SStructVector     e;                /* temporary vector used for error computations */
	int                     nA;               /* number of discr. matrices that have been created */
	int                     max_nA;           /* max nA value allowed */
	HYPRE_SStructMatrix    *A;                /* nA sized array of discr. matrices (one per time level) */
	double                 *dt_A;             /* nA sized array of time step sizes for each  matrix  */
	HYPRE_StructSolver     *solver;           /* nA sized array of solvers (one per time level) */
	int                     use_rand;         /* binary, use random initial guess (1) or zero initial guess (0) */
	int                    *runtime_max_iter; /* runtime info for the max number of spatial solve iterations at each level */
	int                    *max_iter_x;       /* length 2 array of expensive and cheap max PFMG iters for spatial solves*/
} swe_plane_App;

int print_app(my_App * app)
{
	int myid,i;
	MPI_Comm_rank( app->comm, &myid );
	printf("\n\nmyid:  %d,  App contents:\n", myid);
	printf("myid:  %d,  pt:            %d\n", myid, app->pt);
	printf("myid:  %d,  use_rand:      %d\n", myid, app->use_rand);
	printf("myid:  %d,  nA:            %d\n", myid, app->nA);
	printf("myid:  %d,  max_iter_x[0]: %d\n", myid, app->max_iter_x[0]);
	printf("myid:  %d,  max_iter_x[1]: %d\n", myid, app->max_iter_x[1]);
	for(i = 0; i < app->nA; i++){
		printf("myid:  %d,  runtime_max_iter[%d]: %d\n", myid, i, app->runtime_max_iter[i]);
	}
	for(i = 0; i < app->nA; i++){
		printf("myid:  %d,  dt_A[%d]:           %1.2e\n", myid, i, app->dt_A[i]);
	}
	printf("\nmyid:  %d,  Note that some object members like comm, comm_t, comm_x, man, A and solver cannot be printed\n\n", myid);
	return 0;
}

/* --------------------------------------------------------------------
 * XBraid vector 
 * Stores the state of the simulation for a given time step
 * -------------------------------------------------------------------- */
typedef struct _braid_Vector_struct {
	HYPRE_SStructVector   x;
} swe_plane_Vector;

/* --------------------------------------------------------------------
 * Time integrator routine that performs the update
 *   u_i = Phi_i(u_{i-1}) + g_i 
 * 
 * When Phi is called, u is u_{i-1}.
 * The return value is that u is set to u_i upon completion
 * -------------------------------------------------------------------- */
int
swe_plane_Step(	braid_App        app,
		braid_Vector     ustop,
		braid_Vector     fstop,
		braid_Vector     u,
		braid_StepStatus status)
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

   /* Set up a new matrix */
   if( app->dt_A[level] == -1.0 ){
      app->nA++;
      app->dt_A[level] = tstop-tstart;

      setUpImplicitMatrix( app->man );
      app->A[level] = app->man->A;
      
      /* Set up the PFMG solver using u->x as dummy vectors. */
      setUpStructSolver( app->man, u->x, u->x );
      app->solver[level] = app->man->solver;
   } 

   /* Time integration to next time point: Solve the system Ax = b.
    * First, "trick" the user's manager with the right matrix and solver */ 
   app->man->A = app->A[level];
   app->man->solver = app->solver[level];

   /* Use level specific max_iter */
   if( level == 0 )
      app->man->max_iter = app->max_iter_x[0];
   else
      app->man->max_iter = app->max_iter_x[1];

   /* Take step */
   if (fstop == NULL)
   {
      bstop = NULL;
   }
   else
   {
      bstop = fstop->x;
   }
   take_step(app->man, ustop->x, bstop, u->x, tstart, tstop, &iters_taken);

   /* Store iterations taken */
   app->runtime_max_iter[level] = max_i( (app->runtime_max_iter[level]), iters_taken);

   /* Tell XBraid no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

/* --------------------------------------------------------------------
 * -------------------------------------------------------------------- */
int
swe_plane_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
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
swe_plane_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   
   my_Vector * u = (my_Vector *) malloc( sizeof(my_Vector) );
   
   if( t == app->man->tstart ){
      /* Sets u_ptr as the initial condition */
      t = 0.0;
   }
   else if (app->use_rand){
      /* This t-value will tell set_initial_condition() below to make u_ptr uniformly random */
      t = -1.0;
   }
   else{
      /* Sets u_ptr as an all zero vector*/
      t = 1.0;
   }

   set_initial_condition(app->man, &(u->x), t);
   *u_ptr = u;
   return 0;
}

/* --------------------------------------------------------------------
 * Create a copy of a vector object.
 * -------------------------------------------------------------------- */
int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v = (my_Vector *) malloc(sizeof(my_Vector));
   double    *values;
   initialize_vector(app->man, &(v->x));

   /* Set the values. */
   values = (double *) malloc( (app->man->nlx)*(app->man->nly)*sizeof(double) );
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, 0, app->man->ilower, app->man->iupper, 0, values );
   HYPRE_SStructVectorSetBoxValues( v->x, 0, app->man->ilower, app->man->iupper, 0, values );
   free( values );
   HYPRE_SStructVectorAssemble( v->x );

   *v_ptr = v;
   return 0;
}

/* --------------------------------------------------------------------
 * Destroy vector object.
 * -------------------------------------------------------------------- */
int
swe_plane_Free(braid_App    app,
        braid_Vector u)
{
   HYPRE_SStructVectorDestroy( u->x );
   free( u );

   return 0;
}

/* --------------------------------------------------------------------
 * Compute vector sum y = alpha*x + beta*y.
 * -------------------------------------------------------------------- */
int
swe_plane_Sum(braid_App    app,
       double       alpha,
       braid_Vector x,
       double       beta,
       braid_Vector y)
{
   int i;
   double *values_x, *values_y;
   
   values_x = (double *) malloc( (app->man->nlx)*(app->man->nly)*sizeof(double) );
   values_y = (double *) malloc( (app->man->nlx)*(app->man->nly)*sizeof(double) );

   HYPRE_SStructVectorGather( x->x );
   HYPRE_SStructVectorGetBoxValues( x->x, 0, (app->man->ilower), (app->man->iupper), 0, values_x );
   HYPRE_SStructVectorGather( y->x );
   HYPRE_SStructVectorGetBoxValues( y->x, 0, (app->man->ilower), (app->man->iupper), 0, values_y );

   for( i = 0; i < (app->man->nlx)*(app->man->nly); i++ ){
      values_y[i] = alpha*values_x[i] + beta*values_y[i];
   }

   HYPRE_SStructVectorSetBoxValues( y->x, 0, (app->man->ilower), (app->man->iupper), 0, values_y );

   free( values_x );
   free( values_y );
   return 0;
}

/* --------------------------------------------------------------------
 * User access routine to spatial solution vectors and allows for user
 * output.  The default XBraid parameter of access_level=1, calls 
 * my_Access only after convergence and at every time point.
 * -------------------------------------------------------------------- */
int
swe_plane_Access(braid_App           app,
          braid_Vector        u,
          braid_AccessStatus  astatus)
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
      /* Print discretization error to screen for only final time */
      index = ((t - tstart) / ((tstop - tstart)/nt) + 0.1);
      compute_disc_err(app->man, u->x, t, app->e, &disc_err);
      if( (t == app->man->tstop) && myid == 0 ) {
         printf("\n  Discr. error         = %1.5e\n", disc_err);
         printf("\n  my_Access():  Braid iter %d,  discr. error at final time:  %1.4e\n", iter, disc_err);

      }
      
      /* Write the norm of the discretization error to a separate file for each time step */
      if( app->man->output_files ){
         sprintf(filename, "%s.iter%03d.time%07d", "ex-03.error_norm", iter, index);
         output_error_file(app->man, t, disc_err, filename); 
      }
   }

   /* Write THREE GLVIS visualization files for the final time step:
    * (1) the discretization error (2) the true solution (3) the discrete solution */ 
   if( app->man->output_vis && (level == 0) && (t == app->man->tstop) )
   {
      sprintf(filename_mesh, "%s.iter%03d", "ex-03_mesh", iter);
      sprintf(filename_err, "%s.iter%03d", "ex-03_err_tstop", iter);
      sprintf(filename_sol, "%s.iter%03d", "ex-03_sol_tstop", iter);
      output_vis(app->man, u->x, t, filename_mesh, filename_err, filename_sol);
   }
   
   return 0;
}

/* --------------------------------------------------------------------
 * Compute norm of a spatial vector 
 * -------------------------------------------------------------------- */
int
swe_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   norm(u->x, norm_ptr);
   return 0;
}

/* --------------------------------------------------------------------
 * Return buffer size needed to pack one spatial braid_Vector.  Here the
 * vector contains one double at every grid point and thus, the buffer 
 * size is the number of grid points.
 * -------------------------------------------------------------------- */
int
swe_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  status)
{
    *size_ptr = (app->man->nlx)*(app->man->nly)*sizeof(double);
    return 0;
}

/* --------------------------------------------------------------------
 * Pack a braid_Vector into a buffer.
 * -------------------------------------------------------------------- */
int
swe_BufPack(braid_App           app,
           braid_Vector        u,
           void                *buffer,
           braid_BufferStatus  status)
{
   double *dbuffer = buffer;
   
   /* Place the values in u into the buffer */
   HYPRE_SStructVectorGather( u->x );
   HYPRE_SStructVectorGetBoxValues( u->x, 0, app->man->ilower, 
                          app->man->iupper, 0, &(dbuffer[0]) );

   /* Return the number of bytes actually packed */
   braid_BufferStatusSetSize( status, (app->man->nlx)*(app->man->nly)*sizeof(double) );
   return 0;
}

/* --------------------------------------------------------------------
 * Unpack a buffer and place into a braid_Vector
 * -------------------------------------------------------------------- */
int
swe_BufUnpack(braid_App           app,
             void                *buffer,
             braid_Vector        *u_ptr,
             braid_BufferStatus  status)
{
   double    *dbuffer = buffer;
   my_Vector *u       = (my_Vector *) malloc( sizeof(my_Vector) );
   
   /* Set the values in u based on the values in the buffer */
   initialize_vector(app->man, &(u->x));
   HYPRE_SStructVectorSetBoxValues( u->x, 0, app->man->ilower, 
                           app->man->iupper, 0, &(dbuffer[0]) );
   HYPRE_SStructVectorAssemble( u->x );
   *u_ptr = u;

   return 0;
}
