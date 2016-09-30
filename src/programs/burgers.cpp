/*
 * Burgers equation
 *
 */

#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/BurgersValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/Sampler2D.hpp>
#include <sweet/SemiLagrangian.hpp>
#include <programs/burgers_HelmholtzSolver.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 0
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal.hpp>
#endif

// Input parameters (cmd line)

//general parameters
SimulationVariables simVars;

//specific parameters
int param_time_scheme = -1;
bool param_compute_error = 0;
bool param_use_staggering = 0;
bool param_semilagrangian = 0;
int param_time_scheme_coarse = -1;


class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	DataArray<2> prog_u, prog_v;

	// Prognostic variables at time step t-dt
	DataArray<2> prog_u_prev, prog_v_prev;

	// temporary variables - may be overwritten, use locally
	DataArray<2> tmp;

	// Points mapping [0,simVars.sim.domain_size[0])x[0,simVars.sim.domain_size[1])
	// with resolution simVars.sim.resolution
	DataArray<2> pos_x, pos_y;

	// Arrival points for semi-lag
	DataArray<2> posx_a, posy_a;

	// Departure points for semi-lag
	DataArray<2> posx_d, posy_d;

	//Staggering displacement array (use 0.5 for each displacement)
	// [0] - delta x of u variable
	// [1] - delta y of u variable
	// [2] - delta x of v variable
	// [3] - delta y of v variable
	// Default - A grid (there is shift in x,y of 1/2 to all vars)
	// For C grid use {0,-0.5,-0.5,0}
	double stag_displacement[4] = {-0.5,-0.5,-0.5,-0.5};
	double stag_h[2] = {-0.5,-0.5};
	double stag_u[2] = {-0.5,-0.5};
	double stag_v[2] = {-0.5,-0.5};

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	// Max difference to initial conditions
	double benchmark_diff_u;
	double benchmark_diff_v;

	DataArray<2> benchmark_analytical_error;

	// Error measures L2 norm
	double benchmark_analytical_error_rms_u;
	double benchmark_analytical_error_rms_v;

	// Error measures max norm
	double benchmark_analytical_error_maxabs_u;
	double benchmark_analytical_error_maxabs_v;

	// Finite difference operators
	Operators2D op;

	// Runge-Kutta stuff
	TimesteppingRK timestepping;

	// Interpolation stuff
	Sampler2D sampler2D;

	// Semi-Lag stuff
	SemiLagrangian semiLagrangian;


	/**
	 * Two dimensional Burgers equation
	 *
	 * Equations:
	 *
	 *     \f$ u_t + uu_x + vu_y = \nu (u_xx + u_yy) \f$
	 *     \f$ v_t + uv_x + vv_y = \nu (v_xx + v_yy) \f$
	 *
	 *   ______________
	 *   |            |
	 *   |    u0,1    |
	 *   v0,0 P0,0 v1,0
	 *   |    u0,0    |
	 *   |____________|
	 */
public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup
		prog_u(simVars.disc.res),	// velocity (x-direction)
		prog_v(simVars.disc.res),	// velocity (y-direction)

		prog_u_prev(simVars.disc.res),
		prog_v_prev(simVars.disc.res),

		tmp(simVars.disc.res),

		pos_x(simVars.disc.res),
		pos_y(simVars.disc.res),

		posx_a(simVars.disc.res),
		posy_a(simVars.disc.res),

		posx_d(simVars.disc.res),
		posy_d(simVars.disc.res),

		benchmark_analytical_error(simVars.disc.res),

		// Init operators
		op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_u(simVars.disc.res), _parareal_data_start_v(simVars.disc.res),
		_parareal_data_fine_u(simVars.disc.res), _parareal_data_fine_v(simVars.disc.res),
		_parareal_data_coarse_u(simVars.disc.res), _parareal_data_coarse_v(simVars.disc.res),
		_parareal_data_output_u(simVars.disc.res), _parareal_data_output_v(simVars.disc.res),
		_parareal_data_error_u(simVars.disc.res), _parareal_data_error_v(simVars.disc.res)
#endif
	{
		// Calls initialization of the run (e.g. sets u, v)
		reset();

#if SWEET_PARAREAL
		parareal_setup();

#endif
	}

	virtual ~SimulationInstance()
	{
	}


	void reset()
	{
		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		//TODO: is there a reason, why this is not called in swe_rexi
		simVars.reset();

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_SPECTRAL_SPACE
		prog_u.set_spec_all(0,0);
		prog_v.set_spec_all(0,0);
#endif

		//Setup prog vars
		prog_u.set_all(0);
		prog_v.set_all(0);

		//Check if input parameters are adequate for this simulation
		if (param_use_staggering && simVars.disc.use_spectral_basis_diffs)
		{
			std::cerr << "Staggering and spectral basis not supported!" << std::endl;
			exit(1);
		}

		if (param_use_staggering && param_compute_error)
			std::cerr << "Warning: Staggered data will be interpolated to/from A-grid for exact linear solution" << std::endl;

		if (param_use_staggering)
		{
			/*
			 *              ^
			 *              |
			 *       ______v0,1_____
			 *       |             |
			 *       |			   |
			 *       |   (0.5,0.5) |
			 *  u0,0 |->  H/P0,0   |u1,0 ->
			 *(0,0.5)|	           |
			 *       |      ^      |
			 *   q0,0|______|______|
			 * (0,0)      v0,0
			 *           (0.5,0)
			 *
			 *
			 * These staggering should be used when interpolating from a staggered variable
			 * If interpolating from A grid to C staggered, use negative of displacements.
			 *
			 */
			stag_displacement[0] = 0.0;  // u_dx
			stag_displacement[1] = -0.5; // u_dy
			stag_displacement[2] = -0.5; // v_dx
			stag_displacement[3] = 0.0;  // v_dy
			stag_h[0] = -0.5;
			stag_h[1] = -0.5;
			stag_u[0] = -0.0;
			stag_u[1] = -0.5;
			stag_v[0] = -0.5;
			stag_v[1] = -0.0;
		}

		//Setup sampler for future interpolations
		sampler2D.setup(simVars.sim.domain_size, simVars.disc.res);

		//Setup semi-lag
		semiLagrangian.setup(simVars.sim.domain_size, simVars.disc.res);

		//Setup general (x,y) grid with position points
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{
		    	/* Equivalent to q position on C-grid */
				pos_x.set(j, i, ((double)i)*simVars.sim.domain_size[0]/simVars.disc.res[0]); //*simVars.sim.domain_size[0];
				pos_y.set(j, i, ((double)j)*simVars.sim.domain_size[1]/simVars.disc.res[1]); //*simVars.sim.domain_size[1];
				//std::cout << i << " " << j << " " << pos_x.get(j,i) << std::endl;
			}
		}

		//Initialize arrival points with h position
		posx_a = pos_x+0.5*simVars.disc.cell_size[0];
		posy_a = pos_y+0.5*simVars.disc.cell_size[1];

		// Set initial conditions
		for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
		{
			for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
			{

				if (param_use_staggering)
				{
					{
						// u space
						double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
						prog_u.set(j,i, BurgersValidationBenchmarks::return_u(simVars, x, y));
					}

					{
						// v space
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
						prog_v.set(j, i, BurgersValidationBenchmarks::return_v(simVars, x, y));
					}
				}
				else
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

					prog_u.set(j, i, BurgersValidationBenchmarks::return_u(simVars, x, y));
					prog_v.set(j, i, BurgersValidationBenchmarks::return_v(simVars, x, y));

				}

			}
		}


		//Initialize t-dt time step with initial condition
		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
		{
			prog_u.file_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		if (simVars.setup.input_data_filenames.size() > 1)
		{
			prog_v.file_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

	// Print output info (if gui is disabled, this is done in main
		if (simVars.misc.gui_enabled)
			timestep_output();
	}


	// Calculate the model diagnostics
	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;


		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res[0]*(double)simVars.disc.res[1]);

		// energy
		simVars.diag.total_energy =
			0.5*((
					prog_u*prog_u +
					prog_v*prog_v
				).reduce_sum_quad()) * normalization;
		/*
		 * Not used with the Burgers equation
		 *
		 * // mass
		 * simVars.diag.total_mass = -1;
		 *
		 * // potential enstropy
		 * simVars.diag.total_potential_enstrophy = -1;
		 */
	}


	void set_source( DataArray<2> &o_u_t )
	{
		double t = simVars.timecontrol.current_simulation_time;
		double tp = 2.0*M_PIl;

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 1/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 57)
		{
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = tp * std::sin(tp*k*x)*std::cos(tp*k*t)
								  + tp*std::sin(tp*k*x)*std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
								  + simVars.sim.viscosity * (tp*tp*k* std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+4*2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - 2*nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 2/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 59)
		{
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = 2*tp * std::sin(tp*k*x)*std::cos(tp*k*t)
								  + 4*tp*std::sin(tp*k*x)*std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
								  + 2*simVars.sim.viscosity * (tp*tp*k* std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*x)*cos(2*PI*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)
		 *			+ [sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 *			* [2*PI*cos(2*PI*x)*sin(2*PI*t)+2*PI*cos(2*PI*k*x)*sin(2*PI*k*t)]
		 *          - NU*[-4*PI*PI*sin(2*PI*x)*sin(2*PI*t)
		 *          - 4*PI*PI*k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 58)
		{
			double k = simVars.sim.f0;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = tp*std::sin(tp*x)*std::cos(tp*t)+tp*std::sin(tp*k*x)*std::cos(tp*k*t)
					              + (std::sin(tp*x)*std::sin(tp*t)+1/k*std::sin(tp*k*x)*std::sin(tp*k*t))
					              * (tp*std::cos(tp*x)*std::sin(tp*t)+tp*std::cos(tp*k*x)*std::sin(tp*k*t))
								  - simVars.sim.viscosity*(-tp*tp*std::sin(tp*x)*std::sin(tp*t)
								  - tp*tp*k*std::sin(tp*k*x)*std::sin(tp*k*t));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 1
		 * matching to:
		 * u(t,x,y) = t
		 */
		if (simVars.setup.scenario == 51)
		{
			o_u_t.set_all(1.0);
		}

		/*
		 * f(t,x,y) = 2*t
		 * matching to:
		 * u(t,x,y) = t^2
		 */
		if (simVars.setup.scenario == 52)
		{
			o_u_t.set_all(2.0*t);
		}

		/*
		 * f(t,x,y) = 3*t^2
		 * matching to:
		 * u(t,x,y) = t^3
		 */
		if (simVars.setup.scenario == 53)
		{
			o_u_t.set_all(3.0*t*t);
		}

		/*
		 * f(t,x,y) = 1000*sin(2*PI*x) + 1000^2*t*sin(2*PI*x)*t*cos(2*PI*x)*2*PI - 1000*NU*(-4*PI*PI*t*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = 1000*t*sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 54)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = 1000*std::sin(tp*x)+1000*1000*t*std::sin(tp*x)*t*std::cos(tp*x)*tp
							      - 1000*simVars.sim.viscosity*(-tp*tp*t*std::sin(tp*x));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*t)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*t)
		 */
		if (simVars.setup.scenario == 55)
		{
			o_u_t.set_all(tp*std::cos(tp*t));
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*k*t)
		 * matching to:
		 * u(t,x,y) = 1/k*sin(2*PI*k*t)
		 */
		if (simVars.setup.scenario == 56)
		{
			double k=simVars.sim.f0;
			o_u_t.set_all(tp*std::cos(tp*k*t));
		}

		/*
		 * f(t,x,y) = sin(2*PI*x)*cos(2*PI*x)*2*PI - NU*(-4*PI*PI*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 60)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = std::sin(tp*x)*std::cos(tp*x)*tp
							      - simVars.sim.viscosity*(-tp*tp*std::sin(tp*x));

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/* Test for 2u-grad(u)=F
		 * f(t,x,y) = (1+4*PI*PI)sin(2*PI*x)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (simVars.setup.scenario == 61)
		{
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = (1+tp*tp)*std::sin(tp*x);

					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		/*
		 * f(t,x,y) = 0.25*[SUM (cos(**)*(-PI*k)+sin(**)*(2*PI*k)^2)*# + 0.25*SUM sin(**)*# * SUM cos(**)*2*PI*k*#]
		 * mit: ** = 2*PI*k*x-PI*k*t+PI*k
		 * und: #  = EPS/sinh(0.5*PI*k*EPS)
		 * matching to:
		 * u(t,x,y) = 0.5*SUM_(k=1)^(k_max) sin(2*PI*k*x-PI*k*t+PI*k)*EPS/sinh(0.5*PI*k*EPS)
		 */
		if (simVars.setup.scenario == 62)
		{
			int kmax = 100;
			double eps = 0.1;
#if SWEET_THREADING
#pragma omp parallel for OPENMP_PAR_SIMD
#endif
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					// u space
					double x = 0.0;
					if (param_use_staggering)
					{
						x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
					}
					double tmpvar = 0;
					for (int k = 1; k < kmax; k++)
					{
						double argument = tp*k*x - M_PIl*k*t + M_PIl*k;
						tmpvar += cos(argument)*(-M_PIl*k) + sin(argument)*(tp*k)*(tp*k)*eps/sinh(M_PIl*k*eps/2);
						for (int l = 1; l < kmax; l++)
						{
							double argument2 = tp*l*x - M_PIl*l*t + M_PIl*l;
							tmpvar += sin(argument)*eps/sinh(M_PIl*k*eps/2)*cos(argument2)*eps/sinh(M_PIl*l*eps/2)*tp*l;
						}
					}

					tmpvar *= 0.25;
					o_u_t.set(j,i, tmpvar);
				}
			}
		}

		if (simVars.setup.scenario == 63)
			o_u_t.set_all(0.0);

	}

	/**
	 * Compute derivative for time stepping and store it to
	 * u_t and v_t
	 */
	void p_run_euler_timestep_update(
			const DataArray<2> &i_tmp,	///< prognostic variables
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_tmp_t,		///< time updates
			DataArray<2> &o_u_t,		///< time updates
			DataArray<2> &o_v_t,		///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		/* 2D Burgers equation [with source term]
		 * u_t + u*u_x + v*u_y = nu*(u_xx + u_yy) [+ f(t,x,y)]
		 * v_t + u*v_x + v*v_y = nu*(v_xx + v_yy) [+ g(t,x,y)]
		 */

		//TODO: staggering vs. non staggering

		DataArray<2> f(i_u.resolution);
		set_source(f);

		/*
		 * u and v updates
		 */

		if (param_semilagrangian)
		{
			o_u_t = simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
			// Delete this line if no source is used.
			o_u_t += f;
			o_v_t = simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));
		}else{
			o_u_t = -(i_u*op.diff_c_x(i_u) + i_v*op.diff_c_y(i_u));
			o_u_t += simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
			// Delete this line if no source is used.
			o_u_t += f;
			o_v_t = -(i_u*op.diff_c_x(i_v) + i_v*op.diff_c_y(i_v));
			o_v_t += simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));
		}

		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt > 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			std::cout << "Only fixed time step size supported" << std::endl;
			assert(false);
			exit(1);
		}

	}



	void run_timestep()
	{
		double dt = 0.0;

		// Only fixed time stepping supported with the Burgers equation
		assert(simVars.sim.CFL < 0);
		simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;


		if (param_semilagrangian)
		{
			dt = -simVars.sim.CFL;
			//Padding for last time step
			if (simVars.timecontrol.current_simulation_time+dt > simVars.timecontrol.max_simulation_time)
				dt = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

			//Calculate departure points
			semiLagrangian.semi_lag_departure_points_settls(
							prog_u_prev, prog_v_prev,
							prog_u,	prog_v,
							posx_a,	posy_a,
							dt,
							posx_d,	posy_d,
							stag_displacement
					);


			DataArray<2> U = prog_u;
			DataArray<2> V = prog_v;

			// Save old velocities
			prog_u_prev = prog_u;
			prog_v_prev = prog_v;

			//Now interpolate to the the departure points
			//Departure points are set for physical space
			sampler2D.bicubic_scalar( U, posx_d, posy_d, prog_u, stag_u[0], stag_u[1]);
			sampler2D.bicubic_scalar( V, posx_d, posy_d, prog_v, stag_v[0], stag_v[1]);

			if (simVars.disc.timestepping_runge_kutta_order>=0)
			{
				// run standard Runge Kutta
				timestepping.run_rk_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						tmp, prog_u, prog_v, ///< tmp is used to make use of the swe version of run_rk_timestep
						dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_runge_kutta_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
			}
			else
			{
				// run IMEX Runge Kutta
				run_timestep_imex(
						prog_u, prog_v,
						dt,
						simVars.timecontrol.current_timestep_size,
						op,
						simVars,
						simVars.timecontrol.max_simulation_time
				);
			}

		}else{

			if (simVars.disc.timestepping_runge_kutta_order>=0)
			{
				// run standard Runge Kutta
				timestepping.run_rk_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						tmp, prog_u, prog_v, ///< tmp is used to make use of the swe version of run_rk_timestep
						dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_runge_kutta_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
			}
			else
			{
				// run IMEX Runge Kutta
				run_timestep_imex(
						prog_u, prog_v,
						dt,
						simVars.timecontrol.current_timestep_size,
						op,
						simVars,
						simVars.timecontrol.max_simulation_time
				);
			}
		}

		//dt = simVars.timecontrol.current_timestep_size;

		// provide information to parameters
		simVars.timecontrol.current_timestep_size = dt;
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}

	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "TIME\t\t\t\tTOT_ENERGY";
				if (param_compute_error){
					o_ostream << "\tMAX_ABS_U\tMAX_RMS_U";
				}

				o_ostream << std::endl;

			}

			// Print timestep data to given output stream
			o_ostream << std::setprecision(8) << std::fixed << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_energy;

			if (param_compute_error){
				compute_errors();
				o_ostream << std::setprecision(8) << "\t" << benchmark_analytical_error_maxabs_u << "\t" << benchmark_analytical_error_rms_u << "\t" << prog_u.reduce_max();
			}

			if ((simVars.misc.output_file_name_prefix.size() > 0) && !simVars.parareal.enabled)
			{
				// output each time step
				if (simVars.misc.output_each_sim_seconds < 0)
					simVars.misc.output_next_sim_seconds = 0;

				if ((simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time) ||
						(simVars.timecontrol.current_simulation_time == simVars.timecontrol.max_simulation_time))
				{
					if (simVars.misc.output_each_sim_seconds > 0)
					{
						while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
							simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
					}

					double secs = simVars.timecontrol.current_simulation_time;
					double msecs = 1000000.*(simVars.timecontrol.current_simulation_time - floor(simVars.timecontrol.current_simulation_time));
					char t_buf[256];
					sprintf(	t_buf,
								"%08d.%06d",
								(int)secs, (int)msecs
						);

					std::string ss = simVars.misc.output_file_name_prefix+"_t"+t_buf;


					// write velocity field u to file
					prog_u.file_saveData_ascii((ss+"_u.csv").c_str());
					prog_u.file_saveSpectralData_ascii((ss+"_u_spec.csv").c_str());
					// write velocity field v to file
					//prog_v.file_saveData_ascii((ss+"_v.csv").c_str());

					if (param_compute_error)
					{
						benchmark_analytical_error.file_saveData_ascii((ss+".err").c_str());
						benchmark_analytical_error.file_saveSpectralData_ascii((ss+"_spec.err").c_str());
					}

				}
			}

		}
		o_ostream << std::endl;
	}


public:
	void compute_errors()
	{
			// Compute exact solution for linear part and compare with numerical solution

			// Only possible for manufactured solutions
			if (simVars.setup.scenario<51 && simVars.setup.scenario>59)
				return;

			//Analytical solution at specific time on original grid (stag or not)
			DataArray<2> ts_u(simVars.disc.res);
			DataArray<2> ts_v(simVars.disc.res);

			// Set solution
			for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
			{
				for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
				{
					if (param_use_staggering)
					{
						{
							// u space
							double x = (((double)i)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
							ts_u.set(j,i, BurgersValidationBenchmarks::return_u(simVars, x, y));
						}

						{
							// v space
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
							ts_v.set(j, i, BurgersValidationBenchmarks::return_v(simVars, x, y));
						}
					}
					else
					{
						double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
						double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];

						ts_u.set(j, i, BurgersValidationBenchmarks::return_u(simVars, x, y));
						ts_v.set(j, i, BurgersValidationBenchmarks::return_v(simVars, x, y));

					}

				}
			}

			benchmark_analytical_error = ts_u-prog_u;

			benchmark_analytical_error_rms_u = (ts_u-prog_u).reduce_rms_quad();
			benchmark_analytical_error_rms_v = (ts_v-prog_v).reduce_rms_quad();

			benchmark_analytical_error_maxabs_u = (ts_u-prog_u).reduce_maxAbs();
			benchmark_analytical_error_maxabs_v = (ts_v-prog_v).reduce_maxAbs();
		}


	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		if (simVars.timecontrol.max_simulation_time != -1 && simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time)
			return true;

		return false;
	}



	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	struct VisStuff
	{
		const DataArray<2>* data;
		const char *description;
	};

	/**
	 * Arrays for online visualisation and their textual description
	 */
	VisStuff vis_arrays[2] =
	{
			{&prog_u,	"u"},
			{&prog_v,	"v"}
	};

	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
	}


	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
		update_diagnostics();

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[1024];
		sprintf(title_string, "Time: %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				vis_arrays[id].description,
				simVars.diag.total_mass, simVars.diag.total_energy, simVars.diag.total_potential_enstrophy);
		return title_string;
	}



	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			if (simVars.misc.vis_id >= 2)
				simVars.misc.vis_id = 0;
			break;

		case 'V':
			simVars.misc.vis_id--;
			if (simVars.misc.vis_id < 0)
				simVars.misc.vis_id = 1;
			break;
		}
	}


	bool instability_detected()
	{
		return !(prog_u.reduce_all_finite() && prog_v.reduce_all_finite());
	}


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	DataArray<2> _parareal_data_start_u, _parareal_data_start_v;
	Parareal_Data_DataArrays<2> parareal_data_start;

	DataArray<2> _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_DataArrays<2> parareal_data_fine;

	DataArray<2> _parareal_data_coarse_u, _parareal_data_coarse_v;
	Parareal_Data_DataArrays<2> parareal_data_coarse;

	DataArray<2> _parareal_data_output_u, _parareal_data_output_v;
	Parareal_Data_DataArrays<2> parareal_data_output;

	DataArray<2> _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_DataArrays<2> parareal_data_error;

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		{
			DataArray<2>* data_array[2] = {&_parareal_data_start_u, &_parareal_data_start_v};
			parareal_data_start.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_coarse_u, &_parareal_data_coarse_v};
			parareal_data_coarse.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_output_u, &_parareal_data_output_v};
			parareal_data_output.setup(data_array);
		}

		{
			DataArray<2>* data_array[2] = {&_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

		output_data_valid = false;
	}



	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;

		timeframe_start = i_timeframe_start;
		timeframe_end = i_timeframe_end;
	}



	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data(
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		reset();


		*parareal_data_start.data_arrays[0] = prog_u;
		*parareal_data_start.data_arrays[1] = prog_v;

	}

	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	void sim_set_data(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		parareal_data_start = i_pararealData;

		// cast to pararealDataArray stuff
	}

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	void sim_set_mpi_comm(
			int i_mpi_comm
	)
	{
		// NOTHING TO DO HERE
	}

	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		simVars.disc.timestepping_runge_kutta_order = param_time_scheme;

		bool was_sl = false;
		if (param_semilagrangian)
		{
			param_semilagrangian = false;
			was_sl = true;
		}

		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		if (was_sl)
			param_semilagrangian = true;

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_u;
		*parareal_data_fine.data_arrays[1] = prog_v;
	}


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_Data& get_data_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_fine()" << std::endl;

		return parareal_data_fine;
	}

#endif

	/**
	 * IMEX time stepping for the coarse timestepping method
	 *
	 * The IMEX RK schemes combine an implicit and explicit time discretization.
	 * The diffusive, stiff term gets treated implicitly and the convective, non-
	 * stiff term gets treated explicitly. This results in the following system, 
	 * which is solved by this routine:
	 * (I-\nu\Delta) u(t+\tau) = u(t) - \tau (u(t)*\nabla(u(t))
	 * for u(t+\tau).
	 */
	bool run_timestep_imex(
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double& o_dt,			///< return time step size for the computed time step
			double i_timestep_size,	///< timestep size

			Operators2D &op,
			const SimulationVariables &i_simVars,

			double i_max_simulation_time = std::numeric_limits<double>::infinity()	///< limit the maximum simulation time
			)
	{
		DataArray<2> u=io_u;
		DataArray<2> v=io_v;

		// Initialize and set timestep dependent source for manufactured solution
		DataArray<2> f(io_u.resolution);
		set_source(f);
		f.requestDataInSpectralSpace();

		// Modify timestep to final time if necessary
		double& t = o_dt;
		if (simVars.timecontrol.current_simulation_time+i_timestep_size < i_max_simulation_time)
			t = i_timestep_size;
		else
			t = i_max_simulation_time-simVars.timecontrol.current_simulation_time;

		// Setting explicit right hand side and operator of the left hand side
		DataArray<2> rhs_u = u;
		DataArray<2> rhs_v = v;

		if (param_semilagrangian)
		{
			rhs_u += t*f;
		}else{
			rhs_u += - t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + t*f;
			rhs_v += - t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
		}

		if (simVars.disc.use_spectral_basis_diffs) //spectral
		{
			DataArray<2> lhs = u;
			if (param_semilagrangian)
			{
				lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spec_addScalarAll(1.0);
			}else{
				lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spec_addScalarAll(1.0);
			}

#if 1   // solving the system directly by inverting the left hand side operator
			io_u = rhs_u.spec_div_element_wise(lhs);
			io_v = rhs_v.spec_div_element_wise(lhs);
		} else { //Jacobi
			/*
			 * TODO:
			 * set these values non manually
			 */

			bool retval=false;
			int max_iters = 26000;
			double eps = 1e-7;
			double* domain_size = simVars.sim.domain_size;
			double omega = 1.0;
			retval = burgers_HelmholtzSolver::smoother_jacobi( // Velocity u
										simVars.sim.viscosity*t,
										rhs_u,
										io_u,
										domain_size,
										eps,
										max_iters,
										omega,
										0
									);
			if (!retval)
			{
				std::cout << "Did not converge!!!" << std::endl;
				exit(-1);
			}
			retval = burgers_HelmholtzSolver::smoother_jacobi( // Velocity v
										simVars.sim.viscosity*t,
										rhs_v,
										io_v,
										domain_size,
										eps,
										max_iters,
										omega,
										0
									);
			if (!retval)
			{
				std::cout << "Did not converge!!!" << std::endl;
				exit(-1);
			}

		}


#else	// making the second step of the IMEX-RK1 scheme
		DataArray<2> u1 = rhs_u.spec_div_element_wise(lhs);
		DataArray<2> v1 = rhs_v.spec_div_element_wise(lhs);

		io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
				- t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) +f*t;
		io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
				- t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
#endif

		return true;
	}

#if SWEET_PARAREAL

	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// set Runge-Kutta scheme to the chosen one for coarse time stepping
		simVars.disc.timestepping_runge_kutta_order=param_time_scheme_coarse;
		// save the fine delta t to restore it later
		double tmpCFL = simVars.sim.CFL;
		simVars.sim.CFL=timeframe_start-timeframe_end;

		// make multiple time steps in the coarse solver possible
		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// restore fine delta t
		simVars.sim.CFL = tmpCFL;
		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_u;
		*parareal_data_coarse.data_arrays[1] = prog_v;
	}



	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_Data& get_data_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_coarse()" << std::endl;

		return parareal_data_coarse;
	}



	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_difference()" << std::endl;

		for (int k = 0; k < 2; k++)
			*parareal_data_error.data_arrays[k] = *parareal_data_fine.data_arrays[k] - *parareal_data_coarse.data_arrays[k];
	}



	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: Error indicator based on the computed error norm between the
	 * old values and new values
	 */
	double compute_output_data(
			bool i_compute_convergence_test
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_output_data()" << std::endl;

		double convergence = -1;

		if (!i_compute_convergence_test || !output_data_valid)
		{
			for (int k = 0; k < 2; k++)
				*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			output_data_valid = true;
			return convergence;
		}



		for (int k = 0; k < 2; k++)
		{
			tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_u = *parareal_data_output.data_arrays[0];
		prog_v = *parareal_data_output.data_arrays[1];

		if (param_compute_error){
			compute_errors();
			std::cout << "maxabs error compared to analytical solution: " << benchmark_analytical_error_maxabs_u << std::endl;
		}

		output_data_valid = true;
		return convergence;
	}



	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_Data& get_output_data()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_output_data()" << std::endl;

		return parareal_data_output;
	}


	void output_data_file(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		Parareal_Data_DataArrays<2>& data = (Parareal_Data_DataArrays<2>&)i_data;

		std::ostringstream ss;
		ss << simVars.misc.output_file_name_prefix << "_iter" << iteration_id << "_slice" << time_slice_id << ".csv";

		std::string filename = ss.str();

//		std::cout << "filename: " << filename << std::endl;
		data.data_arrays[0]->file_saveData_ascii(filename.c_str());
		//data.data_arrays[0]->file_saveSpectralData_ascii(filename.c_str());

		std::ostringstream ss2;
		ss2 << simVars.misc.output_file_name_prefix << "_iter" << iteration_id << "_slice" << time_slice_id << ".err";

		std::string filename2 = ss2.str();

		compute_errors();

		benchmark_analytical_error.file_saveSpectralData_ascii(filename2.c_str());

		//data.data_arrays[0]->file_saveData_vtk(filename.c_str(), filename.c_str());
	}


	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		update_diagnostics();
		// Print timestep data to console
		std::cout << std::setprecision(8) << "Total energy: " << simVars.diag.total_energy << std::endl;
	}

#endif
};




int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

	MemBlockAlloc::setup();

	// program specific input parameter names
	const char *bogus_var_names[] = {
			"timestepping-scheme",
			"compute-error",
			"staggering",
			"semi-lagrangian",
			"ts2",
			nullptr
	};

	// default values for program specific parameter
	simVars.bogus.var[0] = param_time_scheme; // time stepping scheme
	simVars.bogus.var[1] = param_compute_error; // compute-error
	simVars.bogus.var[2] = param_use_staggering; // staggering
	simVars.bogus.var[3] = param_semilagrangian; // semi-lagrangian
	simVars.bogus.var[4] = param_time_scheme_coarse; //timestepping scheme for parareal coarse

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--timestepping-scheme [int]	-n: Use IMEX time stepping of order n (Right now only -1)" << std::endl;
		std::cout << "					 n: Use explicit Runge-Kutta time stepping of order n (n in 1..4), default:-1" << std::endl;
		std::cout << "	--compute-error [0/1]		Compute the errors, default:0" << std::endl;
		std::cout << "	--staggering [0/1]		Use staggered grid, default:0" << std::endl;
		std::cout << "	--semi-lagrangian [0/1]		Use semi-lagrangian formulation, default:0" << std::endl;
		std::cout << std::endl;

#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	//Burgers parameters
	param_time_scheme = simVars.bogus.var[0];
	param_compute_error = simVars.bogus.var[1];
	param_use_staggering = simVars.bogus.var[2];
	param_semilagrangian = simVars.bogus.var[3];
	param_time_scheme_coarse = simVars.bogus.var[4];


	std::ostringstream buf;
	buf << std::setprecision(14);

	{

#if SWEET_PARAREAL
		if (simVars.parareal.enabled)
		{
			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller_Serial<SimulationInstance> parareal_Controller_Serial;

			// setup controller. This initializes several simulation instances
			parareal_Controller_Serial.setup(&simVars.parareal);

			// execute the simulation
			parareal_Controller_Serial.run();
		}
		else
#endif

#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationBurgers = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationBurgers);
			delete simulationBurgers;
		}
		else
#endif
		{
			SimulationInstance *simulationBurgers = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			simulationBurgers->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
				simulationBurgers->update_diagnostics();
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
			}

			//Start counting time
			time.reset();

			//Main time loop
			while(true)
			{
				//Output data
				if (simVars.misc.verbosity > 1)
				{
					simulationBurgers->timestep_output(buf);

					std::string output = buf.str();
					buf.str("");

					std::cout << output << std::flush;
				}

				//Stop simulation if requested
				if (simulationBurgers->should_quit())
					break;

				//Main call for timestep run
				simulationBurgers->run_timestep();

				//Instability
				if (simulationBurgers->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}
			}

			//Stop counting time
			time.stop();

			double seconds = time();

			//End of output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;


			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs(simVars.diag.total_energy-diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs(simVars.diag.total_mass-diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs(simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start) << std::endl;

				if (simVars.setup.scenario == 2 || simVars.setup.scenario == 3 || simVars.setup.scenario == 4)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationBurgers->benchmark_diff_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationBurgers->benchmark_diff_v << std::endl;
				}
			}

			delete simulationBurgers;
		}
	}

	return 0;
}
