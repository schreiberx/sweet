/*
 * main.cpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */


#include <sweet/FatalError.hpp>
#include <sweet/SimulationVariables.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneData.hpp>



SimulationVariables simVars;

bool param_use_coriolis_formulation = true;
bool param_compute_error = false;
double param_geostr_balance_freq_multiplier = 1.0;

int param_max_levels = 0;


class LevelSingletons
{
public:
	PlaneDataConfig planeDataConfigArray;
	PlaneOperators op;
};

std::vector<LevelInstances> levelSingletons;


class SimulationInstance
{
public:
	// MRes level (0 is the highest resolution)
	int level;

	// geopotential and velocities
	PlaneData prog_phi;
	PlaneData prog_u;
	PlaneData prog_v;

	PlaneDataTimesteppingRK timestepping;

	LevelSingleton *levelSingleton;

#if 0
	inline
	PlaneData f(PlaneData i_sphData)
	{
		return levelSingletons[level].op.mu(i_sphData*2.0*simVars.sim.coriolis_omega);
	}
#endif
	// Main routine for method to be used in case of finite differences
	void p_run_euler_timestep_update(
			const PlaneData &i_h,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_h_t,	///< time updates
			PlaneData &o_u_t,	///< time updates
			PlaneData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (!simVars.pde.use_nonlinear_equations)
		{
			// use Robert functions for velocity
			// linear equations
			o_h_t = -(op.robert_div_lon(i_u)+op.robert_div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

			o_u_t = -op.robert_grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
			o_v_t = -op.robert_grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

			if (simVars.sim.coriolis_omega != 0)
			{
				o_u_t += f(i_v);
				o_v_t -= f(i_u);
			}
		}
		else
		{
			assert(simVars.sim.earth_radius > 0);
			assert(simVars.sim.gravitation);

			if (simVars.misc.plane_use_robert_functions)
			{
				FatalError("Only non-robert formulation is supported so far for non-linear SWE on plane!");
				// TODO: rewrite for robert functions
				// TODO: Also initialize velocities correctly
			}

			/*
			 * Height
			 */
			// non-linear equations
			o_h_t = -(op.div_lon(i_h*i_u)+op.div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);

			/*
			 * Velocity
			 */
			// linear terms
			o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
			o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

			if (simVars.sim.coriolis_omega != 0)
			{
				o_u_t += f(i_v);
				o_v_t -= f(i_u);
			}

			// non-linear terms
			o_u_t -= (i_u*op.grad_lon(i_u) + i_v*op.grad_lat(i_u))*(1.0/simVars.sim.earth_radius);
			o_v_t -= (i_u*op.grad_lon(i_v) + i_v*op.grad_lat(i_v))*(1.0/simVars.sim.earth_radius);
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			double scalar = simVars.sim.viscosity/(simVars.sim.earth_radius*simVars.sim.earth_radius);

			o_h_t += op.laplace(i_h)*scalar;
			o_u_t += op.laplace(i_u)*scalar;
			o_v_t += op.laplace(i_v)*scalar;
		}
	}
};



#if 0
int main(int argc, char *argv[])
{
  auto const quad_type = pfasst::quadrature::QuadratureType::GaussLobatto;

  auto transfer = make_shared<BilinearTransfer1D<>>();

  auto const nlevels   = pfasst::config::get_value<int>("nlevels", 1);
  auto const nnodes    = pfasst::config::get_value<int>("nnodes", 3);
  auto const nspace    = pfasst::config::get_value<int>("nspace", 8193);
  auto const nsteps    = pfasst::config::get_value<int>("nsteps", 16);
  auto const niters    = pfasst::config::get_value<int>("niters", 4);
  auto const dt        = pfasst::config::get_value<double>("dt", 0.1);

  auto quad    = pfasst::quadrature::quadrature_factory(nnodes, quad_type);
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(nspace);
  auto sweeper = make_shared<ImplicitHeatSweeper<double>>();

  sweeper->set_quadrature(quad);
  sweeper->set_factory(factory);

  pf.set_comm(&comm);
  pf.add_level(sweeper, transfer);

  if (nlevels > 1) {
    auto quad    = pfasst::quadrature::quadrature_factory((nnodes-1)/2+1, quad_type);
    auto factory = make_shared<pfasst::encap::VectorFactory<double>>((nspace-1)/2+1);
    auto sweeper = make_shared<ImplicitHeatSweeper<double>>();
    sweeper->set_quadrature(quad);
    sweeper->set_factory(factory);
    pf.add_level(sweeper, transfer);
  }

  pf.set_duration(0.0, nsteps*dt, dt, niters);
  pf.setup();

  auto q0 = sweeper->get_start_state();
  sweeper->exact(q0, 0.0);

  pf.run();

  return 0;
}
#endif


class SimulationPFASST	:
		public pfasst
{

};


int main(int i_argc, char *i_argv[])
{
	MPI_Init(&argc, &argv);
	pfasst::init(argc, argv);

	pfasst::mpi::MPICommunicator comm(MPI_COMM_WORLD);
	pfasst::PFASST<> pf;

#if 0
	// TODO: implement transfer class
	auto transfer = make_shared<BilinearTransfer1D<>>();
#endif

#if 0
	auto quad    = pfasst::quadrature::quadrature_factory(nnodes, quad_type);
	auto factory = make_shared<pfasst::encap::VectorFactory<double>>(nspace);
	auto sweeper = make_shared<ImplicitHeatSweeper<double>>();
#endif

//	  pf.set_comm(&comm);
	// TODO: implement sweeper and transfer
	pf.add_level(sweeper, transfer);

#if 0
	  if (nlevels > 1) {
	    auto factory = make_shared<pfasst::encap::VectorFactory<double>>((nspace-1)/2+1);
	    auto sweeper = make_shared<ImplicitHeatSweeper<double>>();
	    sweeper->set_quadrature(quad);
	    sweeper->set_factory(factory);
	    pf.add_level(sweeper, transfer);
	  }

	  pf.set_duration(0.0, nsteps*dt, dt, niters);
	  pf.setup();

	  auto q0 = sweeper->get_start_state();
	  sweeper->exact(q0, 0.0);
#endif


	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"rexi-use-coriolis-formulation",
			"compute-error",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 1;
	simVars.bogus.var[1] = 1;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << "	--compute-error [0/1]	Output errors (if available, default: 1)" << std::endl;
		std::cout << "	--rexi-use-coriolis-formulation [0/1]	Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
		return -1;
	}

	param_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_use_coriolis_formulation == 0 || param_use_coriolis_formulation == 1);
	param_compute_error = simVars.bogus.var[1];


	if (simVars.setup.benchmark_scenario_id == 1)
	{
		/// Setup Galewski parameters
		simVars.sim.coriolis_omega = 7.292e-5;
		simVars.sim.gravitation = 9.80616;
		simVars.sim.earth_radius = 6.37122e6;
		simVars.sim.h0 = 10000.0;

		simVars.misc.output_time_scale = 1.0/(60.0*60.0);
	}


	/**********************************************************
	 * SETUP SPH for all levels
	 **********************************************************/

	param_max_levels = 3;

	planeDataConfigArray.resize(3);

	planeDataConfigArray[0].setupAutoPhysicalSpace(
			simVars.disc.res_spectral[0],
			simVars.disc.res_spectral[1],
			&simVars.disc.res_physical[0],
			&simVars.disc.res_physical[1]
	);

	for (int i = 1; i < param_max_levels; i++)
	{
		int Nx, Ny;
		planeDataConfigArray[i].setupAutoPhysicalSpace(
				planeDataConfigArray[i-1].spectral_modes_m_max,
				planeDataConfigArray[i-1].spectral_modes_n_max,
				&Nx,
				&Ny
		);
	}


	/**********************************************************
	 * SETUP top level data
	 **********************************************************/

	SimulationInstance topLevelData;
	topLevelData.level = 0;
	topLevelData.prog_phi.setup(&planeDataConfigArray[0]);
	topLevelData.prog_u.setup(&planeDataConfigArray[0]);
	topLevelData.prog_v.setup(&planeDataConfigArray[0]);

	SphereBenchmarkCombined::setupInitialConditions(topLevelData.prog_phi, topLevelData.prog_u, topLevelData.prog_v, simVars, op);
	topLevelData.prog_phi = topLevelData.prog_phi*simVars.sim.gravitation;


#if 0
	// TODO: setup the rest correctly
	pf.run();
#endif
}


