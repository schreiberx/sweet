/*
 * main.cpp
 *
 * PFASST SWE on the plane implementation
 *
 *  Created on: 30 Nov 2016
 *      Author: martin
 */


/*
 * See http://www.parallelintime.org/PFASST/index.html
 * for PFASST documentation
 */

#include <sweet/FatalError.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <rexi/swe_plane_rexi/SWE_Plane_REXI.hpp>
#include <benchmarks_plane/PlaneBenchmarksCombined.hpp>
#include <mpi.h>

#define WITH_MPI

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/mpi_communicator.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/implicit_sweeper.hpp>
#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/poly_interp.hpp>

#include <pfasst/interfaces.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/encap/encapsulation.hpp>
#include <pfasst/encap/imex_sweeper.hpp>

#if !SWEET_PFASST_CPP
#	error "Use --pfasst=enable for this file"
#endif

/**
 * Global variables which are shared between everything
 */
SimulationVariables simVars;

bool param_rexi_use_coriolis_formulation = true;
bool param_compute_error = false;
double param_geostr_balance_freq_multiplier = 1.0;
int param_timestepping_mode = 0.0;
int param_max_levels = 0;


class LevelSingleton
{
public:
	int level;
	PlaneDataConfig dataConfig;
	PlaneOperators op;

	// Rexi stuff
	SWE_Plane_REXI swe_plane_rexi;
};

std::vector<LevelSingleton> levelSingletons;


/**
 * Extend ImplicitSweeper in PFASST++
 * TODO: Why do we use IMEX and not Implicit?
 */
template<typename time = pfasst::time_precision>
class SimulationInstance
		: public pfasst::encap::IMEXSweeper<time>
{
    double xstart = 0.0;
    double xstop = M_PI;

public:
	// geopotential and velocities
	PlaneData prog_phi;
	PlaneData prog_u;
	PlaneData prog_v;

	PlaneDataTimesteppingRK timestepping;

	LevelSingleton *levelSingleton;

	SimulationInstance()	:
		prog_phi(1),
		prog_u(1),
		prog_v(1),
		levelSingleton(nullptr)
	{

	}


    /*
     * compute exact solution
     */
    void p_exact(
    		shared_ptr<pfasst::encap::Encapsulation<time>> u_encap,
			time t
	)
    {
//    	auto& u = pfasst::encap::as_vector<double, time>(u_encap);
//      this->exact(u, t);
    }


    void p_echo_error(time t)
    {
#if 0
      auto& qend = pfasst::encap::as_vector<double, time>(this->get_end_state());
      pfasst::encap::VectorEncapsulation<double> qex(qend.size());

      this->exact(qex, t);

      double max = 0.0;
      for (size_t i = 0; i < qend.size(); i++) {
        double d = abs(qend[i] - qex[i]);
      if (d > max) { max = d; }
      }

      auto n = this->get_controller()->get_step();
      auto k = this->get_controller()->get_iteration();


      ML_CLOG(INFO, "User", "step: " << n << " iter: " << k << " err: " << std::scientific << max);
#endif
    }


    /*
     * Override the post_sweep hook and echo the error.
     */
    void post_sweep() override
    {
#if 0
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_step_size();
      p_echo_error(t+dt);
#endif
    }



	/**
	 * PFASST
	 *
	 * Evaluate the implicit part of the ODE.
	 *
	 * This is typically called to compute the implicit part of the right hand side at the first
	 * collocation node, and on all nodes after restriction or interpolation.
	 *
	 * @param[in,out] f_impl_encap Encapsulation to store the implicit function evaluation.
	 * @param[in] u_encap Encapsulation storing the solution state at which to evaluate the
	 *     implicit part of the ODE.
	 * @param[in] t Time point of the evaluation.
	 *
	 * @note This method must be implemented in derived sweepers.
	 */
    void f_impl_eval(shared_ptr<pfasst::encap::Encapsulation<time>> io_f_impl_encap,
                     shared_ptr<pfasst::encap::Encapsulation<time>> i_u_encap,
                     time i_t
	) override
    {
    }


	/**
	 * PFASST
	 *
	 * Evaluate the explicit part of the ODE.
	 *
	 * @param[in,out] f_expl_encap Encapsulation to store the explicit function evaluation.
	 * @param[in] u_encap Encapsulation that stores the solution state at which to evaluate the
	 *     explicit part of the ODE.
	 * @param[in] t Time point of the evaluation.
	 *
	 * @note This method must be implemented in derived sweepers.
	 */
    void f_expl_eval(shared_ptr<pfasst::encap::Encapsulation<time>> io_f_expl_encap,
                     shared_ptr<pfasst::encap::Encapsulation<time>> i_u_encap,
                     time i_t
	)
    {
    }



	/**
	 * Main routine for method to be used in case of finite differences
	 */
private:
	void p_run_euler_timestep_update(
			const PlaneData &i_phi,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_phi_t,	///< time updates
			PlaneData &o_u_t,	///< time updates
			PlaneData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		/**
		 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		 * TODO: rearrange equations to phi formulation
		 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		 */
		o_dt = simVars.timecontrol.current_timestep_size;

		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * phi_t = -(h0*g)*u_x - (h0*g)*v_ym
		 * u_t = -phi_x + f*v
		 * v_t = -phi_y - f*u
		 */
		PlaneOperators &op = levelSingleton->op;

		o_u_t = -op.diff_c_x(i_phi);
		o_v_t = -op.diff_c_y(i_phi);

		o_u_t += simVars.sim.f0*i_v;
		o_v_t -= simVars.sim.f0*i_u;

		// standard update
		o_phi_t = -(op.diff_c_x(i_u) + op.diff_c_y(i_v))*(simVars.sim.h0*simVars.sim.gravitation);

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			o_phi_t += op.laplace(i_phi);
			o_u_t += op.laplace(i_u);
			o_v_t += op.laplace(i_v);
		}
	}



	/**
	 * Execute a single simulation time step
	 */
	void p_run_timestep(
			PlaneData &io_prog_phi,
			PlaneData &io_prog_u,
			PlaneData &io_prog_v,
			int i_method	/* 0: explicit, 1: REXI, 3: implicit */
	)
	{
		double o_dt;

		if (i_method == 0)
		{
			// either set time step size to 0 for autodetection or to
			// a positive value to use a fixed time step size
			simVars.timecontrol.current_timestep_size = (simVars.sim.CFL < 0 ? -simVars.sim.CFL : 0);

			// standard time stepping
			timestepping.run_rk_timestep(
					this,
					&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
					io_prog_phi, io_prog_u, io_prog_v,
					o_dt,
					simVars.timecontrol.current_timestep_size,
					simVars.disc.timestepping_runge_kutta_order,
					simVars.timecontrol.current_simulation_time,
					simVars.timecontrol.max_simulation_time
				);
		}
		else if (i_method == 1) //REXI
		{
			FatalError("TODO: implement");
#if 0
			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;

			timestepping.run_timestep_rexi( prog_h, io_prog_u, io_prog_v, o_dt, op,	simVars);
#endif
		}
		else if (param_timestepping_mode == 2) //Direct solution
		{
			if (simVars.misc.use_nonlinear_equations>0)
				FatalError("Direct solution on staggered grid not supported!");

			// Analytical solution
			assert(simVars.sim.CFL < 0);
			o_dt = -simVars.sim.CFL;
			levelSingleton->swe_plane_rexi.run_timestep_direct_solution(
					prog_phi, io_prog_u, io_prog_v,
					-simVars.sim.CFL,
					levelSingleton->op,
					simVars
			);
		}
		else if (param_timestepping_mode == 3)
		{   //  Implicit time step - Backward Euler - checked - linear only
			assert(simVars.sim.CFL < 0);

			o_dt = -simVars.sim.CFL;
			levelSingleton->swe_plane_rexi.run_timestep_implicit_ts(
					prog_phi, io_prog_u, io_prog_v,
					o_dt,
					levelSingleton->op,
					simVars
			);
		}
	}


};




int main(int i_argc, char *i_argv[])
{
	MPI_Init(&i_argc, &i_argv);
	pfasst::init(i_argc, i_argv);

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

#if 0
//	  pf.set_comm(&comm);
	// TODO: implement sweeper and transfer
	pf.add_level(sweeper, transfer);
#endif

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


	// input parameter names (specific ones for this program)
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

	param_rexi_use_coriolis_formulation = simVars.bogus.var[0];
	assert (param_rexi_use_coriolis_formulation == 0 || param_rexi_use_coriolis_formulation == 1);
	param_compute_error = simVars.bogus.var[1];


	if (simVars.setup.benchmark_scenario_id == 1)
	{
		std::cout << "WARNING: OVERRIDING PARAMETERS FOR GALEWSKI BENCHMARK" << std::endl;

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

	levelSingletons.resize(3);


	levelSingletons[0].dataConfig.setupAutoPhysicalSpace(
			simVars.disc.res_spectral[0],
			simVars.disc.res_spectral[1]
	);

	levelSingletons[0].level = 0;

	levelSingletons[0].op.setup(
			&(levelSingletons[0].dataConfig),
			simVars.sim.domain_size,
			simVars.disc.use_spectral_basis_diffs
		);

	for (int i = 1; i < param_max_levels; i++)
	{
		levelSingletons[i].dataConfig.setupAdditionalModes(
				&(levelSingletons[i-1].dataConfig),
				-1,
				-1
		);
	}




	/**********************************************************
	 * SETUP top level data
	 **********************************************************/

	SimulationInstance<pfasst::time_precision> topLevelData;
	topLevelData.prog_phi.setup(&levelSingletons[0].dataConfig);
	topLevelData.prog_u.setup(&levelSingletons[0].dataConfig);
	topLevelData.prog_v.setup(&levelSingletons[0].dataConfig);

	PlaneBenchmarksCombined::setupInitialConditions(
			topLevelData.prog_phi,
			topLevelData.prog_u,
			topLevelData.prog_v,
			simVars,
			levelSingletons[0].op
		);

	topLevelData.prog_phi = topLevelData.prog_phi*simVars.sim.gravitation;


#if 0
	// TODO: setup the rest correctly
	pf.run();
#endif
}

