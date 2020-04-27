/*
 * SWE_Sphere_TS_l_erk_na_sl.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk_na_sl.hpp"
#include <sweet/SimulationBenchmarkTiming.hpp>



void SWE_Sphere_TS_l_erk_na_sl::run_timestep_pert(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation*simVars.sim.h0;
	io_phi_pert += gh0;
	run_timestep_nonpert(io_phi_pert, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
	io_phi_pert -= gh0;
}


void SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_linear(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vort,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_t,	///< time updates
		SphereData_Spectral &o_vort_t,	///< time updates
		SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	double avgphi = simVars.sim.gravitation*simVars.sim.h0;

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vort.getSphereDataPhysical();
	SphereData_Physical divg = i_div.getSphereDataPhysical();
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vort, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vort, i_div, ug, vg);

	SphereData_Physical phig = i_phi.getSphereDataPhysical();

	SphereData_Physical tmpg1 = ug*(/*vrtg+*/op.fg);
	SphereData_Physical tmpg2 = vg*(/*vrtg+*/op.fg);

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	SphereData_Physical tmpg = o_div_t.getSphereDataPhysical();

	/*
	tmpg1 = ug*phig;
	tmpg2 = vg*phig;
	*/
	tmpg1 = ug*avgphi;
	tmpg2 = vg*avgphi;

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);
	else
		op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

/*
	SphereDataPhysical tmpg = 0.5*(ug*ug+vg*vg);

	if (simVars.misc.sphere_use_robert_functions)
		tmpg = tmpg.robert_convertToNonRobertSquared();
*/
	tmpspec = (phig/*+tmpg*/);

	o_div_t += -op.laplace(tmpspec);
}



void SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_n_erk_stage(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vrt,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_dt,	///< time updates
		SphereData_Spectral &o_vrt_dt,	///< time updates
		SphereData_Spectral &o_div_dt,	///< time updates

		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.start();
#endif

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	double gh0 = simVars.sim.gravitation*simVars.sim.h0;

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vrt.getSphereDataPhysical();
	SphereData_Physical divg = i_div.getSphereDataPhysical();
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_vortdiv_to_uv(i_vrt, i_div, ug, vg);
	else
		op.vortdiv_to_uv(i_vrt, i_div, ug, vg);

	SphereData_Physical phig = i_phi.getSphereDataPhysical();

	SphereData_Physical tmpg1 = ug*(vrtg/*+fg*/);
	SphereData_Physical tmpg2 = vg*(vrtg/*+fg*/);

	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_dt, o_vrt_dt);
	else
		op.uv_to_vortdiv(tmpg1, tmpg2, o_div_dt, o_vrt_dt);

	o_vrt_dt *= -1.0;

	/*
	 * SL VORT TREATMENT!!!
	 *
	 * See also Andre Robert, 1980,
	 * "A stable numerical integration scheme for the primitive meteorological equations"
	 * Eq. (31)
	 */
	o_vrt_dt -= SphereData_Spectral(i_vrt.getSphereDataPhysical()*i_div.getSphereDataPhysical());

	tmpg1 = ug*(phig-gh0);
	tmpg2 = vg*(phig-gh0);

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	if (simVars.misc.sphere_use_robert_functions)
		op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_dt);
	else
		op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_dt);

	o_phi_dt *= -1.0;

	SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	if (simVars.misc.sphere_use_robert_functions)
		tmpg = tmpg.robert_convertToNonRobertSquared();

	tmpspec = /*phig+*/tmpg;

	o_div_dt += -op.laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.stop();
#endif
}


void SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_n(
		SphereData_Spectral &io_U_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_U_vrt,		///< prognostic variables
		SphereData_Spectral &io_U_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	SphereData_Spectral phi_pert_tmp = io_U_phi_pert;
	SphereData_Spectral vrt_tmp = io_U_vrt;
	SphereData_Spectral div_tmp = io_U_div;


	const SphereData_Config *sphereDataConfig = io_U_phi_pert.sphereDataConfig;
//	double dt_div_radius = i_dt/simVars.sim.sphere_radius;

	const SphereData_Spectral &U_phi_pert = io_U_phi_pert;
	const SphereData_Spectral &U_vrt = io_U_vrt;
	const SphereData_Spectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 * TODO: Maybe we should replace this with some 1st order SL time step or some subcycling
		 */
		U_phi_pert_prev = U_phi_pert;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}


	/*
	 * Compute Lagrangian trajectories
	 */

	SphereData_Physical U_u_lon_prev(sphereDataConfig);
	SphereData_Physical U_v_lat_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev, false);

	SphereData_Physical U_u_phys(sphereDataConfig);
	SphereData_Physical U_v_phys(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt, U_div, U_u_phys, U_v_phys, false);

	// Initialize departure points with arrival points for iterative corrections
	ScalarDataArray pos_lon_d(sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray pos_lat_d(sphereDataConfig->physical_array_data_number_of_elements);

	// Calculate departure points
	semiLagrangian.semi_lag_departure_points_settls(
			U_u_lon_prev, U_v_lat_prev,
			U_u_phys, U_v_phys,

			i_dt,
			i_simulation_timestamp,
			simVars.sim.sphere_radius,
			nullptr,

			pos_lon_d, pos_lat_d,		// OUTPUT

			op,

			simVars.disc.timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry,
			simVars.disc.semi_lagrangian_interpolation_limiter
	);

	/*
	 * Backup previous points
	 */
	U_phi_pert_prev = io_U_phi_pert;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;


	/*
	 * Apply SL method
	 */
	// phi
	SphereData_Physical R_phi_pert_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar_new(
			io_U_phi_pert.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			R_phi_pert_D_phys,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);
	io_U_phi_pert.loadSphereDataPhysical(R_phi_pert_D_phys);



	// vrt/div in u/v space
	SphereData_Physical U_u_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar_new(
			U_u_phys,
			pos_lon_d, pos_lat_d,
			U_u_D_phys,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical U_v_D_phys(sphereDataConfig);
	sphereSampler.bicubic_scalar_new(
			U_v_phys,
			pos_lon_d, pos_lat_d,
			U_v_D_phys,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	op.uv_to_vortdiv(U_u_D_phys, U_v_D_phys, io_U_vrt, io_U_div, false);


	timestepping_rk_linear.run_timestep(
			this,
			&SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_n_erk_stage,	///< pointer to function to compute euler time step updates
			phi_pert_tmp, vrt_tmp, div_tmp,
			i_dt,
			timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
			i_simulation_timestamp
		);

//	io_U_phi_pert = phi_pert_tmp;
//	io_U_vrt = vrt_tmp;
	io_U_div = div_tmp;

}





void SWE_Sphere_TS_l_erk_na_sl::run_timestep_nonpert(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vort,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		FatalError("TODO");
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		euler_timestep_update_n(
			io_phi, io_vort, io_div,
			i_dt,
			i_simulation_timestamp
		);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		euler_timestep_update_n(
			io_phi, io_vort, io_div,
			i_dt,
			i_simulation_timestamp
		);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_na_sl::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi, io_vort, io_div,
				i_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 0)
	{
		FatalError("Please specify the timestepping order via --timestepping-order=[int]");
	}
	else
	{
		FatalError("programs/swe_sphere/SWE_Sphere_TS_l_erk_na_sl.cpp: This order is not yet supported!");
	}
}



/*
 * Setup
 */
void SWE_Sphere_TS_l_erk_na_sl::setup(
		int i_order,	///< order of RK time stepping method for non-linear parts
		int i_order2	///< order of RK time stepping method for non-linear parts
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	// Setup semi-lag
	semiLagrangian.setup(op.sphereDataConfig, simVars);
}

void SWE_Sphere_TS_l_erk_na_sl::setup_auto()
{
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2);
}


SWE_Sphere_TS_l_erk_na_sl::SWE_Sphere_TS_l_erk_na_sl(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		sphereSampler(semiLagrangian.sampler2D),

		U_phi_pert_prev(i_op.sphereDataConfig),
		U_vrt_prev(i_op.sphereDataConfig),
		U_div_prev(i_op.sphereDataConfig)
{
}



SWE_Sphere_TS_l_erk_na_sl::~SWE_Sphere_TS_l_erk_na_sl()
{
}

