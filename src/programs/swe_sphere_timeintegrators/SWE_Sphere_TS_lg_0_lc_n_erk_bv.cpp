/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 *         
 */
#include "SWE_Sphere_TS_lg_0_lc_n_erk_bv.hpp"


/*
 * Barotropic vorticity equation implementation
 * Details from here: https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/barotropic.pdf
 *
 *  The main prognostic variable will be vorticity (vrt)
 *  We are basically solving only the vorticity equation of the SWE equation written in vort-div formulation
 *  The flow is non-divergent, so div=0 everywhere, therefore, also, we don't need phi (fuild depth)
 *  We only need the streamfunction $\psi$, 
 *  but not the the velocity potential (which is zero for this equation).
 * In summary:
 *  - Main prognostic: vrt
 *  - Zero variables: chi, div, phi_pert
 *  - Constant variables: phi_pert
 *  - Diagnostic variables: u, v (obtained from psi and chi)
 */

void SWE_Sphere_TS_lg_0_lc_n_erk_bv::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{

	/* Calculate velocities and stream function */
	SphereData_Physical ug(io_phi_pert.sphereDataConfig);
	SphereData_Physical vg(io_phi_pert.sphereDataConfig);

	SphereData_Physical vrtg = io_vrt.toPhys();
	SphereData_Physical divg = io_div.toPhys(); /* this should be zero! */

	SphereData_Spectral psi = op.inv_laplace(io_vrt)/simVars.sim.sphere_radius;
	SphereData_Spectral chi = op.inv_laplace(io_div)/simVars.sim.sphere_radius; /*this should be zero! */


	op.vrtdiv_to_uv(io_vrt, io_div, ug, vg);

	op.uv_to_vort(ug, vg);

	// standard time stepping RK
	timestepping_rk.run_timestep(
			this,
			&SWE_Sphere_TS_lg_0_lc_n_erk_bv::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);

}

void SWE_Sphere_TS_lg_0_lc_n_erk_bv::euler_timestep_update(
		const SphereData_Spectral &i_phi, //prog
		const SphereData_Spectral &i_vrt, //prog
		const SphereData_Spectral &i_div, //prog

		SphereData_Spectral &o_phi_t, //updated with euler
		SphereData_Spectral &o_vrt_t, //updated with euler
		SphereData_Spectral &o_div_t, //updated with euler

		double i_simulation_timestamp
)
{
	//zero tendencies
	o_phi_t.spectral_set_zero();
	o_vrt_t.spectral_set_zero();
	o_div_t.spectral_set_zero();

	//i_vrt.spectral_print(5);
	// Calculate velocities in physical space
	SphereData_Physical u_phys, v_phys;
	op.vrtdiv_to_uv(i_vrt, i_div, u_phys, v_phys);
	//v_phys.physical_print();
	/*
	 * Calculate absolute vorticity in physical space (vrt+f)
	 */
	SphereData_Physical abs_vrtg = i_vrt.toPhys()+op.fg;

	// Nonlinear product (velocity * abs_vort)
	SphereData_Physical u_nl = u_phys*abs_vrtg;
	SphereData_Physical v_nl = v_phys*abs_vrtg;

	//nonlinear vort and divergence of (velocity * abs_vort)
	SphereData_Spectral vrt, div; 
	op.uv_to_vrtdiv(u_nl, v_nl, vrt, div);
	o_vrt_t -= div; //This is basically the tendency in the Barotropic Vorticity Eq.

	//Keep div constant
	//o_div_t = o_div_t;

	//Phi stays constant
	//o_phi_t = o_phi_t; 

}

/*
 * Setup
 */
void SWE_Sphere_TS_lg_0_lc_n_erk_bv::setup(
		int i_timestepping_order	///< order of RK time stepping method
)
{
	timestepping_order = i_timestepping_order;
	timestep_size = simVars.timecontrol.current_timestep_size;
}



void SWE_Sphere_TS_lg_0_lc_n_erk_bv::setup_auto()
{
	int version = 0;
	if (simVars.disc.timestepping_method == "SWE_Sphere_TS_lg_0_lc_n_erk_bv")
		version = 1;

	setup(
		simVars.disc.timestepping_order
		);
}

void SWE_Sphere_TS_lg_0_lc_n_erk_bv::print_help()
{
	std::cout << "	Borotropic Vorticity Equation:" << std::endl;
	std::cout << "		+ lg_0_lc_n_erk_bv" << std::endl;
}

SWE_Sphere_TS_lg_0_lc_n_erk_bv::SWE_Sphere_TS_lg_0_lc_n_erk_bv(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1)
		
{
}



SWE_Sphere_TS_lg_0_lc_n_erk_bv::~SWE_Sphere_TS_lg_0_lc_n_erk_bv()
{
}

