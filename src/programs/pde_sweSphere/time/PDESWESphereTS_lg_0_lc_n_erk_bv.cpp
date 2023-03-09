/*
 * Author: Pedor Peixoto <ppeixoto@usp.br>
 * based on stuff from:
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *         
 */
#include "PDESWESphereTS_lg_0_lc_n_erk_bv.hpp"


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

void PDESWESphereTS_lg_0_lc_n_erk_bv::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{

	/* Calculate velocities and stream function */
	//sweet::SphereData_Physical ug(io_phi_pert.sphereDataConfig);
	//sweet::SphereData_Physical vg(io_phi_pert.sphereDataConfig);

	//sweet::SphereData_Physical vrtg = io_vrt.toPhys();
	//sweet::SphereData_Physical divg = io_div.toPhys(); /* this should be zero! */

	//SphereData_Spectral psi = ops->inv_laplace(io_vrt);
	//SphereData_Spectral chi = ops->inv_laplace(io_div); /*this should be zero! */

	//ops->vrtdiv_to_uv(io_vrt, io_div, ug, vg);

	//ops->uv_to_vort(ug, vg);

	// standard time stepping RK
	timestepping_rk.runTimestep(
			this,
			&PDESWESphereTS_lg_0_lc_n_erk_bv::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);

}

void PDESWESphereTS_lg_0_lc_n_erk_bv::euler_timestep_update(
		const sweet::SphereData_Spectral &i_phi, //prog
		const sweet::SphereData_Spectral &i_vrt, //prog
		const sweet::SphereData_Spectral &i_div, //prog

		sweet::SphereData_Spectral &o_phi_t, //updated with euler
		sweet::SphereData_Spectral &o_vrt_t, //updated with euler
		sweet::SphereData_Spectral &o_div_t, //updated with euler

		double i_simulation_timestamp
)
{
	//zero tendencies
	o_phi_t.spectral_set_zero();
	o_vrt_t.spectral_set_zero();
	o_div_t.spectral_set_zero();


	// Calculate velocities in physical space
	sweet::SphereData_Physical u_phys, v_phys;
	ops->vrtdiv_to_uv(i_vrt, i_div, u_phys, v_phys);

	/*
	 * Calculate absolute vorticity in physical space (vrt+f)
	 */
	sweet::SphereData_Physical abs_vrtg = i_vrt.toPhys()+fg;
	//std::cout << "Vort" << std::endl;
	//i_vrt.spectral_print(6);

	// Nonlinear product (velocity * abs_vort)
	sweet::SphereData_Physical u_nl = u_phys*abs_vrtg;
	sweet::SphereData_Physical v_nl = v_phys*abs_vrtg;

	//nonlinear vort and divergence of (velocity * abs_vort)
	sweet::SphereData_Spectral vrt, div; 
	ops->uv_to_vrtdiv(u_nl, v_nl, vrt, div);

	o_vrt_t -= div; //This is basically the tendency in the Barotropic Vorticity Eq.
	//std::cout << "Vort tendency" << std::endl;
	//o_vrt_t.spectral_print(6);
	
	//Keep div constant
	//o_div_t = o_div_t;

	//Phi stays constant
	//o_phi_t = o_phi_t; 

}



bool PDESWESphereTS_lg_0_lc_n_erk_bv::setup_auto(
		sweet::SphereOperators *io_ops
)
{
	return setup(
		io_ops,
		timestepping_order
		);
}

bool PDESWESphereTS_lg_0_lc_n_erk_bv::setup(
		sweet::SphereOperators *io_ops,
		int i_timestepping_order	///< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_timestepping_order;
	timestep_size = shackTimestepControl->current_timestep_size;

	setupFG();

	return true;
}

void PDESWESphereTS_lg_0_lc_n_erk_bv::printHelp()
{
	std::cout << "	Barotropic Vorticity Equation:" << std::endl;
	std::cout << "		+ lg_0_lc_n_erk_bv" << std::endl;
}

PDESWESphereTS_lg_0_lc_n_erk_bv::PDESWESphereTS_lg_0_lc_n_erk_bv()	:
		timestepping_order(-1)
{
}



PDESWESphereTS_lg_0_lc_n_erk_bv::~PDESWESphereTS_lg_0_lc_n_erk_bv()
{
}

