#include <iomanip>
#include <math.h>
#include <string>

#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>

#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_t_erk.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_irk_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_ln_erk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"

#include "ceval.hpp"
#include "cencap.hpp"

#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>


/**
 * Write file to data and return string of file name
 */
std::string write_file(
		SphereDataCtx  &i_ctx,
		const SphereData &i_sphereData,
		const char* i_name	///< name of output variable
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	SimulationVariables* simVars = i_ctx.get_simulation_variables();

	// create copy
	SphereData sphereData(i_sphereData);

	// Write the data into the file
	const char* filename_template = simVars->iodata.output_file_name.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			simVars->timecontrol.current_timestep_size);
	sphereData.physical_file_write(buffer);

	return buffer;
}

/**
 *  Write the spectrum to file and return string of file name
 **/
std::string write_spectrum_to_file(
		SphereDataCtx &i_ctx,
		const SphereData &i_sphereData,
		const char* i_name
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	SimulationVariables* simVars = i_ctx.get_simulation_variables();

	// create copy
	SphereData sphereData(i_sphereData);

	// Write the spectrum into the file
	const char* filename_template = simVars->iodata.output_file_name.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			simVars->timecontrol.current_timestep_size);
	sphereData.spectrum_file_write(buffer);

	return buffer;
}

/** 
 *  Write the physical invariants to file
 **/
std::string write_physical_invariants_to_file(
		SphereDataCtx &i_ctx,
		const char* i_name
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	SimulationVariables* simVars = i_ctx.get_simulation_variables();

	// get the vector of time and invariants
	const std::vector<double>& time               = i_ctx.get_time();
	const std::vector<double>& mass               = i_ctx.get_mass();
	const std::vector<double>& energy             = i_ctx.get_energy();
	const std::vector<double>& potentialEnstrophy = i_ctx.get_potential_enstrophy();

	// Write the spectrum into the file
	const char* filename_template = simVars->iodata.output_file_name.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			simVars->timecontrol.current_timestep_size);

	std::ofstream file(buffer, std::ios_base::trunc);

	file << std::setprecision(20);

	for (unsigned int i = 0; i < time.size(); ++i)
		file << time[i] << " "
		<< mass[i] << " "
		<< energy[i] << " "
		<< potentialEnstrophy[i]
							  << std::endl;

	file.close();

	return buffer;
}

/**
 * Write the residuals to file
 **/
std::string write_residuals_to_file(
		SphereDataCtx &i_ctx,
		int i_rank,
		const char* i_name
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	SimulationVariables* simVars = i_ctx.get_simulation_variables();

	// get the vector of residuals
	const std::vector<std::vector<double> >& residuals = i_ctx.get_residuals();

	// Write the spectrum into the file
	const char* filename_template = simVars->iodata.output_file_name.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			simVars->timecontrol.current_timestep_size);

	std::ofstream file(buffer, std::ios_base::trunc);

	file << std::setprecision(20);

	for (unsigned int i = 0; i < residuals[i_rank].size(); ++i)
	{
		file << i << " "
				<< residuals[i_rank][i] << " "
				<< std::endl;

	}

	file.close();

	return buffer;
}



extern "C"
{
// initialization of the variables (initial condition)
void cinitial(
		SphereDataCtx *i_ctx, 
		double i_t,
		double i_dt,
		SphereDataVars *o_Y
)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	SphereData& phi_Y  = o_Y->get_phi();
	SphereData& vort_Y = o_Y->get_vort();
	SphereData& div_Y  = o_Y->get_div();

	// set the time stepping params
	i_ctx->setup_time_steps(
			0,
			i_dt
	);

	// get the SimulationVariables object from context
	SimulationVariables* simVars(i_ctx->get_simulation_variables());

	if (simVars->benchmark.use_topography)
		write_file(*i_ctx, simVars->benchmark.h_topo,  "prog_h_topo");

	// get the configuration for this level
	SphereDataConfig* data_config_nodealiasing = i_ctx->get_sphere_data_config_nodealiasing();

	// get the operator for this level
	SphereOperators* op_nodealiasing = i_ctx->get_sphere_operators_nodealiasing();

	SWESphereBenchmarksCombined *benchmarks = i_ctx->get_swe_benchmark(o_Y->get_level());

	// // instantiate phi, vort, and div without dealiasing to get the initial condition
	SphereData phi_Y_nodealiasing(data_config_nodealiasing);
	SphereData vort_Y_nodealiasing(data_config_nodealiasing);
	SphereData div_Y_nodealiasing(data_config_nodealiasing);

	phi_Y_nodealiasing.physical_set_zero();
	vort_Y_nodealiasing.physical_set_zero();
	div_Y_nodealiasing.physical_set_zero();

	// get the initial condition in phi, vort, and div
	benchmarks->setup(*simVars,
			*op_nodealiasing);
	benchmarks->setupInitialConditions(phi_Y_nodealiasing,
			vort_Y_nodealiasing,
			div_Y_nodealiasing);

	phi_Y.load_nodealiasing(phi_Y_nodealiasing);
	vort_Y.load_nodealiasing(vort_Y_nodealiasing);
	div_Y.load_nodealiasing(div_Y_nodealiasing);


	// output the configuration
	// simVars->outputConfig();

	if (rank == 0)
	{
		write_file(*i_ctx, phi_Y,  "prog_phi_init");
		write_file(*i_ctx, vort_Y, "prog_vort_init");
		write_file(*i_ctx, div_Y,  "prog_div_init");
	}
	// write_spectrum_to_file(*i_ctx, phi_Y,  "init_spectrum_phi");
	// write_spectrum_to_file(*i_ctx, vort_Y, "init_spectrum_vort");
	// write_spectrum_to_file(*i_ctx, div_Y,  "init_spectrum_div");

	SphereData phi_Y_init(phi_Y);
	SphereData phi_Y_final(phi_Y);
	phi_Y_final.request_data_spectral();
	phi_Y_final.request_data_physical();
	phi_Y_init -= phi_Y_final;
	//std::cout << "Geopotential error during conversion (infty norm) = " << phi_Y_init.physical_reduce_max_abs() << std::endl;

	SphereData div_Y_init(div_Y);
	SphereData div_Y_final(div_Y);
	div_Y_final.request_data_spectral();
	div_Y_final.request_data_physical();
	div_Y_init -= div_Y_final;
	//std::cout << "Divergence error during conversion (infty norm) = " << div_Y_init.physical_reduce_max_abs() << std::endl;

	SphereData vort_Y_init(vort_Y);
	SphereData vort_Y_final(vort_Y);
	vort_Y_final.request_data_spectral();
	vort_Y_final.request_data_physical();
	vort_Y_init -= vort_Y_final;
	//std::cout << "Vorticity error during conversion (infty norm) = " << vort_Y_init.physical_reduce_max_abs() << std::endl;

	// // get the timestepper
	// SWE_Sphere_TS_lg_irk_lc_n_erk* timestepper = i_ctx->get_lg_irk_lc_n_erk_timestepper();
	// //SWE_Sphere_TS_ln_erk* timestepper = i_ctx->get_ln_erk_timestepper();

	// double current_simulation_time = 0;
	// int nsteps                     = 0;

	// std::cout << "current_simulation_time = " << current_simulation_time
	//  	      << " i_t = " << i_t
	//  	      << std::endl;

	// std::cout << "i_dt = " << i_dt << std::endl;

	// // compute the reference solution (i.e., obtained with the reference time stepper)
	//  while (current_simulation_time < i_t)
	//    {
	//  	if (nsteps%120 == 0)
	//  	  {
	//  	    std::cout << "current_simulation_time = "
	//  		      << current_simulation_time
	//  		      << std::endl;

	// // 	    // std::string name = "prog_phi_ref_"+std::to_string(nsteps)+"_";
	// // 	    // write_file(*i_ctx, phi_Y, name.c_str());
	// // 	    // name = "prog_vort_ref_"+std::to_string(nsteps)+"_";
	// // 	    // write_file(*i_ctx, vort_Y,name.c_str());
	// // 	    // name = "prog_div_ref_"+std::to_string(nsteps)+"_";
	// // 	    // write_file(*i_ctx, div_Y, name.c_str());

	// 	  }

	// 	// solve the implicit system using the Helmholtz solver
	// 	timestepper->run_timestep(
	// 				  phi_Y,
	// 				  vort_Y,
	// 				  div_Y,
	// 				  i_dt,
	// 				  current_simulation_time
	// 				  );

	// 	if (simVars->sim.viscosity != 0)
	// 	  {
	// 	    double scalar = simVars->sim.viscosity*i_dt;
	// 	    double r      = simVars->sim.earth_radius;

	// 	    phi_Y  = phi_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
	// 	    vort_Y = vort_Y.spectral_solve_helmholtz(1.0, -scalar, r);
	// 	    div_Y  = div_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
	// 	  }

	// 	current_simulation_time += i_dt;
	// 	nsteps += 1;
	//    }

	//  std::cout << "nsteps = " << nsteps << std::endl;

	//  write_spectrum_to_file(*i_ctx, phi_Y,  "spectrum_phi_ref");
	//  write_spectrum_to_file(*i_ctx, vort_Y, "spectrum_vort_ref");
	//  write_spectrum_to_file(*i_ctx, div_Y,  "spectrum_div_ref");

	//  write_file(*i_ctx, phi_Y,  "prog_phi_ref");
	//  write_file(*i_ctx, vort_Y, "prog_vort_ref");
	//  write_file(*i_ctx, div_Y,  "prog_div_ref");

	//  FatalError("stop the simulation");

}

// finalizes the time step when libpfasst is done
// currently does nothing else than outputting the solution
void cfinal(
		SphereDataCtx *i_ctx,
		SphereDataVars *i_Y,
		int i_nnodes,
		int i_niters
)
{
	int rank   = 0;
	int nprocs = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	const SphereData& phi_Y  = i_Y->get_phi();
	const SphereData& vort_Y = i_Y->get_vort();
	const SphereData& div_Y  = i_Y->get_div();

	const int& level_id = i_Y->get_level();

	// get the SimulationVariables object from context
	SimulationVariables* simVars(i_ctx->get_simulation_variables());

	SphereData phi_Y_init(phi_Y);
	SphereData phi_Y_final(phi_Y);
	phi_Y_final.request_data_spectral();
	phi_Y_final.request_data_physical();
	phi_Y_init -= phi_Y_final;
	//std::cout << "Geopotential error during conversion (infty norm) = " << phi_Y_init.physical_reduce_max_abs() << std::endl;

	SphereData div_Y_init(div_Y);
	SphereData div_Y_final(div_Y);                                                                      div_Y_final.request_data_spectral();
	div_Y_final.request_data_physical();
	div_Y_init -= div_Y_final;
	//std::cout << "Divergence error during conversion (infty norm) = " << div_Y_init.physical_reduce_max_abs() << std::endl;

	SphereData vort_Y_init(vort_Y);
	SphereData vort_Y_final(vort_Y);
	vort_Y_final.request_data_spectral();
	vort_Y_final.request_data_physical();
	vort_Y_init -= vort_Y_final;
	//std::cout << "Vorticity error during conversion (infty norm) = " << vort_Y_init.physical_reduce_max_abs() << std::endl;

	if (nprocs == 1)
	{
		std::string filename = "prog_phi_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, phi_Y, filename.c_str());

		filename = "prog_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, vort_Y, filename.c_str());

		filename = "prog_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, div_Y, filename.c_str());

		// filename = "spectrum_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, vort_Y, filename.c_str());

		// filename = "spectrum_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, div_Y, filename.c_str());

		// filename = "spectrum_phi_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, phi_Y, filename.c_str());

		// filename = "invariants_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_physical_invariants_to_file(*i_ctx, filename.c_str());

	}
	else if (rank == 0)
	{

		std::string filename = "prog_phi_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, phi_Y, filename.c_str());

		// filename = "prog_conversion_phi_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_file(*i_ctx, phi_Y_init, filename.c_str());

		filename = "prog_vort_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, vort_Y, filename.c_str());

		filename = "prog_div_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_file(*i_ctx, div_Y, filename.c_str());

		// filename = "spectrum_vort_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, vort_Y, filename.c_str());

		// filename = "spectrum_div_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, div_Y, filename.c_str());

		// filename = "spectrum_phi_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_spectrum_to_file(*i_ctx, phi_Y, filename.c_str());

		// filename = "invariants_nprocs_"+std::to_string(nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		// write_physical_invariants_to_file(*i_ctx, filename.c_str());

	}

	if (level_id == simVars->libpfasst.nlevels-1)
	{
		std::string filename = "";
		if (nprocs == 1)
			filename = "residuals_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		else
			filename = "residuals_nprocs_"+std::to_string(nprocs)+"_current_proc_"+std::to_string(rank)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
		write_residuals_to_file(*i_ctx, rank, filename.c_str());
	}
}

// evaluates the explicit (nonlinear) piece
void ceval_f1(SphereDataVars *i_Y,
		double i_t,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F1
)
{
	const SphereData& phi_Y  = i_Y->get_phi();
	const SphereData& vort_Y = i_Y->get_vort();
	const SphereData& div_Y  = i_Y->get_div();

	SphereData& phi_F1  = o_F1->get_phi();
	SphereData& vort_F1 = o_F1->get_vort();
	SphereData& div_F1  = o_F1->get_div();

	// get the time step parameters
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

#if 0
	// return immediately if no nonlinear terms
	if (simVars->pde.use_only_linear_divergence == 1)
	{
		c_sweet_data_setval(o_F1, 0.0);
		return;
	}
#endif


	if (simVars->libpfasst.implicit_coriolis_force)
	{

		SWE_Sphere_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());

		// compute the explicit nonlinear right-hand side
		timestepper->euler_timestep_update_n(
				phi_Y,
				vort_Y,
				div_Y,
				phi_F1,
				vort_F1,
				div_F1,
				simVars->timecontrol.max_simulation_time
		);
	}
	else
	{

		if (!simVars->benchmark.use_topography)
		{

			SWE_Sphere_TS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(i_Y->get_level());

			// compute the explicit nonlinear right-hand side
			timestepper->euler_timestep_update_lc_n(
					phi_Y,
					vort_Y,
					div_Y,
					phi_F1,
					vort_F1,
					div_F1,
					simVars->timecontrol.max_simulation_time
			);
		}
		else
		{

			SWE_Sphere_TS_lg_erk_lc_n_t_erk* timestepper = i_ctx->get_lg_erk_lc_n_t_erk_timestepper(i_Y->get_level());

			// compute the explicit nonlinear right-hand side
			timestepper->euler_timestep_update_coriolis_and_nonlinear(
					phi_Y,
					vort_Y,
					div_Y,
					phi_F1,
					vort_F1,
					div_F1,
					simVars->timecontrol.max_simulation_time
			);


		}

	}
}

// evaluates the first implicit piece o_F2 = F2(i_Y)
void ceval_f2 (
		SphereDataVars *i_Y,
		double i_t,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F2
)
{
	const SphereData& phi_Y  = i_Y->get_phi();
	const SphereData& vort_Y = i_Y->get_vort();
	const SphereData& div_Y  = i_Y->get_div();

	SphereData& phi_F2  = o_F2->get_phi();
	SphereData& vort_F2 = o_F2->get_vort();
	SphereData& div_F2  = o_F2->get_div();

	// initialize the right-hand side
	c_sweet_data_setval(o_F2, 0.0);

	// get the simulation variables
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

	if (simVars->libpfasst.implicit_coriolis_force)
	{

		// get the explicit timestepper
		SWE_Sphere_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());

		// compute the linear right-hand side
		timestepper->euler_timestep_update_linear(
				phi_Y,
				vort_Y,
				div_Y,
				phi_F2,
				vort_F2,
				div_F2,
				simVars->timecontrol.max_simulation_time
		);
	}
	else
	{

		if (!simVars->benchmark.use_topography)
		{

			// get the explicit timestepper
			SWE_Sphere_TS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(i_Y->get_level());

			// compute the linear right-hand side
			timestepper->euler_timestep_update_linear(
					phi_Y,
					vort_Y,
					div_Y,
					phi_F2,
					vort_F2,
					div_F2,
					simVars->timecontrol.max_simulation_time
			);

		}
		else
		{

			// get the explicit timestepper
			SWE_Sphere_TS_lg_erk_lc_n_t_erk* timestepper = i_ctx->get_lg_erk_lc_n_t_erk_timestepper(i_Y->get_level());

			// compute the linear right-hand side
			timestepper->euler_timestep_update_linear(
					phi_Y,
					vort_Y,
					div_Y,
					phi_F2,
					vort_F2,
					div_F2,
					simVars->timecontrol.max_simulation_time
			);


		}
	}
}

// solves the first implicit system for io_Y
// then updates o_F2 with the new value of F2(io_Y)
void ccomp_f2 (
		SphereDataVars *io_Y,
		double i_t,
		double i_dt,
		SphereDataVars *i_Rhs,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F2
)
{
	SphereData& phi_Y  = io_Y->get_phi();
	SphereData& vort_Y = io_Y->get_vort();
	SphereData& div_Y  = io_Y->get_div();

	const SphereData& phi_Rhs  = i_Rhs->get_phi();
	const SphereData& vort_Rhs = i_Rhs->get_vort();
	const SphereData& div_Rhs  = i_Rhs->get_div();

	// get the simulation variables
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

	// first copy the rhs into the solution vector
	// this is needed to call the SWEET function run_ts_timestep
	phi_Y  = phi_Rhs;
	vort_Y = vort_Rhs;
	div_Y  = div_Rhs;


	if (simVars->libpfasst.use_rexi)
	{
		// get the rexi time stepper
		SWE_Sphere_TS_l_rexi* timestepper = i_ctx->get_l_rexi_timestepper(io_Y->get_level());

		// solve the implicit system using the Helmholtz solver
		timestepper->run_timestep(
				phi_Y,
				vort_Y,
				div_Y,
				i_dt,
				simVars->timecontrol.max_simulation_time
		);

	}
	else 
	{
		if (simVars->libpfasst.implicit_coriolis_force)
		{

			// get the irk timestepper
			SWE_Sphere_TS_l_irk* timestepper = i_ctx->get_l_irk_timestepper(io_Y->get_level());

			// solve the implicit system using the Helmholtz solver
			timestepper->run_timestep(
					phi_Y,
					vort_Y,
					div_Y,
					i_dt,
					simVars->timecontrol.max_simulation_time
			);

		}
		else
		{

			// get the irk timestepper
			SWE_Sphere_TS_lg_irk* timestepper = i_ctx->get_lg_irk_timestepper(io_Y->get_level());

			// solve the implicit system using the Helmholtz solver
			timestepper->run_timestep(
					phi_Y,
					vort_Y,
					div_Y,
					i_dt,
					simVars->timecontrol.max_simulation_time
			);


		}
	}


	SphereData& phi_F2  = o_F2->get_phi();
	SphereData& vort_F2 = o_F2->get_vort();
	SphereData& div_F2  = o_F2->get_div();

	phi_F2  = (phi_Y  - phi_Rhs)  / i_dt;
	vort_F2 = (vort_Y - vort_Rhs) / i_dt;
	div_F2  = (div_Y  - div_Rhs)  / i_dt;

}

// evaluates the second implicit piece o_F3 = F3(i_Y)
// currently contains the implicit artificial viscosity term
void ceval_f3 (
		SphereDataVars *i_Y,
		double i_t,
		int i_level,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F3
)
{
	const SphereData& phi_Y  = i_Y->get_phi();
	const SphereData& vort_Y = i_Y->get_vort();
	const SphereData& div_Y  = i_Y->get_div();

	SphereData& phi_F3  = o_F3->get_phi();
	SphereData& vort_F3 = o_F3->get_vort();
	SphereData& div_F3  = o_F3->get_div();

	// initialize F3 to zero in case no artificial viscosity
	c_sweet_data_setval(o_F3, 0.0);

	// get the simulation variables
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

	// no need to do anything if no artificial viscosity
	if (simVars->sim.viscosity == 0)
		return;

	// get the parameters used to apply diffusion
	const double r    = simVars->sim.sphere_radius;
	const double visc = simVars->sim.viscosity;

	phi_F3 = phi_Y;
	phi_F3.spectral_update_lambda(
			[&](
					int n, int m,
					std::complex<double> &io_data
			)
			{
		io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
			}
	);

	vort_F3 = vort_Y;
	vort_F3.spectral_update_lambda(
			[&](
					int n, int m,
					std::complex<double> &io_data
			)
			{
		io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
			}
	);

	div_F3 = div_Y;
	div_F3.spectral_update_lambda(
			[&](
					int n, int m,
					std::complex<double> &io_data
			)
			{
		io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
			}
	);
}

// solves the second implicit system for io_Y
// then updates o_F3 with the new value o_F3 = F3(io_Y)
// currently solve the implicit system formed with artificial viscosity term
void ccomp_f3 (
		SphereDataVars *io_Y,
		double i_t,
		double i_dt,
		int i_level,
		SphereDataVars *i_Rhs,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F3
)
{
	SphereData& phi_Y  = io_Y->get_phi();
	SphereData& vort_Y = io_Y->get_vort();
	SphereData& div_Y  = io_Y->get_div();

	SphereData& phi_Rhs  = i_Rhs->get_phi();
	SphereData& vort_Rhs = i_Rhs->get_vort();
	SphereData& div_Rhs  = i_Rhs->get_div();

	// initialize F3 to zero in case no artificial viscosity
	c_sweet_data_setval(o_F3, 0.0);

	// get the simulation variables
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

	// no need to do anything if no artificial viscosity
	if (simVars->sim.viscosity == 0)
		return;

	// get the parameters used to apply diffusion
	const double scalar = simVars->sim.viscosity*i_dt;
	const double r      = simVars->sim.sphere_radius;

	// solve (1-dt*visc*diff_op)*rhs = y
	phi_Y  = phi_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r);
	vort_Y = vort_Rhs.spectral_solve_helmholtz(1.0, -scalar, r);
	div_Y  = div_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r);

	// now recompute F3 with the new value of Y
	SphereData& phi_F3  = o_F3->get_phi();
	SphereData& vort_F3 = o_F3->get_vort();
	SphereData& div_F3  = o_F3->get_div();

	phi_F3  = (phi_Y  - phi_Rhs)  / i_dt;
	vort_F3 = (vort_Y - vort_Rhs) / i_dt;
	div_F3  = (div_Y  - div_Rhs)  / i_dt;

}

// applies artificial diffusion to the system
void cfinalize(
		SphereDataVars *io_Y,
		double i_t,
		double i_dt,
		SphereDataCtx *i_ctx
)
{
	// get the simulation variables
	SimulationVariables* simVars = i_ctx->get_simulation_variables();

	if (simVars->sim.viscosity == 0)
		return;

	SphereData& phi_Y  = io_Y->get_phi();
	SphereData& vort_Y = io_Y->get_vort();
	SphereData& div_Y  = io_Y->get_div();

	const double scalar = simVars->sim.viscosity*i_dt;
	const double r      = simVars->sim.sphere_radius;

	phi_Y  = phi_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
	vort_Y = vort_Y.spectral_solve_helmholtz(1.0, -scalar, r);
	div_Y  = div_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
}


}
