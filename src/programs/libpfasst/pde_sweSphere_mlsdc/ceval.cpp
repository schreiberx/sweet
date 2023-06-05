#include <iomanip>
#include <math.h>
#include <string>
#include <utility>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>

#include "../../pde_sweSphere/PDESWESphere_BenchmarksCombined.hpp"
#include "../../pde_sweSphere/time/PDESWESphereTS_l_irk.hpp"
#include "../../pde_sweSphere/time/PDESWESphereTS_lg_irk.hpp"
#include "../../pde_sweSphere/time/PDESWESphereTS_ln_erk.hpp"

#include "ceval.hpp"
#include "cencap.hpp"

extern "C"
{
#include <mpi.h>
}


/**
 * Write data to file and return string of file name
 */
std::string write_file(
		SphereDataCtx  *i_ctx,
		const sweet::SphereData_Spectral &i_sphereData,
		const char* i_name	///< name of output variable
)
{
	char buffer[1024];

	// get the pointer to the Simulation Variables object
	//sweet::ShackDictionary* shackDict = i_ctx->get_simulation_variables();

	// create copy
	sweet::SphereData_Spectral sphereData(i_sphereData);

	// Write the data into the file
	const char* filename_template = i_ctx->shackIOData->output_file_name.c_str();
	sprintf(buffer,
			filename_template,
			i_name,
			i_ctx->shackTimestepControl->current_simulation_time*i_ctx->shackIOData->output_time_scale);

    sphereData.file_write_binary_spectral(buffer);

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

	sweet::SphereData_Spectral& phi_pert_Y  = o_Y->get_phi_pert();
	sweet::SphereData_Spectral& vrt_Y = o_Y->get_vrt();
	sweet::SphereData_Spectral& div_Y  = o_Y->get_div();

	PDESWESphere_BenchmarksCombined *benchmarks = i_ctx->get_swe_benchmark(o_Y->get_level());

	sweet::SphereOperators* ops = i_ctx->get_sphere_operators(i_ctx->get_number_of_levels() - 1);
	benchmarks->setup_1_registerAllBenchmark();
	benchmarks->setup_2_shackRegistration(i_ctx->shackDict);
	benchmarks->setup_3_benchmarkDetection();
	benchmarks->setup_4_benchmarkSetup_1_withoutOps();
	benchmarks->setup_5_benchmarkSetup_2_withOps(ops);

	benchmarks->benchmark->getInitialState(phi_pert_Y, vrt_Y, div_Y);


	// output the configuration
	i_ctx->shackDict->printShackData();

	if (rank == 0)
	{
		if (i_ctx->shackIOData->output_each_sim_seconds >= 0) {
			write_file(i_ctx, phi_pert_Y, "prog_phi_pert");
			write_file(i_ctx, vrt_Y, "prog_vrt");
			write_file(i_ctx, div_Y, "prog_div");
		}
		if (i_ctx->shackIOData->output_each_sim_seconds < 0) {
		    // do not write output
		    i_ctx->shackIOData->output_next_sim_seconds = i_ctx->shackTimestepControl->max_simulation_time;
		}
		else if (i_ctx->shackIOData->output_each_sim_seconds > 0) {
		    // write output every output_each_sim_seconds
		    // next output time is thus equal to output_each_sim_seconds
		    i_ctx->shackIOData->output_next_sim_seconds = i_ctx->shackIOData->output_each_sim_seconds;
		}
		else {
		    // output at every time step
		    i_ctx->shackIOData->output_next_sim_seconds = i_ctx->shackTimestepControl->current_timestepSize;
		}
	}

	sweet::SphereData_Spectral phi_pert_Y_init(phi_pert_Y);
	sweet::SphereData_Spectral phi_pert_Y_final(phi_pert_Y);
	phi_pert_Y_init -= phi_pert_Y_final;

	sweet::SphereData_Spectral div_Y_init(div_Y);
	sweet::SphereData_Spectral div_Y_final(div_Y);
	div_Y_init -= div_Y_final;

	sweet::SphereData_Spectral vrt_Y_init(vrt_Y);
	sweet::SphereData_Spectral vrt_Y_final(vrt_Y);
	vrt_Y_init -= vrt_Y_final;

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

	const sweet::SphereData_Spectral& phi_pert_Y  = i_Y->get_phi_pert();
	const sweet::SphereData_Spectral& vrt_Y = i_Y->get_vrt();
	const sweet::SphereData_Spectral& div_Y  = i_Y->get_div();

	const int& level_id = i_Y->get_level();

	sweet::SphereData_Spectral phi_pert_Y_init(phi_pert_Y);
	sweet::SphereData_Spectral phi_pert_Y_final(phi_pert_Y);
	phi_pert_Y_init -= phi_pert_Y_final;
	
	sweet::SphereData_Spectral div_Y_init(div_Y);
	sweet::SphereData_Spectral div_Y_final(div_Y);
	div_Y_init -= div_Y_final;

	sweet::SphereData_Spectral vrt_Y_init(vrt_Y);
	sweet::SphereData_Spectral vrt_Y_final(vrt_Y);
	vrt_Y_init -= vrt_Y_final;

	if (i_ctx->shackIOData->output_each_sim_seconds < 0) {
		// do not write output
		return;
	}

	if (level_id == i_ctx->shackLibPFASST->nlevels-1)
	{
		if (rank == 0)
		{
			std::string filename = "prog_phi_pert";
			write_file(i_ctx, phi_pert_Y, filename.c_str());

			filename = "prog_vrt";
			write_file(i_ctx, vrt_Y, filename.c_str());

			filename = "prog_div";
			write_file(i_ctx, div_Y, filename.c_str());

		}
	}
}

// evaluates the explicit (nonlinear) piece
void ceval_f1(SphereDataVars *i_Y,
		double i_t,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F1
)
{
	const sweet::SphereData_Spectral& phi_pert_Y  = i_Y->get_phi_pert();
	const sweet::SphereData_Spectral& vrt_Y = i_Y->get_vrt();
	const sweet::SphereData_Spectral& div_Y  = i_Y->get_div();

	int level = i_Y->get_level();

	sweet::SphereData_Spectral& phi_pert_F1  = o_F1->get_phi_pert();
	sweet::SphereData_Spectral& vrt_F1 = o_F1->get_vrt();
	sweet::SphereData_Spectral& div_F1  = o_F1->get_div();

	PDESWESphereTS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(level);
	// compute the explicit nonlinear right-hand side
	timestepper->euler_timestep_update_lc_n(
			phi_pert_Y,
			vrt_Y,
			div_Y,
			phi_pert_F1,
			vrt_F1,
			div_F1,
			i_ctx->shackTimestepControl->current_simulation_time
	);
}


// evaluates the implicit (linear) piece
void ceval_f2(SphereDataVars *i_Y,
		double i_t,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F2
)
{
	const sweet::SphereData_Spectral& phi_pert_Y  = i_Y->get_phi_pert();
	const sweet::SphereData_Spectral& vrt_Y = i_Y->get_vrt();
	const sweet::SphereData_Spectral& div_Y  = i_Y->get_div();

	int level = i_Y->get_level();

	sweet::SphereData_Spectral& phi_pert_F2  = o_F2->get_phi_pert();
	sweet::SphereData_Spectral& vrt_F2 = o_F2->get_vrt();
	sweet::SphereData_Spectral& div_F2  = o_F2->get_div();

	PDESWESphereTS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(level);
	// compute the linear right-hand side
	timestepper->euler_timestep_update_linear(
			phi_pert_Y,
			vrt_Y,
			div_Y,
			phi_pert_F2,
			vrt_F2,
			div_F2,
			i_ctx->shackTimestepControl->current_simulation_time
	);
}

// solves the first implicit system for io_Y
// then updates o_F2 with the new value of F2(io_Y)
void ccomp_f2(
		SphereDataVars *io_Y,
		double i_t,
		double i_dtq,
		SphereDataVars *i_Rhs,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F2
)
{
	sweet::SphereData_Spectral& phi_pert_Y = io_Y->get_phi_pert();
	sweet::SphereData_Spectral& vrt_Y = io_Y->get_vrt();
	sweet::SphereData_Spectral& div_Y = io_Y->get_div();

	const sweet::SphereData_Spectral& phi_pert_Rhs  = i_Rhs->get_phi_pert();
	const sweet::SphereData_Spectral& vrt_Rhs = i_Rhs->get_vrt();
	const sweet::SphereData_Spectral& div_Rhs  = i_Rhs->get_div();

	int level = io_Y->get_level();

	// first copy the rhs into the solution vector
	// this is needed to call the SWEET function run_timestep
	phi_pert_Y = phi_pert_Rhs;
	vrt_Y = vrt_Rhs;
	div_Y = div_Rhs;

	if (i_dtq == 0)
	{
		// quadrature weight is zero -> return trivial solution
		// y = rhs (already done), f = 0.0
		c_sweet_data_setval(o_F2, 0.0);
		return;
	}

	PDESWESphereTS_lg_irk* timestepper = i_ctx->get_lg_irk_timestepper(level);
	// solve the implicit system using the Helmholtz solver
	timestepper->runTimestep(
					phi_pert_Y,
					vrt_Y,
					div_Y,
					i_dtq,
					i_ctx->shackTimestepControl->max_simulation_time
					);

	sweet::SphereData_Spectral& phi_pert_F2  = o_F2->get_phi_pert();
	sweet::SphereData_Spectral& vrt_F2 = o_F2->get_vrt();
	sweet::SphereData_Spectral& div_F2  = o_F2->get_div();
	
	phi_pert_F2 = (phi_pert_Y - phi_pert_Rhs) / i_dtq;
	vrt_F2      = (vrt_Y - vrt_Rhs) / i_dtq;
	div_F2      = (div_Y - div_Rhs) / i_dtq;

	return;
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
	const sweet::SphereData_Spectral& phi_pert_Y  = i_Y->get_phi_pert();
	const sweet::SphereData_Spectral& vrt_Y = i_Y->get_vrt();
	const sweet::SphereData_Spectral& div_Y  = i_Y->get_div();

	sweet::SphereData_Spectral& phi_pert_F3  = o_F3->get_phi_pert();
	sweet::SphereData_Spectral& vrt_F3 = o_F3->get_vrt();
	sweet::SphereData_Spectral& div_F3  = o_F3->get_div();

	// initialize F3 to zero in case of no artificial viscosity
	c_sweet_data_setval(o_F3, 0.0);

	// get the parameters used to apply diffusion
	const double visc2 = i_ctx->shackLibPFASST->hyperviscosity_2[i_level];
	const double visc4 = i_ctx->shackLibPFASST->hyperviscosity_4[i_level];
	const double visc6 = i_ctx->shackLibPFASST->hyperviscosity_6[i_level];
	const double visc8 = i_ctx->shackLibPFASST->hyperviscosity_8[i_level];
	// no need to do anything if no artificial viscosity
	if (!(visc2 || visc4 || visc6 || visc8))
	{
		return;
	}
	const double r = i_ctx->shackSphereDataOps->sphere_radius;

	// lambda for applying viscosity
	auto viscosity_applier = [&](int n, int m, std::complex<double> &io_data) 
	    {
			double laplace_op_2 = -(double)n*((double)n+1.0)/(r*r);
			double laplace_op_4 = laplace_op_2 * laplace_op_2;
			double laplace_op_6 = laplace_op_4 * laplace_op_2;
			double laplace_op_8 = laplace_op_4 * laplace_op_4;

			io_data *= (
				visc2 * laplace_op_2 + 
				(-visc4) * laplace_op_4 +
			    visc6 * laplace_op_6 +
				(-visc8) * laplace_op_8);
		};
	
	// only do this for the fields that should get viscosity applied
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("phi_pert"))
	{
		phi_pert_F3 = phi_pert_Y;
		phi_pert_F3.spectral_update_lambda(viscosity_applier);
	}
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("vrt"))
	{
		vrt_F3 = vrt_Y;
		vrt_F3.spectral_update_lambda(viscosity_applier);
	}
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("div"))
	{
		div_F3 = div_Y;
		div_F3.spectral_update_lambda(viscosity_applier);
	}	
}

// solves the second implicit system for io_Y
// then updates o_F3 with the new value o_F3 = F3(io_Y)
// currently solve the implicit system formed with artificial viscosity term
void ccomp_f3 (
		SphereDataVars *io_Y,
		double i_t,
		double i_dtq,
		int i_level,
		SphereDataVars *i_Rhs,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F3
)
{
	// set y = rhs, f = 0.0
	sweet::SphereData_Spectral& phi_pert_Y  = io_Y->get_phi_pert();
	sweet::SphereData_Spectral& vrt_Y = io_Y->get_vrt();
	sweet::SphereData_Spectral& div_Y  = io_Y->get_div();

	sweet::SphereData_Spectral& phi_pert_Rhs  = i_Rhs->get_phi_pert();
	sweet::SphereData_Spectral& vrt_Rhs = i_Rhs->get_vrt();
	sweet::SphereData_Spectral& div_Rhs  = i_Rhs->get_div();

	phi_pert_Y = phi_pert_Rhs;
	vrt_Y = vrt_Rhs;
	div_Y = div_Rhs;

	// initialize F3 to zero in case of no artificial viscosity
	c_sweet_data_setval(o_F3, 0.0);

	// get the parameters used to apply diffusion
	const double visc2  = i_ctx->shackLibPFASST->hyperviscosity_2[i_level] * i_dtq;
	const double visc4  = i_ctx->shackLibPFASST->hyperviscosity_4[i_level] * i_dtq;
	const double visc6  = i_ctx->shackLibPFASST->hyperviscosity_6[i_level] * i_dtq;
	const double visc8  = i_ctx->shackLibPFASST->hyperviscosity_8[i_level] * i_dtq;
	// no need to do anything if no artificial viscosity
	if (!(visc2 || visc4 || visc6 || visc8))
	{
		return;
	}
	const double r = i_ctx->shackSphereDataOps->sphere_radius;

	const std::array<double, 4> visc_factors {visc2, -visc4, +visc6, -visc8};

	// for all values where viscosity should get applied
	// 1. solve (1-dt*visc*diff_op)*rhs = y
	// 2. recompute F3 with the new value of Y
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("phi_pert"))
	{
		phi_pert_Y  = phi_pert_Rhs.spectral_solve_helmholtz_higher_order(1.0, visc_factors, r);
		sweet::SphereData_Spectral& phi_pert_F3 = o_F3->get_phi_pert();
		phi_pert_F3 = (phi_pert_Y - phi_pert_Rhs) / i_dtq;
	}
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("vrt"))
	{
		vrt_Y = vrt_Rhs.spectral_solve_helmholtz_higher_order(1.0, visc_factors, r);
		sweet::SphereData_Spectral& vrt_F3 = o_F3->get_vrt();
		vrt_F3 = (vrt_Y - vrt_Rhs) / i_dtq;
	}
	if (i_ctx->shackLibPFASST->hyperviscosity_on_field.at("div"))
	{
		div_Y  = div_Rhs.spectral_solve_helmholtz_higher_order(1.0, visc_factors, r);
		sweet::SphereData_Spectral& div_F3 = o_F3->get_div();
		div_F3 = (div_Y - div_Rhs) / i_dtq;
	}
}

}
