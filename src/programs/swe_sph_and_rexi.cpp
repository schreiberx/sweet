/*
 * AppTestSWE.hpp
 *
 *  Created on: 15 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
	#include <sweet/plane/PlaneDataConfig.hpp>
	#include <sweet/plane/PlaneData.hpp>
	#include <sweet/Convert_SphereData_To_PlaneData.hpp>
#endif

#include <benchmarks_sphere/SphereBenchmarksCombined.hpp>

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataPhysical.hpp>
#include <sweet/sphere/SphereDiagnostics.hpp>

// explicit time stepping
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>
#include <sweet/sphere/SphereDataPhysicalTimesteppingExplicitRK.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitLeapfrog.hpp>

// implicit time stepping for SWE
#include <sweet/sphere/app_swe/SWEImplicit_SPHRobert.hpp>

#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/sphere/SphereOperatorsComplex.hpp>
#include <sweet/sphere/SphereDataComplex.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalComplex.hpp>
#include <sweet/sphere/app_swe/SWESphBandedMatrixPhysicalReal.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>


#include <rexi/swe_sphere_rexi/SWE_Sphere_REXI.hpp>



SimulationVariables simVars;

//Diagnostic measures at initial stage
//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

// Plane data config
SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

// Plane data config
SphereDataConfig sphereDataConfigExtInstance;
SphereDataConfig *sphereDataConfigExt = &sphereDataConfigExtInstance;



#if SWEET_GUI
	PlaneDataConfig planeDataConfigInstance;
	PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;
#endif


/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */
bool param_use_coriolis_formulation = true;
bool param_compute_error = false;



class SimulationInstance
{
public:
	SphereOperators op;
	SphereOperatorsComplex opComplex;

	// Runge-Kutta stuff
	SphereDataTimesteppingExplicitRK spheredata_timestepping_explicit_rk;
	SphereDataPhysicalTimesteppingExplicitRK spheredataphysical_timestepping_explicit_rk;

	// Leapfrog
	SphereDataTimesteppingExplicitLeapfrog spheredata_timestepping_explicit_leapfrog;

	// Implicit time stepping solver
	SWEImplicit_SPHRobert timestepping_implicit_swe;

	SWE_Sphere_REXI swe_sphere_rexi;

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	SphereData prog_h;
	SphereData prog_u;
	SphereData prog_v;

	SphereData &prog_phi = prog_h;
	SphereData &prog_vort = prog_u;
	SphereData &prog_div = prog_v;

	SphereDataPhysical prog_phys_phi;
	SphereDataPhysical prog_phys_u;
	SphereDataPhysical prog_phys_v;


	SphereDataPhysical fg;


	REXI<> rexi;

#if SWEET_GUI
	PlaneData viz_plane_data;
#endif

	int render_primitive_id = 1;

	SphereDiagnostics sphereDiagnostics;

public:
	SimulationInstance()	:
		op(sphereDataConfig, simVars.sim.earth_radius),
		opComplex(sphereDataConfig, simVars.sim.earth_radius),
		timestepping_implicit_swe(op),
		prog_h(sphereDataConfig),
		prog_u(sphereDataConfig),
		prog_v(sphereDataConfig),
		prog_phys_phi(sphereDataConfig),
		prog_phys_u(sphereDataConfig),
		prog_phys_v(sphereDataConfig),
		fg(sphereDataConfig),

#if SWEET_GUI
		viz_plane_data(planeDataConfig),
#endif
		sphereDiagnostics(sphereDataConfig, simVars)
	{
#if 0
		SphereData tmp(sphereDataConfig);
		tmp.physical_set_all_value(2.0);

		double a = sphereDiagnostics.compute_zylinder_integral_quad(tmp);
		std::cout << a << std::endl;
		exit(1);
#endif
		reset();

	}



#if 1

	inline
	SphereData f(const SphereData &i_sphData)	const
	{
		if (simVars.sim.f_sphere)
			return i_sphData*simVars.sim.f0;
		else
			return op.mu(i_sphData)*(2.0*simVars.sim.coriolis_omega);
	}

#else
	inline
	SphereData f(const SphereData &i_sphData)	const
	{
		SphereData out(i_sphData);

		out.physical_update_lambda_gaussian_grid(
				[&](double lam, double sin_phi, double &io_data)
				{
					io_data *= sin_phi*2.0*simVars.sim.coriolis_omega;
				}
		);

		return out;
//		return op.mu(i_sphData*(2.0*simVars.sim.coriolis_omega));
	}
#endif


	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		if (simVars.pde.id == 1 || simVars.pde.id == 14)
		{
#if 0
			sphereDiagnostics.update_phi_vort_div_2_mass_energy_enstrophy_4_sphere(
					op,
					prog_phi,
					prog_vort,
					prog_div,
					simVars
			);
#else
			sphereDiagnostics.update_phi_vort_div_2_mass_energy_enstrophy_4_zylinder(
					op,
					prog_phi,
					prog_vort,
					prog_div,
					simVars
			);
#endif
		}
		else if (simVars.pde.id == 2)
		{
#if 0
			sphereDiagnostics.update_h_u_v_2_mass_energy_enstrophy_4_sphere(
					op,
					(prog_phi*(1.0/simVars.sim.gravitation)).getSphereDataPhysical(),
					prog_u.getSphereDataPhysical(),
					prog_v.getSphereDataPhysical(),
					simVars
			);
#else
			sphereDiagnostics.update_phi_vort_div_2_mass_energy_enstrophy_4_zylinder(
					op,
					(prog_phi*(1.0/simVars.sim.gravitation)).getSphereDataPhysical(),
					prog_u.getSphereDataPhysical(),
					prog_v.getSphereDataPhysical(),
					simVars
			);
#endif
		}
		else
		{
#if 0
			sphereDiagnostics.update_h_u_v_2_mass_energy_enstrophy_4_sphere(
					op,
					(simVars.pde.id == 2 ? prog_phi*(1.0/simVars.sim.gravitation) : prog_h),
					prog_u,
					prog_v,
					simVars
			);
#else
			sphereDiagnostics.update_phi_vort_div_2_mass_energy_enstrophy_4_zylinder(
					op,
					(simVars.pde.id == 2 ? prog_phi*(1.0/simVars.sim.gravitation) : prog_h),
					prog_u,
					prog_v,
					simVars
			);
#endif
		}
	}



	void reset()
	{
		simVars.reset();

		if (simVars.sim.f_sphere)
		{
			fg.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = simVars.sim.f0;
				}
			);
		}
		else
		{
			fg.physical_update_lambda_gaussian_grid(
				[&](double lon, double mu, double &o_data)
				{
					o_data = mu*2.0*simVars.sim.coriolis_omega;
				}
			);
		}

		// reset the RK time stepping buffers
		switch (simVars.disc.timestepping_method)
		{
		case simVars.disc.RUNGE_KUTTA_EXPLICIT:
			spheredata_timestepping_explicit_rk.resetAndSetup(sphereDataConfig, simVars.disc.timestepping_order);
			spheredataphysical_timestepping_explicit_rk.resetAndSetup(sphereDataConfig, simVars.disc.timestepping_order);
			break;

		case simVars.disc.LEAPFROG_EXPLICIT:
			spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);
			break;
		}


		switch (simVars.pde.id)
		{
		case 12:
			std::cout << "OVERRIDING BENCHMARK SCENARIO ID TO 11 to match PDE: Advection DIV(U.phi)" << std::endl;
			simVars.setup.benchmark_scenario_id = 11;
			break;

		case 13:
			std::cout << "OVERRIDING BENCHMARK SCENARIO ID TO 12 to match PDE: Advection U.GRAD(phi)" << std::endl;
//			simVars.setup.benchmark_scenario_id = 12;
			break;
		}

		// one month runtime
		if (simVars.timecontrol.max_simulation_time == -1)
		{
			if (simVars.setup.benchmark_scenario_id == 10)
			{
				simVars.timecontrol.max_simulation_time = 31*60*60*24;

				// 200 h
				simVars.timecontrol.max_simulation_time = 200*60*60;
			}
		}

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		simVars.misc.output_next_sim_seconds = 0;


		if (simVars.sim.CFL < 0)
			simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;


		if (simVars.timecontrol.current_timestep_size <= 0)
		{
			// TRY to guess optimal time step size

			// time step size
			if (sphereDataConfig->physical_num_lat < 256)
				simVars.timecontrol.current_timestep_size = 0.002*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
			else
				simVars.timecontrol.current_timestep_size = 0.001*simVars.sim.earth_radius/(double)sphereDataConfig->physical_num_lat;
		}


		if (simVars.misc.output_each_sim_seconds <= 0)
		{
			simVars.misc.output_each_sim_seconds = 60*30;	// output every 1/2 hour
		}


		if (simVars.disc.timestepping_method == simVars.disc.REXI)
		{
			// Override for REXI
			if (simVars.timecontrol.current_timestep_size <= 0)
			{
				std::cout << "Timestep size not positive" << std::endl;
				assert(false);
				exit(1);
			}

			simVars.misc.use_nonlinear_equations = 0;
		}


		SphereBenchmarksCombined::setupInitialConditions(prog_h, prog_u, prog_v, simVars, op);


		simVars.outputConfig();

		std::cout << std::endl;
		std::cout << "LOCAL PARAMETERS:" << std::endl;
		std::cout << " + param_compute_error: " << param_compute_error << std::endl;
		std::cout << " + param_use_coriolis_formulation: " << param_use_coriolis_formulation << std::endl;
		std::cout << std::endl;


		switch (simVars.pde.id)
		{
		case 0:
			std::cout << "PDE: SWE U-V formulation" << std::endl;
			break;

		case 1:
			std::cout << "PDE: SWE VORT/DIV formulation" << std::endl;
			{
				SphereData tmp_vort(prog_vort);
				SphereData tmp_div(prog_div);

				if (simVars.misc.sphere_use_robert_functions)
					op.robert_uv_to_vortdiv(prog_u.getSphereDataPhysical(), prog_v.getSphereDataPhysical(), tmp_vort, tmp_div);
				else
					op.uv_to_vortdiv(prog_u.getSphereDataPhysical(), prog_v.getSphereDataPhysical(), tmp_vort, tmp_div);

				prog_vort = tmp_vort;
				prog_div = tmp_div;
				prog_phi = prog_h*simVars.sim.gravitation;
			}
			break;

		case 2:
			std::cout << "PDE: SWE Advective formulation" << std::endl;
			prog_phys_phi = prog_h.getSphereDataPhysical()*simVars.sim.gravitation;
			prog_phys_u = prog_u.getSphereDataPhysical();
			prog_phys_v = prog_v.getSphereDataPhysical();

			prog_phi = prog_h*simVars.sim.gravitation;
			break;

		case 12:
			std::cout << "PDE: Advection DIV(U.phi)" << std::endl;
			simVars.setup.benchmark_scenario_id = 11;
			break;

		case 13:
			std::cout << "PDE: Advection U.GRAD(phi)" << std::endl;
			simVars.setup.benchmark_scenario_id = 12;
			break;
		}


		if (simVars.disc.timestepping_method == simVars.disc.REXI)
		{
			swe_sphere_rexi.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_M,
					simVars.rexi.rexi_L,

					sphereDataConfig,
					&simVars.sim,
					simVars.timecontrol.current_timestep_size,

					simVars.rexi.rexi_use_half_poles,
					simVars.misc.sphere_use_robert_functions,
					simVars.pde.id,

					simVars.rexi.rexi_use_extended_modes,
					simVars.rexi.rexi_normalization,

					param_use_coriolis_formulation,
					simVars.sim.f_sphere,
					simVars.rexi.rexi_sphere_solver_preallocation
				);
		}

		if (	simVars.disc.timestepping_method == simVars.disc.IMPLICIT_TIMESTEP ||
				simVars.disc.timestepping_method == simVars.disc.CRANK_NICOLSON
		)
		{
			if (simVars.sim.CFL >= 0)
				FatalError("CFL >= 0: Set negative CFL for constant time step size");

			if (simVars.disc.timestepping_order != 1 && simVars.disc.timestepping_order != 2)
				FatalError("Only order 1 implicit Euler (backward) and Crank Nicolson (order 2) supported");

			// setup CN damping factor in case CN is used
			timestepping_implicit_swe.setCrankNicolsonDampingFactor(
					simVars.disc.crank_nicolson_filter
			);

			timestepping_implicit_swe.setup(
					sphereDataConfig,
					sphereDataConfigExt,

					simVars.sim.earth_radius,
					simVars.sim.coriolis_omega,
					simVars.sim.gravitation*simVars.sim.h0,
					-simVars.sim.CFL,

					param_use_coriolis_formulation,

					simVars.disc.timestepping_order
				);
		}

		update_diagnostics();

		simVars.diag.backup_reference();
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const SphereData &i_sphereData,
			const char* i_name,	///< name of output variable
			bool i_phi_shifted
		)
	{
		char buffer[1024];

		// create copy
		SphereData sphereData(i_sphereData);

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;
	}


	void write_file_output()
	{
		if (simVars.misc.output_file_name_prefix.length() == 0)
			return;

		std::string output_filename;

		std::cout << "Simulation time: " << simVars.timecontrol.current_simulation_time << std::endl;

		switch(simVars.pde.id)
		{
		case 1:
			{
				SphereData h = prog_phi*(1.0/simVars.sim.gravitation);
				output_filename = write_file(h, "prog_h", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << " (min: " << SphereData(h).physical_reduce_min() << ", max: " << SphereData(h).physical_reduce_max() << ")" << std::endl;

				SphereDataPhysical u(sphereDataConfig);
				SphereDataPhysical v(sphereDataConfig);

				if (simVars.misc.sphere_use_robert_functions)
					op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
				else
					op.vortdiv_to_uv(prog_vort, prog_div, u, v);

				output_filename = write_file(u, "prog_u", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;

				output_filename = write_file(v, "prog_v", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;

				output_filename = write_file(prog_vort, "prog_eta", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;
			}
			break;

		default:
			{
				SphereData tmp_vort(prog_vort);
				SphereData tmp_div(prog_div);

				output_filename = write_file(prog_h, "prog_h", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << " (min: " << SphereData(prog_h).physical_reduce_min() << ", max: " << SphereData(prog_h).physical_reduce_max() << ")" << std::endl;

				output_filename = write_file(prog_u, "prog_u", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;

				output_filename = write_file(prog_v, "prog_v", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;


				if (simVars.misc.sphere_use_robert_functions)
					output_filename = write_file(op.robert_uv_to_vort(prog_u.getSphereDataPhysical(), prog_v.getSphereDataPhysical())/simVars.sim.earth_radius, "prog_eta", simVars.setup.benchmark_scenario_id == 0);
				else
					output_filename = write_file(op.vort(SphereData(prog_u), SphereData(prog_v))/simVars.sim.earth_radius, "prog_eta", simVars.setup.benchmark_scenario_id == 0);
				std::cout << output_filename << std::endl;
			}
		}



		/*
		 * write reference solution
		 */
#if 0
		if (	param_compute_error &&
				simVars.misc.use_nonlinear_equations == 0 &&
				simVars.setup.benchmark_scenario_id == 10
		)
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);


			output_filename = write_file(test_h, "reference_h", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file(test_u, "reference_u", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file(test_v, "reference_v", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;



			output_filename = write_file((SphereData(prog_h)-SphereData(test_h)), "ref_diff_h", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file((SphereData(prog_u)-SphereData(test_u)), "ref_diff_u", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;

			output_filename = write_file((SphereData(prog_v)-SphereData(test_v)), "ref_diff_v", simVars.setup.benchmark_scenario_id == 0);
			std::cout << output_filename << std::endl;
		}
#endif
	}


	void timestep_do_output()
	{
		if (param_compute_error)
		{
			if (
					simVars.setup.benchmark_scenario_id != 10 &&
					simVars.setup.benchmark_scenario_id != 101
			)
			{
				FatalError("Analytical solution not available for this benchmark");
			}

			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);

			double error_h = -1;
			double error_u = -1;
			double error_v = -1;

			if (simVars.misc.sphere_use_robert_functions && simVars.pde.id == 1)
			{
				static int blarg = 0;
				if (blarg == 0)
				{
					std::cerr << "WARNING!!!" << std::endl;
					std::cerr << "WARNING!!!" << std::endl;
					std::cerr << "WARNING: USING NON-ROBERT FORMULATION for PDE.id == 1 LEADS TO STRONG ERRORS AT THE POLES!!!" << std::endl;
					std::cerr << "WARNING!!!" << std::endl;
					std::cerr << "WARNING!!!" << std::endl;
					blarg = 1;
				}
			}

			switch(simVars.pde.id)
			{
			case 1:
			case 14:
				{
					SphereData h = prog_phi*(1.0/simVars.sim.gravitation);
					SphereDataPhysical u(sphereDataConfig);
					SphereDataPhysical v(sphereDataConfig);

					if (simVars.misc.sphere_use_robert_functions)
						op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
					else
						op.vortdiv_to_uv(prog_vort, prog_div, u, v);

					SphereDataPhysical test_ug = test_u.getSphereDataPhysical();
					SphereDataPhysical test_vg = test_v.getSphereDataPhysical();

					error_h = h.physical_reduce_max_abs(test_h);
					error_u = u.physical_reduce_max_abs(test_ug);
					error_v = v.physical_reduce_max_abs(test_vg);
				}
				break;


			case 2:
				{
					SphereData h = prog_phi*(1.0/simVars.sim.gravitation);

					error_h = h.physical_reduce_max_abs(test_h);
					error_u = prog_u.physical_reduce_max_abs(test_u);
					error_v = prog_v.physical_reduce_max_abs(test_v);
				}
				break;

			default:
				{
					error_h = (prog_h-test_h).physical_reduce_max_abs();
					error_u = (prog_u-test_u).physical_reduce_max_abs();
					error_v = (prog_v-test_v).physical_reduce_max_abs();
				}
			}

			std::cerr << "error (h,u,v):\t" << error_h << "\t" << error_u << "\t" << error_v << std::endl;
		}
		write_file_output();

		// output line break
		std::cout << std::endl;

		if (simVars.misc.verbosity > 1)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				std::cerr << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
				std::cerr << std::endl;
			}

			//Print simulation time, energy and pot enstrophy
			std::cerr << simVars.timecontrol.current_simulation_time << "\t";
			std::cerr << simVars.diag.total_mass << "\t";
			std::cerr << simVars.diag.potential_energy << "\t";
			std::cerr << simVars.diag.kinetic_energy << "\t";
			std::cerr << simVars.diag.total_energy << "\t";
			std::cerr << simVars.diag.total_potential_enstrophy << "\t";

			std::cerr << (simVars.diag.total_mass-simVars.diag.ref_total_mass)/simVars.diag.total_mass << "\t";
			std::cerr << (simVars.diag.potential_energy-simVars.diag.ref_potential_energy)/simVars.diag.potential_energy << "\t";
			std::cerr << (simVars.diag.kinetic_energy-simVars.diag.ref_kinetic_energy)/simVars.diag.kinetic_energy << "\t";
			std::cerr << (simVars.diag.total_energy-simVars.diag.total_energy)/simVars.diag.total_energy << "\t";
			std::cerr << (simVars.diag.total_potential_enstrophy-simVars.diag.total_potential_enstrophy)/simVars.diag.total_potential_enstrophy << std::endl;

			static double start_tot_energy = -1;
			if (start_tot_energy == -1)
				start_tot_energy = simVars.diag.total_energy;
		}

		if (simVars.misc.verbosity > 0)
		{
			std::cout << "prog_h min/max:\t" << SphereData(prog_h).physical_reduce_min() << ", " << SphereData(prog_h).physical_reduce_max() << std::endl;
		}


		if (	param_compute_error &&
				simVars.misc.use_nonlinear_equations == 0 &&
				simVars.setup.benchmark_scenario_id == 10
		)
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);

			std::cout << "ERRORS - time, RMS(h,u,v), MAXABS(h,u,v):\t";
			std::cout << simVars.timecontrol.current_simulation_time << "\t";

			std::cout << (SphereData(test_h)-SphereData(prog_h)).physical_reduce_rms() << "\t";
			std::cout << (SphereData(test_u)-SphereData(prog_u)).physical_reduce_rms() << "\t";
			std::cout << (SphereData(test_v)-SphereData(prog_v)).physical_reduce_rms() << "\t";

			std::cout << (SphereData(test_h)-SphereData(prog_h)).physical_reduce_max_abs() << "\t";
			std::cout << (SphereData(test_u)-SphereData(prog_u)).physical_reduce_max_abs() << "\t";
			std::cout << (SphereData(test_v)-SphereData(prog_v)).physical_reduce_max_abs() << std::endl;
		}


		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
	}



public:
	bool timestep_check_output()
	{
		std::cout << "." << std::flush;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		timestep_do_output();

		return true;
	}



public:
	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		double diff = std::abs(simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time);

		if (	simVars.timecontrol.max_simulation_time != -1 &&
				(
						simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time
						||
						diff/simVars.timecontrol.max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
				)
			)
			return true;

		return false;
	}



	bool detect_instability()
	{
		double max_abs_value = std::abs(simVars.sim.h0)*2.0;

		switch (simVars.pde.id)
		{
		case 1:
		case 2:
			max_abs_value *= simVars.sim.gravitation;
		}

		if (
				SphereData(prog_h).physical_reduce_max_abs() > max_abs_value &&
				simVars.setup.benchmark_scenario_id != 4
		)
		{
			std::cerr << "Instability detected (max abs value of h > " << max_abs_value << ")" << std::endl;
			return true;
		}

		if (SphereData(prog_h).physical_isAnyNaNorInf())
		{
			std::cerr << "Inf value detected" << std::endl;
			return true;
		}

		return false;
	}



	void outInfo(const std::string &i_string, const SphereData &i_data)
	{
		std::cout << i_string << ": " << i_data.physical_reduce_min() << ", " << i_data.physical_reduce_max() << std::endl;
	}


/*

	SphereData add_f(const SphereData &i_data)
	const
	{
		SphereData d(i_data);

		double two_omega = 2.0*simVars.sim.coriolis_omega;
		d.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data += mu*two_omega;
			}
		);

		return d;
	}
*/



	void normal_mode_analysis()
	{
		// dummy time step to get time step size
		reset();
		run_timestep();

		/*
		 * Do a normal mode analysis, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		char buffer_real[1024];
		const char* filename = simVars.misc.output_file_name_prefix.c_str();

		sprintf(buffer_real, filename, "normal_modes_physical", simVars.timecontrol.current_timestep_size*simVars.misc.output_time_scale);
		std::ofstream file(buffer_real, std::ios_base::trunc);
		std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

		std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

		// use very high precision
		file << std::setprecision(20);

		SphereData* prog[3] = {&prog_h, &prog_u, &prog_v};

		int max_prog_id = -1;

		switch(simVars.pde.id)
		{
		case 0:
		case 1:
		case 2:
			max_prog_id = 3;
			prog_h.physical_set_zero();
			prog_u.physical_set_zero();
			prog_v.physical_set_zero();
			break;

		default:
			max_prog_id = 1;
			prog_h.physical_set_zero();
			break;
		}

		int leapfrog_start_num_timesteps = 16;
		double leapfrog_start_timesteps_size;
		double leapfrog_end_timestep_size;
		double leapfrog_original_cfl;
		if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
		{
			std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
			std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
			std::cout << "Therefore, we also halve the time step size here" << std::endl;

			leapfrog_original_cfl = simVars.sim.CFL;
			// do 2 leapfrog time steps -> half time step size
			simVars.sim.CFL *= 0.5;
			simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

			leapfrog_start_timesteps_size = simVars.timecontrol.current_timestep_size/(double)leapfrog_start_num_timesteps;
			leapfrog_end_timestep_size = simVars.timecontrol.current_timestep_size;
		}

		int num_timesteps = 1;
		if (simVars.disc.normal_mode_analysis_generation >= 10)
		{
			if (simVars.timecontrol.max_timesteps_nr > 0)
				num_timesteps = simVars.timecontrol.max_timesteps_nr;
		}

		file << "# " << (num_timesteps*(-simVars.sim.CFL)) << std::endl;

		// iterate over all prognostic variables
		for (int outer_prog_id = 0; outer_prog_id < max_prog_id; outer_prog_id++)
		{
			std::cout << "normal mode analysis for prog " << outer_prog_id << std::endl;

			if (simVars.disc.normal_mode_analysis_generation == 1 || simVars.disc.normal_mode_analysis_generation == 11)
			{
				// iterate over physical space
				for (int outer_i = 0; outer_i < sphereDataConfig->physical_array_data_number_of_elements; outer_i++)
				{
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->physical_set_zero();

					// activate mode
					prog[outer_prog_id]->request_data_physical();
					prog[outer_prog_id]->physical_space_data[outer_i] = 1;

					/*
					 * RUN timestep
					 */

					if (simVars.pde.id == 2)
					{
						prog_phys_phi = prog_phi.getSphereDataPhysical();
						prog_phys_u = prog_u.getSphereDataPhysical();
						prog_phys_v = prog_v.getSphereDataPhysical();
					}

					// In case of a multi-step scheme, reset it!
					if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
					{
						spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);

						simVars.sim.CFL = -leapfrog_start_timesteps_size;
						simVars.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

						for (int i = 0; i < leapfrog_start_num_timesteps; i++)
							run_timestep();

						simVars.sim.CFL = -leapfrog_end_timestep_size;
						simVars.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

						run_timestep();

						simVars.sim.CFL = leapfrog_original_cfl;
						simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;
					}
					else
					{
						run_timestep();
					}

					for (int i = 1; i < num_timesteps; i++)
					{
//						std::cout << "Timestep " << i << std::endl;
						run_timestep();
					}


					if (simVars.disc.normal_mode_analysis_generation == 1)
					{
						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_physical();
						for (int k = 0; k < sphereDataConfig->physical_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->physical_space_data[k];
							if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->physical_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
			else if (simVars.disc.normal_mode_analysis_generation == 2 || simVars.disc.normal_mode_analysis_generation == 12)
			{

				// iterate over physical space
				for (int outer_i = 0; outer_i < sphereDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					for (int imag_i = 0; imag_i < 2; imag_i++)
					{
						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->spectral_set_zero();

						// activate mode
						prog[outer_prog_id]->request_data_spectral();
						if (imag_i)
							prog[outer_prog_id]->spectral_space_data[outer_i].imag(1);
						else
							prog[outer_prog_id]->spectral_space_data[outer_i].real(1);


						if (simVars.pde.id == 2)
						{
							prog_phys_phi = prog_phi.getSphereDataPhysical();
							prog_phys_u = prog_u.getSphereDataPhysical();
							prog_phys_v = prog_v.getSphereDataPhysical();
						}


						// In case of a multi-step scheme, reset it!
						if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
						{
							spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);

							simVars.sim.CFL = -leapfrog_start_timesteps_size;
							simVars.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

							for (int i = 0; i < leapfrog_start_num_timesteps; i++)
								run_timestep();

							simVars.sim.CFL = -leapfrog_end_timestep_size;
							simVars.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

							run_timestep();

							simVars.sim.CFL = leapfrog_original_cfl;
							simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;
						}
						else
						{
							run_timestep();
						}

						for (int i = 1; i < num_timesteps; i++)
						{
//							std::cout << "Timestep " << i << std::endl;
							run_timestep();
						}




						if (simVars.disc.normal_mode_analysis_generation == 2)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog[outer_prog_id]->request_data_spectral();
							prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
								prog[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							prog[inner_prog_id]->request_data_spectral();
							for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}
						}

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						{
							prog[inner_prog_id]->request_data_spectral();

							for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->spectral_space_data[k].imag();
								if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->spectral_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
			}
			else if (simVars.disc.normal_mode_analysis_generation == 3 || simVars.disc.normal_mode_analysis_generation == 13)
			{
				// iterate over physical space
				for (int outer_i = 0; outer_i < sphereDataConfig->spectral_array_data_number_of_elements; outer_i++)
				{
					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
						prog[inner_prog_id]->spectral_set_zero();

					// activate mode via real coefficient
					prog[outer_prog_id]->request_data_spectral();
					prog[outer_prog_id]->spectral_space_data[outer_i].real(1);


					if (simVars.pde.id == 2)
					{
						prog_phys_phi = prog_phi.getSphereDataPhysical();
						prog_phys_u = prog_u.getSphereDataPhysical();
						prog_phys_v = prog_v.getSphereDataPhysical();

						/*
						 * This is the vort/div formulation which uses the velocities as
						 * prognostic values
						 *
						 * We assume that the delta function was applied to the vorticity/divergence
						 * and convert this to the velocity
						 */
						if (simVars.misc.sphere_use_robert_functions)
							op.robert_vortdiv_to_uv(prog_phys_u, prog_phys_v, prog_phys_u, prog_phys_v);
						else
							op.vortdiv_to_uv(prog_phys_u, prog_phys_v, prog_phys_u, prog_phys_v);
					}

					// In case of a multi-step scheme, reset it!
					if (simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
					{
						spheredata_timestepping_explicit_leapfrog.resetAndSetup(prog_h, simVars.disc.timestepping_order, simVars.disc.leapfrog_robert_asselin_filter);

						simVars.sim.CFL = -leapfrog_start_timesteps_size;
						simVars.timecontrol.current_timestep_size = leapfrog_start_timesteps_size;

						for (int i = 0; i < leapfrog_start_num_timesteps; i++)
							run_timestep();

						simVars.sim.CFL = -leapfrog_end_timestep_size;
						simVars.timecontrol.current_timestep_size = leapfrog_end_timestep_size;

						run_timestep();

						simVars.sim.CFL = leapfrog_original_cfl;
						simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

						if (num_timesteps > 1)
							FatalError("Doesn't make sense because the previous time step is half the time step size in advance");
					}
					else
					{
						run_timestep();
					}

					for (int i = 1; i < num_timesteps; i++)
					{
//						std::cout << "Timestep " << i << std::endl;
						run_timestep();
					}

					if (simVars.pde.id == 2)
					{
						/*
						 * This is the vort/div formulation which uses the velocities as
						 * prognostic values
						 *
						 * We assume that the delta function was applied to the vorticity/divergence
						 * and convert this to the velocity
						 */
						if (simVars.misc.sphere_use_robert_functions)
							op.robert_uv_to_vortdiv(prog_phys_u, prog_phys_v, prog_u, prog_v);
						else
							op.uv_to_vortdiv(prog_phys_u, prog_phys_v, prog_u, prog_v);
					}

					if (simVars.disc.normal_mode_analysis_generation == 3)
					{
						/*
						 * Compute
						 *    1/dt * (U(t+1) - U(t))
						 * for linearization
						 */
						prog[outer_prog_id]->request_data_spectral();
						prog[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

						for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
							prog[inner_prog_id]->operator*=(1.0/simVars.timecontrol.current_timestep_size);
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_spectral();

						// eliminate shift for zero mode since this is non-sense information
						// TODO: Are we really allowed to do this?
						for (int n = 0; n <= sphereDataConfig->spectral_modes_n_max; n++)
						{
							if (std::abs(prog[inner_prog_id]->spectral_space_data[n].imag()) > 1e-10)
							{
								std::cerr << "Phase shift at mode " << n << " too large: " << prog[inner_prog_id]->spectral_space_data[n].imag() << std::endl;
								std::cerr << "This is not an error" << std::endl;
								FatalError("Phase shift");
							}
						}

						for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].real();
							file << "\t";
						}
					}

					for (int inner_prog_id = 0; inner_prog_id < max_prog_id; inner_prog_id++)
					{
						prog[inner_prog_id]->request_data_spectral();

						for (int k = 0; k < sphereDataConfig->spectral_array_data_number_of_elements; k++)
						{
							file << prog[inner_prog_id]->spectral_space_data[k].imag();

							if (inner_prog_id != max_prog_id-1 || k != sphereDataConfig->spectral_array_data_number_of_elements-1)
								file << "\t";
							else
								file << std::endl;
						}
					}
				}
			}
		}
	}



	void run_timestep()
	{

#if SWEET_GUI
			if (simVars.misc.gui_enabled)
				timestep_check_output();
#endif

#if 0
		std::cout
			<< prog_h.physical_reduce_min() << ", " << prog_h.physical_reduce_max() << "   "
			<< prog_u.physical_reduce_min() << ", " << prog_u.physical_reduce_max() << "   "
			<< prog_v.physical_reduce_min() << ", " << prog_v.physical_reduce_max() << std::endl;
#endif

		// output of time step size
		double o_dt;

		if (simVars.disc.timestepping_method == simVars.disc.RUNGE_KUTTA_EXPLICIT)
		{
			switch (simVars.pde.id)
			{
			case 0:
				// SWE, UV
				spheredata_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_advective_proguv,	///< pointer to function to compute euler time step updates
						prog_h, prog_u, prog_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 1:
				// SWE, VORT/DIV
				spheredata_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_vectorinvariant_progvortdiv,	///< pointer to function to compute euler time step updates
						prog_phi, prog_vort, prog_div,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;


			case 2:
			{
				// Advective formulation
				spheredataphysical_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_vectorinvariant_proguv,	///< pointer to function to compute euler time step updates
						prog_phys_phi, prog_phys_u, prog_phys_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);

				prog_phi = prog_phys_phi;
				prog_u = prog_phys_u;
				prog_v = prog_phys_v;
			}
				break;

			case 12:
				spheredata_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 13:
				spheredata_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection_div_free,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 14:
				// ADVECTION, VORT/DIV
				spheredata_timestepping_explicit_rk.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection_vortdiv,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;
			}
		}
		else if (simVars.disc.timestepping_method == simVars.disc.LEAPFROG_EXPLICIT)
		{
			switch (simVars.pde.id)
			{
			case 0:
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_advective_proguv,	///< pointer to function to compute euler time step updates
						prog_h, prog_u, prog_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 1:
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_vectorinvariant_progvortdiv,	///< pointer to function to compute euler time step updates
						prog_phi, prog_vort, prog_div,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 2:
				FatalError("TODO");
#if 0
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_swe_vectorinvariant_proguv,	///< pointer to function to compute euler time step updates
						prog_phi, prog_u, prog_v,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
#endif
				break;

			case 12:
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 13:
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection_div_free,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;

			case 14:
				spheredata_timestepping_explicit_leapfrog.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update_advection_vortdiv,	///< pointer to function to compute euler time step updates
						prog_h,
						o_dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
				break;
			}
		}
		else if (simVars.disc.timestepping_method == simVars.disc.IMPLICIT_TIMESTEP || simVars.disc.timestepping_method == simVars.disc.CRANK_NICOLSON)
		{
			SphereData o_prog_h(sphereDataConfig);
			SphereData o_prog_u(sphereDataConfig);
			SphereData o_prog_v(sphereDataConfig);

			o_dt = -simVars.sim.CFL;

			// padding to max simulation time if exceeding the maximum
			if (simVars.timecontrol.max_simulation_time >= 0)
			{
				if (o_dt + simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
				{
					o_dt = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

					std::cout << "WARNING: IMPLICIT TS SETUP CALLED DURING SIMULATION TIME FRAME" << std::endl;
					std::cerr << "WARNING: IMPLICIT TS SETUP CALLED DURING SIMULATION TIME FRAME" << std::endl;

					timestepping_implicit_swe.setup(
							sphereDataConfig,
							sphereDataConfigExt,

							simVars.sim.earth_radius,
							simVars.sim.coriolis_omega,
							simVars.sim.gravitation*simVars.sim.h0,
							o_dt,

							param_use_coriolis_formulation,

							simVars.disc.timestepping_order
						);
				}
			}

			switch (simVars.pde.id)
			{
			case 0:
				timestepping_implicit_swe.solve(
						prog_h*simVars.sim.gravitation, prog_u, prog_v,
						o_prog_h, o_prog_u, o_prog_v,
						o_dt
					);

				prog_h = o_prog_h / simVars.sim.gravitation;
				prog_u = o_prog_u;
				prog_v = o_prog_v;

				break;

			case 12:
				timestepping_implicit_swe.solve_advection(
						prog_h*simVars.sim.gravitation, prog_u, prog_v,
						o_prog_h,
						o_dt
					);

				prog_h = o_prog_h / simVars.sim.gravitation;
				break;

			default:
				FatalError("PDE not supported for implicit TS");
			}
		}
		else if (simVars.disc.timestepping_method == simVars.disc.REXI)
		{
			o_dt = simVars.timecontrol.current_timestep_size;
			assert(o_dt > 0);

			// padding to max simulation time if exceeding the maximum
			if (simVars.timecontrol.max_simulation_time >= 0)
			{
				if (o_dt + simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
				{
					o_dt = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

					std::cout << "WARNING: REXI SETUP CALLED DURING SIMULATION TIME FRAME" << std::endl;
					std::cerr << "WARNING: REXI SETUP CALLED DURING SIMULATION TIME FRAME" << std::endl;

					swe_sphere_rexi.setup(
							simVars.rexi.rexi_h,
							simVars.rexi.rexi_M,
							simVars.rexi.rexi_L,

							sphereDataConfig,
							&simVars.sim,
							o_dt,

							simVars.rexi.rexi_use_half_poles,
							simVars.misc.sphere_use_robert_functions,
							simVars.pde.id,

							simVars.rexi.rexi_use_extended_modes,
							simVars.rexi.rexi_normalization,

							param_use_coriolis_formulation,
							simVars.sim.f_sphere,
							simVars.rexi.rexi_sphere_solver_preallocation
						);
				}
			}


			switch (simVars.pde.id)
			{
			case 0:
				// REXI time stepping
				prog_h *= simVars.sim.gravitation;
				swe_sphere_rexi.run_timestep_rexi_advection_progphiuv(
						prog_h,	// phi
						prog_u,
						prog_v,
						o_dt,
						simVars
					);
				prog_h *= (1.0/simVars.sim.gravitation);
				break;


			case 1:
				// REXI time stepping
				prog_h *= simVars.sim.gravitation;
				swe_sphere_rexi.run_timestep_rexi_vectorinvariant_progphivortdiv(
						prog_h,	// phi
						prog_vort,
						prog_div,
						o_dt,
						simVars
					);
				prog_h *= (1.0/simVars.sim.gravitation);
				break;

			default:
				FatalError("PDE id not yet implemented");
				break;
			}



			/*
			 * Add implicit viscosity
			 */
			if (simVars.sim.viscosity != 0)
			{
				double scalar = simVars.sim.viscosity*simVars.timecontrol.current_timestep_size;
				double r = simVars.sim.earth_radius;

				/*
				 * (1-dt*visc*D2)p(t+dt) = p(t)
				 */
				prog_h = prog_h.spectral_solve_helmholtz(1.0, -scalar, r);
				prog_u = prog_u.spectral_solve_helmholtz(1.0, -scalar, r);
				prog_v = prog_v.spectral_solve_helmholtz(1.0, -scalar, r);
			}
		}
		else if (simVars.disc.timestepping_method == simVars.disc.CRANK_NICOLSON)
		{
		}
		else
		{
			FatalError("Timestepping method is not supported!");
		}


		// advance time step and provide information to parameters
		simVars.timecontrol.current_timestep_size = o_dt;
		simVars.timecontrol.current_simulation_time += o_dt;
		simVars.timecontrol.current_timestep_nr++;
	}



	/**
	 * Euler time step for advection along the longitude
	 */
	void p_run_euler_timestep_update_advection(
			const SphereData &i_h,	///< prognostic variables
			SphereData &o_h_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (simVars.misc.sphere_use_robert_functions)
		{
			o_h_t =
				-(
					op.robert_div_lon(prog_u*i_h)+
					op.robert_div_lat(prog_v*i_h)
				)*(1.0/simVars.sim.earth_radius);
		}
		else
		{
			o_h_t = -(
					op.div_lon(prog_u*i_h)+
					op.div_lat(prog_v*i_h)
				)*(1.0/simVars.sim.earth_radius);
		}


		assert(simVars.sim.viscosity_order == 2 || simVars.sim.viscosity_order == 4);
		if (simVars.sim.viscosity != 0)
		{
			if (simVars.sim.viscosity_order == 2)
			{
				double scalar = simVars.sim.viscosity;
				o_h_t += op.laplace(i_h)*scalar;
			}
			else if (simVars.sim.viscosity_order == 4)
			{
				double scalar = simVars.sim.viscosity;
				o_h_t += op.laplace(op.laplace(i_h))*scalar;
			}
		}
	}


	void p_run_euler_timestep_update_advection_div_free(
			const SphereData &i_h,	///< prognostic variables
			SphereData &o_h_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		if (simVars.misc.use_nonlinear_equations)
			FatalError("Advection equation is only possible without non-linearities and with robert functions");

		if (simVars.misc.sphere_use_robert_functions)
			o_h_t = -(prog_u*op.robert_grad_lon_M(i_h) + prog_v*op.robert_grad_lat_M(i_h))*(1.0/simVars.sim.earth_radius);
		else
			o_h_t = -(prog_u*op.grad_lon(i_h)+prog_v*op.grad_lat(i_h))*(1.0/simVars.sim.earth_radius);

		assert(simVars.sim.viscosity_order == 2 || simVars.sim.viscosity_order == 4);
		if (simVars.sim.viscosity != 0)
		{
			if (simVars.sim.viscosity_order == 2)
			{
				o_h_t += op.laplace(i_h)*simVars.sim.viscosity;
			}
			else if (simVars.sim.viscosity_order == 4)
			{
				o_h_t += op.laplace(op.laplace(i_h))*simVars.sim.viscosity;
			}
		}
	}




	/**
	 * Euler time step for advection along the longitude
	 */
	void p_run_euler_timestep_update_advection_vortdiv(
			const SphereData &i_h,	///< prognostic variables
			SphereData &o_h_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;


		if (simVars.misc.sphere_use_robert_functions)
		{
			SphereDataPhysical ug = prog_u.getSphereDataPhysical();
			SphereDataPhysical vg = prog_v.getSphereDataPhysical();

			SphereDataPhysical phig = i_h.getSphereDataPhysical();

			SphereData tmpspec(sphereDataConfig);
			op.robert_uv_to_vortdiv(ug*phig,vg*phig, tmpspec, o_h_t);

			o_h_t *= -1.0;
		}
		else
		{
			SphereDataPhysical ug = prog_u.getSphereDataPhysical();
			SphereDataPhysical vg = prog_v.getSphereDataPhysical();

			SphereDataPhysical phig = i_h.getSphereDataPhysical();

			SphereData tmpspec(sphereDataConfig);
			op.uv_to_vortdiv(ug*phig,vg*phig, tmpspec, o_h_t);

			o_h_t *= -1.0;
		}


		assert(simVars.sim.viscosity_order == 2 || simVars.sim.viscosity_order == 4);
		if (simVars.sim.viscosity != 0)
		{
			if (simVars.sim.viscosity_order == 2)
			{
				o_h_t += op.laplace(i_h)*simVars.sim.viscosity;
			}
			else if (simVars.sim.viscosity_order == 4)
			{
				o_h_t += op.laplace(op.laplace(i_h))*simVars.sim.viscosity;
			}
		}
	}


	SphereDataPhysical add_f(const SphereDataPhysical &i_input)
	{
		SphereDataPhysical out(i_input);

		out.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data += mu*2.0*simVars.sim.coriolis_omega;
			}
		);
		return out;
	}



	/*
	 * Shallow water time stepping
	 * (Single stage realized with Euler)
	 */
	void p_run_euler_timestep_update_swe_advective_proguv(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{

		o_dt = simVars.timecontrol.current_timestep_size;


		/*
		 * NON-LINEARITIES
		 */
		assert(simVars.sim.earth_radius > 0);
		assert(simVars.sim.gravitation);

		if (simVars.misc.sphere_use_robert_functions)
		{
			/*
			 * ROBERT
			 */
			if (simVars.misc.use_nonlinear_equations)
			{

				/*
				 * NON-LINEAR
				 */
				SphereData o_phi_t(i_h.sphereDataConfig);

				double i_r = 1.0/simVars.sim.earth_radius;

				SphereData i_phi = i_h*simVars.sim.gravitation;

				o_phi_t =
#if 1
						-i_r*op.robert_div(i_phi*i_u, i_phi*i_v)
#else
						-i_r*i_phi*op.robert_div(i_u, i_v)
						-i_r*op.robert_grad_M(i_phi, i_u, i_v)
#endif
						;

				o_u_t = -i_r*op.ritchie_A(i_u, i_v, i_u)
						-i_r*op.robert_grad_lon(i_phi)
						;

				o_v_t = -i_r*op.ritchie_A(i_u, i_v, i_v)
						-i_r*op.mu(i_u*i_u + i_v*i_v)
						-i_r*op.robert_grad_lat(i_phi)
						;

				if (simVars.sim.coriolis_omega != 0 && param_use_coriolis_formulation)
				{
					if (simVars.sim.f_sphere)
					{
						o_u_t += simVars.sim.f0*i_v;
						o_v_t -= simVars.sim.f0*i_u;
					}
					else
					{
						o_u_t += f(i_v);
						o_v_t -= f(i_u);
					}
				}

				o_h_t = o_phi_t*(1.0/simVars.sim.gravitation);

				o_h_t = -(op.robert_div_lon(i_h*i_u)+op.robert_div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);

			}
			else
			{
				/*
				 * LINEAR
				 */

				o_h_t =
					-(
						op.robert_div_lon(i_u)+
						op.robert_div_lat(i_v)
					)*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t =
						-op.robert_grad_lon(i_h)*
						(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t =
						-op.robert_grad_lat(i_h)*
						(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0 && param_use_coriolis_formulation)
				{
					// This creates an instability!!!
#if 0

					// This doesn't solve the instability
					SphereDataPhysical vg = i_v.getSphereDataPhysical();
					SphereDataPhysical ug = i_u.getSphereDataPhysical();

					o_u_t = o_u_t.getSphereDataPhysical() + vg*fg;
					o_v_t = o_v_t.getSphereDataPhysical() - ug*fg;
#elif 0

					SphereDataPhysical vg = i_v.getSphereDataPhysical();
					SphereDataPhysical ug = i_u.getSphereDataPhysical();

					SphereDataPhysical tmpg1 = ug*(fg);
					SphereDataPhysical tmpg2 = vg*(fg);
					SphereData o_div_t(sphereDataConfig);
					SphereData o_vort_t(sphereDataConfig);

					op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);
					o_vort_t *= -1.0;

					SphereData null(sphereDataConfig);
					null.spectral_set_zero();
					op.robert_vortdiv_to_uv(o_vort_t, null, ug, vg);

					o_u_t += ug;
					o_v_t += vg;
#else
					if (simVars.sim.f_sphere)
					{
						o_u_t += simVars.sim.f0*i_v;
						o_v_t -= simVars.sim.f0*i_u;
					}
					else
					{
						o_u_t += f(i_v);
						o_v_t -= f(i_u);
					}
#endif
				}
			}
		}
		else
		{
			/*
			 * NON-ROBERT
			 */
			if (simVars.misc.use_nonlinear_equations)
			{
				/*
				 * NON-LINEAR
				 */
				o_h_t = -(op.div_lon(i_h*i_u)+op.div_lat(i_h*i_v))*(1.0/simVars.sim.earth_radius);
				o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					if (simVars.sim.f_sphere)
					{
						o_u_t += simVars.sim.f0*i_v;
						o_v_t -= simVars.sim.f0*i_u;
					}
					else
					{
						o_u_t += f(i_v);
						o_v_t -= f(i_u);
					}
				}

				o_u_t -= (i_u*op.grad_lon(i_u) + i_v*op.grad_lat(i_u))*(1.0/simVars.sim.earth_radius);
				o_v_t -= (i_u*op.grad_lon(i_v) + i_v*op.grad_lat(i_v))*(1.0/simVars.sim.earth_radius);
			}
			else
			{
				/*
				 * LINEAR
				 */
				o_h_t = -(op.div_lon(i_u)+op.div_lat(i_v))*(simVars.sim.h0/simVars.sim.earth_radius);

				o_u_t = -op.grad_lon(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);
				o_v_t = -op.grad_lat(i_h)*(simVars.sim.gravitation/simVars.sim.earth_radius);

				if (simVars.sim.coriolis_omega != 0)
				{
					if (simVars.sim.f_sphere)
					{
						o_u_t += simVars.sim.f0*i_v;
						o_v_t -= simVars.sim.f0*i_u;
					}
					else
					{
						o_u_t += f(i_v);
						o_v_t -= f(i_u);
					}
				}

			}
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			o_h_t += op.laplace(i_h)*simVars.sim.viscosity;
			o_u_t += op.laplace(i_u)*simVars.sim.viscosity;
			o_v_t += op.laplace(i_v)*simVars.sim.viscosity;
		}
	}




	/*
	 * Shallow water time stepping
	 * (Single stage realized with Euler)
	 * This formulation uses the vorticity/divergence formulation
	 */
	void p_run_euler_timestep_update_swe_vectorinvariant_progvortdiv(
			const SphereData &i_phispec,		///< prognostic variables
			const SphereData &i_vortspec,	///< prognostic variables
			const SphereData &i_divspec,	///< prognostic variables

			SphereData &o_phi_t,	///< time updates
			SphereData &o_vort_t,	///< time updates
			SphereData &o_div_t,	///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		o_dt = simVars.timecontrol.current_timestep_size;

		/*
		 * NON-LINEARITIES
		 */
		assert(simVars.sim.earth_radius > 0);
		assert(simVars.sim.gravitation);

#if 0
		if (i_simulation_timestamp == 0)
		{
			i_phispec.physical_file_write("o_phi_0.csv");
			i_vortspec.physical_file_write("o_vort_0.csv");
			i_divspec.physical_file_write("o_div_0.csv");
		}
#endif

		if (simVars.misc.sphere_use_robert_functions)
		{
			/*
			 * ROBERT
			 */
			if (simVars.misc.use_nonlinear_equations)
			{
				/*
				 * NON-LINEAR
				 *
				 * Follows Hack & Jakob formulation
				 */

				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.robert_vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*(vrtg+fg);
				SphereDataPhysical tmpg2 = vg*(vrtg+fg);

				op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*phig;
				tmpg2 = vg*phig;

				SphereData tmpspec(sphereDataConfig);
				op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

				o_phi_t *= -1.0;

				tmpspec = (phig+0.5*(ug*ug+vg*vg));
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);
			}
			else
			{
				/*
				 * LINEAR
				 */
				if (simVars.sim.f_sphere)
				{
					double gh = simVars.sim.gravitation * simVars.sim.h0;

					o_phi_t = -gh*i_divspec;
					o_div_t = -op.laplace(i_phispec);

					if (simVars.sim.coriolis_omega != 0 && param_use_coriolis_formulation)
					{
						o_vort_t = -f(i_divspec);
						o_div_t += f(i_vortspec);
					}
					else
					{
						o_vort_t.spectral_set_zero();
					}
				}
				else
				{

					double gh = simVars.sim.gravitation * simVars.sim.h0;

					SphereDataPhysical ug(i_phispec.sphereDataConfig);
					SphereDataPhysical vg(i_phispec.sphereDataConfig);

	//				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
	//				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
					op.robert_vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
					SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

					SphereDataPhysical tmpg1 = ug*fg;
					SphereDataPhysical tmpg2 = vg*fg;

					op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

					o_vort_t *= -1.0;

					tmpg1 = ug*gh;
					tmpg2 = vg*gh;

					SphereData tmpspec(sphereDataConfig);
					op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

					o_phi_t *= -1.0;

					tmpspec = phig;
					tmpspec.request_data_spectral();
					o_div_t += -op.laplace(tmpspec);
				}
			}
		}
		else
		{
			if (simVars.misc.use_nonlinear_equations)
			{
				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*(vrtg+fg);
				SphereDataPhysical tmpg2 = vg*(vrtg+fg);

				op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*phig;
				tmpg2 = vg*phig;

				SphereData tmpspec(sphereDataConfig);
				op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

				o_phi_t *= -1.0;

				tmpspec = (phig+0.5*(ug*ug+vg*vg));
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);

			}
			else
			{
				double gh = simVars.sim.gravitation * simVars.sim.h0;

				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

//				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
//				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*fg;
				SphereDataPhysical tmpg2 = vg*fg;

				op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

//				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*gh;
				tmpg2 = vg*gh;

				SphereData tmpspec(sphereDataConfig);
				op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

				o_phi_t *= -1.0;

				tmpspec = phig;
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);
			}
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			o_phi_t += op.laplace(i_phispec)*simVars.sim.viscosity;
			o_vort_t += op.laplace(i_vortspec)*simVars.sim.viscosity;
			o_div_t += op.laplace(i_divspec)*simVars.sim.viscosity;
		}
	}



	/*
	 * Shallow water time stepping
	 * This formulation uses the vorticity/divergence formulation but the velocities as prognostic variables
	 */
	void p_run_euler_timestep_update_swe_vectorinvariant_proguv(
			const SphereDataPhysical &i_phi,	///< prognostic variables
			const SphereDataPhysical &i_u,		///< prognostic variables
			const SphereDataPhysical &i_v,		///< prognostic variables

			SphereDataPhysical &o_phi_t,	///< time updates
			SphereDataPhysical &o_u_t,		///< time updates
			SphereDataPhysical &o_v_t,		///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		SphereData i_vortspec(i_phi.sphereDataConfig);
		SphereData i_divspec(i_phi.sphereDataConfig);

		SphereData o_vort_t(i_phi.sphereDataConfig);
		SphereData o_div_t(i_phi.sphereDataConfig);

		SphereData i_phispec = i_phi;
		SphereData o_phispec_t(i_phi.sphereDataConfig);

		if (simVars.misc.sphere_use_robert_functions)
			op.robert_uv_to_vortdiv(i_u, i_v, i_vortspec, i_divspec);
		else
			op.uv_to_vortdiv(i_u, i_v, i_vortspec, i_divspec);

		o_dt = simVars.timecontrol.current_timestep_size;

		/*
		 * NON-LINEARITIES
		 */
		assert(simVars.sim.earth_radius > 0);
		assert(simVars.sim.gravitation);


		if (simVars.misc.sphere_use_robert_functions)
		{
			/*
			 * ROBERT
			 */
			if (simVars.misc.use_nonlinear_equations)
			{
				FatalError("TODO");
				/*
				 * NON-LINEAR
				 *
				 * Follows Hack & Jakob formulation
				 */

				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.robert_vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*(vrtg+fg);
				SphereDataPhysical tmpg2 = vg*(vrtg+fg);

				op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*phig;
				tmpg2 = vg*phig;

				SphereData tmpspec(sphereDataConfig);
				op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phispec_t);

				o_phispec_t *= -1.0;

				tmpspec = (phig+0.5*(ug*ug+vg*vg));
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);
			}
			else
			{
				/*
				 * LINEAR
				 */
				double gh = simVars.sim.gravitation * simVars.sim.h0;

#if 1
				SphereDataPhysical tmpg1 = i_u*fg;
				SphereDataPhysical tmpg2 = i_v*fg;

				op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				tmpg1 = i_u*gh;
				tmpg2 = i_v*gh;

				SphereData tmpspec(sphereDataConfig);
				op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phispec_t);

				o_phispec_t *= -1.0;

				tmpspec = i_phi;
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);

#else

				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

//				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
//				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.robert_vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*fg;
				SphereDataPhysical tmpg2 = vg*fg;

				op.robert_uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

//				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*gh;
				tmpg2 = vg*gh;

				SphereData tmpspec(sphereDataConfig);
				op.robert_uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phispec_t);

				o_phispec_t *= -1.0;

				tmpspec = phig;
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);
#endif
			}
		}
		else
		{
			if (simVars.misc.use_nonlinear_equations)
			{
				SphereDataPhysical ug(i_phispec.sphereDataConfig);
				SphereDataPhysical vg(i_phispec.sphereDataConfig);

				SphereDataPhysical vrtg = i_vortspec.getSphereDataPhysical();
				SphereDataPhysical divg = i_divspec.getSphereDataPhysical();
				op.vortdiv_to_uv(i_vortspec, i_divspec, ug, vg);
				SphereDataPhysical phig = i_phispec.getSphereDataPhysical();

				SphereDataPhysical tmpg1 = ug*(vrtg+fg);
				SphereDataPhysical tmpg2 = vg*(vrtg+fg);

				op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				SphereDataPhysical tmpg = o_div_t.getSphereDataPhysical();

				tmpg1 = ug*phig;
				tmpg2 = vg*phig;

				SphereData tmpspec(sphereDataConfig);
				op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phispec_t);

				o_phispec_t *= -1.0;

				tmpspec = (phig+0.5*(ug*ug+vg*vg));
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);

			}
			else
			{
				double gh = simVars.sim.gravitation * simVars.sim.h0;

				SphereDataPhysical tmpg1 = i_u*fg;
				SphereDataPhysical tmpg2 = i_v*fg;

				op.uv_to_vortdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

				o_vort_t *= -1.0;

				tmpg1 = i_u*gh;
				tmpg2 = i_v*gh;

				SphereData tmpspec(sphereDataConfig);
				op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_phispec_t);

				o_phispec_t *= -1.0;

				tmpspec = i_phi;
				tmpspec.request_data_spectral();
				o_div_t += -op.laplace(tmpspec);
			}
		}

		assert(simVars.sim.viscosity_order == 2);
		if (simVars.sim.viscosity != 0)
		{
			o_phispec_t += op.laplace(i_phispec)*simVars.sim.viscosity;
			o_vort_t += op.laplace(i_vortspec)*simVars.sim.viscosity;
			o_div_t += op.laplace(i_divspec)*simVars.sim.viscosity;
		}

		if (simVars.misc.sphere_use_robert_functions)
			op.robert_vortdiv_to_uv(o_vort_t, o_div_t, o_u_t, o_v_t);
		else
			op.vortdiv_to_uv(o_vort_t, o_div_t, o_u_t, o_v_t);

		o_phi_t = o_phispec_t.getSphereDataPhysical();
	}



#if SWEET_GUI

	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data
	)
	{
		// request rendering of sphere
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		int id = simVars.misc.vis_id % 4;
		switch (id)
		{
		default:
		case 0:
			// USE COPY TO AVOID FORWARD/BACKWARD TRANSFORMATION
			switch (simVars.pde.id)
			{
			case 1:
			case 2:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_phi)/simVars.sim.gravitation, planeDataConfig);
				break;

			default:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_h), planeDataConfig);
			}
			break;

		case 1:
			switch (simVars.pde.id)
			{
			case 1:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_vort), planeDataConfig);
				break;

			default:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(
						SphereData(prog_u),
						planeDataConfig);
			}
			break;

		case 2:
			switch (simVars.pde.id)
			{
			case 1:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(SphereData(prog_div), planeDataConfig);
				break;

			default:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(
						SphereData(prog_v),
						planeDataConfig);
			}
			break;

		case 3:
			switch (simVars.pde.id)
			{
			case 1:
				viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(prog_vort, planeDataConfig);
				break;

			default:
				if (simVars.misc.sphere_use_robert_functions)
					viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(op.vort(SphereData(prog_u), SphereData(prog_v)), planeDataConfig);
				else
					viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert(op.robert_vort(SphereData(prog_u), SphereData(prog_v)), planeDataConfig);
			}
			break;

		}


#if 0
		{
			SphereData test_h(sphereDataConfig);
			SphereData test_u(sphereDataConfig);
			SphereData test_v(sphereDataConfig);

			SphereBenchmarksCombined::setupInitialConditions(test_h, test_u, test_v, simVars, op);


			SphereData h = prog_phi*(1.0/simVars.sim.gravitation);
			SphereDataPhysical u(sphereDataConfig);
			SphereDataPhysical v(sphereDataConfig);

			if (simVars.misc.sphere_use_robert_functions)
				op.robert_vortdiv_to_uv(prog_vort, prog_div, u, v);
			else
				op.vortdiv_to_uv(prog_vort, prog_div, u, v);

			SphereDataPhysical test_ug = test_u.getSphereDataPhysical();
			SphereDataPhysical test_vg = test_v.getSphereDataPhysical();

			viz_plane_data = Convert_SphereData_To_PlaneData::physical_convert((h-test_h), planeDataConfig);

//			error_h = h.physical_reduce_max_abs(test_h);
		}
#endif

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
//		update_diagnostics();

		const char* description = "";

		int id = simVars.misc.vis_id % 4;

		switch(simVars.pde.id)
		{
		case 1:
			switch (id)
			{
			default:
			case 0:
				description = "H";
				break;

			case 1:
				description = "vort";
				break;

			case 12:
				description = "div";
				break;

			case 13:
				description = "eta";
				break;
			}
			break;

		default:
			switch (id)
			{
			default:
			case 0:
				description = "H";
				break;

			case 1:
				description = "U";
				break;

			case 2:
				description = "V";
				break;

			case 3:
				description = "eta";
				break;
			}
			break;
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string, "Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
				viz_plane_data.reduce_max(),
				viz_plane_data.reduce_min()
		);

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
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
			break;

		case 'c':
			write_file_output();
			break;

#if 0
case 'C':
			// dump data arrays to VTK
			prog_h.file_physical_saveData_vtk("swe_rexi_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_rexi_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_rexi_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h.file_physical_loadData("swe_rexi_dump_h.csv", simVars.setup.input_data_binary);
			prog_u.file_physical_loadData("swe_rexi_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_physical_loadData("swe_rexi_dump_v.csv", simVars.setup.input_data_binary);
			break;
#endif
		}
	}
#endif
};


int main(int i_argc, char *i_argv[])
{

#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
		{
				std::cerr << "MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye..." << std::endl;
				exit(-1);
		}
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"use-coriolis-formulation",
			"compute-error",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = param_use_coriolis_formulation;
	simVars.bogus.var[1] = param_compute_error;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		std::cout << "	--use-coriolis-formulation [0/1]	Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
		std::cout << "	--compute-error [0/1]	Output errors (if available, default: 1)" << std::endl;
		std::cout << "	--pde-id [0/1]	PDE ID (0: SWE, 1: Advection, 2: Advection divergence free)" << std::endl;
		return -1;
	}


	param_use_coriolis_formulation = simVars.bogus.var[0];
	assert(param_use_coriolis_formulation == 0 || param_use_coriolis_formulation == 1);
//	param_compute_error = simVars.bogus.var[1];

	sphereDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);
/*
	sphereDataConfigInstance.setupAutoPhysicalSpace(
					simVars.disc.res_spectral[0],
					simVars.disc.res_spectral[1],
					&simVars.disc.res_physical[0],
					&simVars.disc.res_physical[1]
			);
*/
	sphereDataConfigExtInstance.setupAdditionalModes(
			&sphereDataConfigInstance,
			simVars.rexi.rexi_use_extended_modes,
			simVars.rexi.rexi_use_extended_modes
		);

#if SWEET_GUI
	planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.res_physical);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);


#if SWEET_MPI

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	std::cout << "MPI RANK: " << mpi_rank << std::endl;

	// only start simulation and time stepping for first rank
	if (mpi_rank == 0)
#endif

	{
		std::cout << "SPH config string: " << sphereDataConfigInstance.getConfigInformationString() << std::endl;

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

#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationSWE);
			delete simulationSWE;
		}
		else
#endif
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			// reset already triggered
			//simulationSWE->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			//double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
//				simulationSWE->update_diagnostics();
#if 0
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
#endif
			}

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			//Start counting time
			time.reset();


			if (simVars.disc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normal_mode_analysis();
			}
			else
			{
				bool output_written = false;

				// Main time loop
				while(true)
				{
					output_written = simulationSWE->timestep_check_output();

					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simVars.misc.stability_checks)
					{
						if (simulationSWE->detect_instability())
						{
							std::cout << "INSTABILITY DETECTED" << std::endl;
							break;
						}
					}
				}

				// Output final time step output!
				if (!output_written)
					simulationSWE->timestep_do_output();
			}

			// Stop counting time
			time.stop();

			double seconds = time();

			// End of run output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;

			delete simulationSWE;
		}
	}
#if SWEET_MPI
	else
	{
		if (simVars.disc.timestepping_method == 100)
		{
			SWE_Sphere_REXI swe_sphere_rexi;

			/*
			 * Setup our little dog REXI
			 */
			swe_sphere_rexi.setup(
					simVars.rexi.rexi_h,
					simVars.rexi.rexi_M,
					simVars.rexi.rexi_L,

					sphereDataConfig,
					&simVars.sim,
					simVars.timecontrol.current_timestep_size,

					simVars.rexi.rexi_use_half_poles,
					simVars.misc.sphere_use_robert_functions,
					simVars.pde.id,

					simVars.rexi.rexi_use_extended_modes,
					simVars.rexi.rexi_normalization,

					param_use_coriolis_formulation,
					simVars.sim.f_sphere,
					simVars.rexi.rexi_sphere_solver_preallocation
				);

			bool run = true;

			SphereData prog_h(sphereDataConfig);
			SphereData prog_u(sphereDataConfig);
			SphereData prog_v(sphereDataConfig);

			// initialize with dummy data
			prog_h.spectral_set_zero();
			prog_u.spectral_set_zero();
			prog_v.spectral_set_zero();

			MPI_Barrier(MPI_COMM_WORLD);

			while (run)
			{
				// REXI time stepping
				prog_h *= simVars.sim.gravitation;
				run = swe_sphere_rexi.run_timestep_rexi(
						prog_h, prog_u, prog_v,
						-simVars.sim.CFL,
						simVars
				);
				prog_h *= (1.0/simVars.sim.gravitation);
			}
		}
	}
#endif


#if SWEET_MPI
	if (simVars.disc.timestepping_method == 100)
	{
		// synchronize REXI
		if (mpi_rank == 0)
			SWE_Sphere_REXI::MPI_quitWorkers(sphereDataConfig);
	}

	MPI_Finalize();
#endif

	return 0;
}
