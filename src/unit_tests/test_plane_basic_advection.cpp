
#include "../include/sweet/plane/PlaneData_Spectral.hpp"
#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include "../programs/swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include "../include/sweet/plane/PlaneOperators.hpp"
#include <sweet/Stopwatch.hpp>

#include <sweet/ProgramArguments.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;


SimulationVariables _simVars;

#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>


/**
 * This class stored the discretization-related parameters
 *
 * resolution / timestepping
 */
class ShackProgramSpecific	:
		public sweet::ClassDictionaryInterface
{
public:
	double velocity_u = 0;
	double velocity_v = 0;
	int advection_scheme = 0;
	bool staggered_use_analytical_solution = 0;
	bool test_mode = 0;


public:
	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << "Program-specific options:" << std::endl;
		std::cout << "	--velocity-u [velocity in u direction]" << std::endl;
		std::cout << std::endl;
		std::cout << "	--velocity-v [velocity in v direction]" << std::endl;
		std::cout << std::endl;
		std::cout << "	--advection-scheme [nr]		Advection scheme" << std::endl;
		std::cout << "	                            0: up/downwinding" << std::endl;
		std::cout << "	                            1: staggered" << std::endl;
		std::cout << "	                            2: non-staggered" << std::endl;
		std::cout << "	                            3: no h update" << std::endl;
		std::cout << std::endl;
		std::cout << "	--staggered-use-analytical-solution [0/1]" << std::endl;
		std::cout << "	                            Use analytical solution for non-staggered advection" << std::endl;
		std::cout << std::endl;
		std::cout << "	--test-mode [nr]	Test mode" << std::endl;
		std::cout << "	                    0: space" << std::endl;
		std::cout << "	                    1: time" << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--velocity-u", velocity_u);
		i_pa.getArgumentValueByKey("--velocity-v", velocity_v);
		i_pa.getArgumentValueByKey("--advection-scheme", advection_scheme);
		i_pa.getArgumentValueByKey("--staggered-use-analytical-solution", staggered_use_analytical_solution);
		i_pa.getArgumentValueByKey("--test-mode", test_mode);

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "ShackProgramSpecific:" << std::endl;
		std::cout << " + velocity_u: " << velocity_u << std::endl;
		std::cout << " + velocity_v: " << velocity_v << std::endl;
		std::cout << " + advection_scheme: " << advection_scheme << std::endl;
		std::cout << " + staggered_use_analytical_solution: " << staggered_use_analytical_solution << std::endl;
		std::cout << " + test_mode: " << test_mode << std::endl;
		std::cout << std::endl;
	}
};

ShackProgramSpecific _shackProgramSpecific;

//
// 0: Gaussian (WARNING! DON'T USE THIS AS INITIAL CONDITIONS!)
// 1: sin curve
//
#define ADV_FUNCTION	1

class SimulationAdvection
{
public:
	PlaneDataConfig *planeDataConfig;

	PlaneData_Spectral prog_h;
	PlaneData_Spectral prog_u;
	PlaneData_Spectral prog_v;

	PlaneData_Spectral hu;
	PlaneData_Spectral hv;

	PlaneData_Spectral tmp;

	PlaneOperators op;

	PlaneDataTimesteppingExplicitRK timestepping;

#if ADV_FUNCTION==1
	double freq_x = 4.0;
	double freq_y = 4.0;
#endif


public:
	SimulationAdvection(PlaneDataConfig *i_planeDataConfig)	:
		planeDataConfig(i_planeDataConfig),
		prog_h(i_planeDataConfig),
		prog_u(i_planeDataConfig),
		prog_v(i_planeDataConfig),

		hu(i_planeDataConfig),
		hv(i_planeDataConfig),

		tmp(i_planeDataConfig),

		op(	i_planeDataConfig, _simVars.sim.plane_domain_size, _simVars.disc.space_use_spectral_basis_diffs)
	{
		reset();
	}


	PlaneData_Spectral get_advected_solution(
		double i_timestamp
	)
	{
		PlaneData_Spectral ret_h(planeDataConfig);
		PlaneData_Physical ret_h_phys(planeDataConfig);

		double adv_x = -_shackProgramSpecific.velocity_u*i_timestamp;
		double adv_y = -_shackProgramSpecific.velocity_v*i_timestamp;

#if ADV_FUNCTION==0
		double radius = _simVars.benchmark.object_scale*
			std::sqrt(
				 (double)_simVars.sim.plane_domain_size[0]*(double)_simVars.sim.plane_domain_size[0]
				+(double)_simVars.sim.plane_domain_size[1]*(double)_simVars.sim.plane_domain_size[1]
			);
#endif

		ret_h_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{

#if ADV_FUNCTION==0

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				double dx = x-_simVars.benchmark.coord_x*_simVars.sim.plane_domain_size[0];
				double dy = y-_simVars.benchmark.coord_y*_simVars.sim.plane_domain_size[1];

				dx /= radius;
				dy /= radius;

				double value = _simVars.benchmark.h0+std::exp(-50.0*(dx*dx + dy*dy));
				io_data = value;

#elif ADV_FUNCTION==1

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				x /= _simVars.sim.plane_domain_size[0];
				y /= _simVars.sim.plane_domain_size[1];

				io_data = std::sin(freq_x*M_PI*x)*std::sin(freq_x*M_PI*y);
#endif
			}
		);

		ret_h.loadPlaneDataPhysical(ret_h_phys);

		return ret_h;
	}



	PlaneData_Spectral get_advected_solution_diffx(
		double i_timestamp
	)
	{
		PlaneData_Spectral ret_h(planeDataConfig);
		PlaneData_Physical ret_h_phys(planeDataConfig);

		double adv_x = -_shackProgramSpecific.velocity_u*i_timestamp;
		double adv_y = -_shackProgramSpecific.velocity_v*i_timestamp;

#if ADV_FUNCTION==0
		double radius_scale = std::sqrt(
				 (double)_simVars.sim.plane_domain_size[0]*(double)_simVars.sim.plane_domain_size[0]
				+(double)_simVars.sim.plane_domain_size[1]*(double)_simVars.sim.plane_domain_size[1]
			);

		double radius = _simVars.benchmark.object_scale*radius_scale;
#endif


		ret_h_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
#if ADV_FUNCTION==0

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				double dx = x-_simVars.benchmark.coord_x*_simVars.sim.plane_domain_size[0];
				double dy = y-_simVars.benchmark.coord_y*_simVars.sim.plane_domain_size[1];

				dx /= radius;
				dy /= radius;

				double value = -50.0*2.0*dx*std::exp(-50.0*(dx*dx + dy*dy));
				value /= object_scale

				io_data = value;

#elif ADV_FUNCTION==1

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				x /= _simVars.sim.plane_domain_size[0];
				y /= _simVars.sim.plane_domain_size[1];

				io_data = freq_x*M_PI*std::cos(freq_x*M_PI*x)*std::sin(freq_y*M_PI*y)/_simVars.sim.plane_domain_size[0];
#endif
			}
		);

		ret_h.loadPlaneDataPhysical(ret_h_phys);

		return ret_h;
	}



	PlaneData_Spectral get_advected_solution_diffy(
		double i_timestamp
	)
	{
		PlaneData_Spectral ret_h(planeDataConfig);
		PlaneData_Physical ret_h_phys(planeDataConfig);

		double adv_x = -_shackProgramSpecific.velocity_u*i_timestamp;
		double adv_y = -_shackProgramSpecific.velocity_v*i_timestamp;

#if ADV_FUNCTION==0
		double radius_scale = std::sqrt(
				 (double)_simVars.sim.plane_domain_size[0]*(double)_simVars.sim.plane_domain_size[0]
				+(double)_simVars.sim.plane_domain_size[1]*(double)_simVars.sim.plane_domain_size[1]
			);

		double radius = _simVars.benchmark.object_scale*radius_scale;
#endif

		ret_h_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
#if ADV_FUNCTION==0

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				double dx = x-_simVars.benchmark.coord_x*_simVars.sim.plane_domain_size[0];
				double dy = y-_simVars.benchmark.coord_y*_simVars.sim.plane_domain_size[1];

				dx /= radius;
				dy /= radius;

				double value = -50.0*2.0*dy*std::exp(-50.0*(dx*dx + dy*dy));
				value /= object_scale;

				io_data = value;

#elif ADV_FUNCTION==1

				double x = (((double)i+0.5)/(double)_simVars.disc.space_res_physical[0])*_simVars.sim.plane_domain_size[0];
				double y = (((double)j+0.5)/(double)_simVars.disc.space_res_physical[1])*_simVars.sim.plane_domain_size[1];

				x += adv_x;
				y += adv_y;

				if (x < 0)	x = _simVars.sim.plane_domain_size[0]-std::fmod(-x, _simVars.sim.plane_domain_size[0]);
				else		x = std::fmod(x, _simVars.sim.plane_domain_size[0]);

				if (y < 0)	y = _simVars.sim.plane_domain_size[1]-std::fmod(-y, _simVars.sim.plane_domain_size[1]);
				else		y = std::fmod(y, _simVars.sim.plane_domain_size[1]);

				x /= _simVars.sim.plane_domain_size[0];
				y /= _simVars.sim.plane_domain_size[1];

				io_data = freq_y*M_PI*std::sin(freq_x*M_PI*x)*std::cos(freq_y*M_PI*y)/_simVars.sim.plane_domain_size[1];
#endif
		});

		ret_h.loadPlaneDataPhysical(ret_h_phys);
		return ret_h;
	}



	void reset()
	{
		_simVars.timecontrol.current_timestep_nr = 0;
		_simVars.timecontrol.current_simulation_time = 0;

		PlaneData_Physical prog_u_phys(planeDataConfig);
		PlaneData_Physical prog_v_phys(planeDataConfig);

		prog_u_phys.physical_set_all_value(_shackProgramSpecific.velocity_u);
		prog_v_phys.physical_set_all_value(_shackProgramSpecific.velocity_v);

		prog_u.loadPlaneDataPhysical(prog_u_phys);
		prog_v.loadPlaneDataPhysical(prog_v_phys);

		prog_h = get_advected_solution(0);

	}



	void p_run_euler_timestep_update(
			const PlaneData_Spectral &i_h,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_t,		///< time updates
			PlaneData_Spectral &o_u_t,		///< time updates
			PlaneData_Spectral &o_v_t,		///< time updates

			double i_simulation_timestamp = -1
	)
	{
		double cell_size_x = _simVars.sim.plane_domain_size[0]/(double)_simVars.disc.space_res_physical[0];
		double cell_size_y = _simVars.sim.plane_domain_size[1]/(double)_simVars.disc.space_res_physical[1];


		if (_shackProgramSpecific.advection_scheme == 0)
		{
			PlaneData_Physical o_h_t_phys = o_h_t.toPhys();
			PlaneData_Physical i_h_phys = i_h.toPhys();
			PlaneData_Physical i_u_phys = i_u.toPhys();
			PlaneData_Physical i_v_phys = i_v.toPhys();

			o_h_t_phys =
				(
					(
						// u is positive
						op.shift_right(i_h_phys)*i_u_phys.physical_query_return_value_if_positive()	// inflow
						-i_h_phys*op.shift_left(i_u_phys.physical_query_return_value_if_positive())	// outflow

						// u is negative
						+(i_h_phys*i_u_phys.physical_query_return_value_if_negative())					// outflow
						-op.shift_left(i_h_phys*i_u_phys.physical_query_return_value_if_negative())	// inflow
					)*(1.0/cell_size_x)				// here we see a finite-difference-like formulation
					+
					(
						// v is positive
						op.shift_up(i_h_phys)*i_v_phys.physical_query_return_value_if_positive()		// inflow
						-i_h_phys*op.shift_down(i_v_phys.physical_query_return_value_if_positive())	// outflow

						// v is negative
						+(i_h_phys*i_v_phys.physical_query_return_value_if_negative())					// outflow
						-op.shift_down(i_h_phys*i_v_phys.physical_query_return_value_if_negative())	// inflow
					)*(1.0/cell_size_y)
				);

			o_h_t.loadPlaneDataPhysical(o_h_t_phys);
		}
		else if (_shackProgramSpecific.advection_scheme == 1)
		{
			// STAGGERED

			//             |                       |                       |
			// --v---------|-----------v-----------|-----------v-----------|
			//   h-1       u0          h0          u1          h1          u2

			PlaneData_Spectral avg_b_x_spec(o_h_t.planeDataConfig);
			PlaneData_Spectral avg_b_y_spec(o_h_t.planeDataConfig);
			avg_b_x_spec.loadPlaneDataPhysical(op.avg_b_x(i_h.toPhys()));
			avg_b_y_spec.loadPlaneDataPhysical(op.avg_b_y(i_h.toPhys()));
			// staggered
			o_h_t = -(
					op.diff_f_x(avg_b_x_spec*i_u) +
					op.diff_f_y(avg_b_y_spec*i_v)
				);
		}
		else  if (_shackProgramSpecific.advection_scheme == 2)
		{
			// NON-STAGGERED

			if (_shackProgramSpecific.staggered_use_analytical_solution == 0)
			{
				// non-staggered
				o_h_t = -(
						op.diff_c_x(i_h*i_u) +
						op.diff_c_y(i_h*i_v)
					);
			}
			else if (_shackProgramSpecific.staggered_use_analytical_solution == 1)
			{
				// non-staggered with analytical solution, only works for constant velocity!
				o_h_t = -(
						get_advected_solution_diffx(i_simulation_timestamp)*i_u +
						get_advected_solution_diffy(i_simulation_timestamp)*i_v
					);
			}
			else
			{
				std::cerr << "Usage of analytical solution not specified, use -d option [0: compute diffs on discrete solution, 1: use analytical diffs]" << std::endl;
				exit(-1);
			}
		}
		else  if (_shackProgramSpecific.advection_scheme == 3)
		{
			// NO H UPDATE
			o_h_t.spectral_set_zero();
		}
		else
		{
			std::cerr << "Advection type not specified, use -c option [0: up/downwinding, 1: staggered, 2: non-staggered]" << std::endl;
			exit(-1);
		}

		o_u_t.spectral_set_zero();
		o_v_t.spectral_set_zero();

		_simVars.timecontrol.current_timestep_nr++;
	}



	void run()
	{
	}



	void run_timestep()
	{
		timestepping.run_timestep(
				this,
				&SimulationAdvection::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
				_simVars.timecontrol.current_timestep_size,
				_simVars.disc.timestepping_order,
				_simVars.timecontrol.current_simulation_time
			);

		// provide information to parameters
		_simVars.timecontrol.current_simulation_time += _simVars.timecontrol.current_timestep_size;
		_simVars.timecontrol.current_timestep_nr++;
	}



	bool should_quit()
	{
		return false;
	}



	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (_simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	void vis_get_vis_data_array(
			const PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		int vis_id = _simVars.misc.vis_id % 6;

		switch (vis_id)
		{
		default:
		{
			PlaneData_Physical prog_h_phys = prog_h.toPhys();
			*o_dataArray = &prog_h_phys;
			break;
		}

		case 1:
		{
			tmp = get_advected_solution(_simVars.timecontrol.current_simulation_time);
			PlaneData_Physical tmp_phys = tmp.toPhys();
			*o_dataArray = &tmp_phys;
			break;
		}

		case 2:
		{
			tmp = op.diff_c_x(get_advected_solution(_simVars.timecontrol.current_simulation_time));
			PlaneData_Physical tmp_phys = tmp.toPhys();
			*o_dataArray = &tmp_phys;
			break;
		}

		case 3:
		{
			tmp = get_advected_solution_diffx(_simVars.timecontrol.current_simulation_time);
			PlaneData_Physical tmp_phys = tmp.toPhys();
			*o_dataArray = &tmp_phys;
			break;
		}

		case 4:
		{
			tmp = op.diff_c_y(get_advected_solution(_simVars.timecontrol.current_simulation_time));
			PlaneData_Physical tmp_phys = tmp.toPhys();
			*o_dataArray = &tmp_phys;
			break;
		}

		case 5:
		{
			tmp = get_advected_solution_diffy(_simVars.timecontrol.current_simulation_time);
			PlaneData_Physical tmp_phys = tmp.toPhys();
			*o_dataArray = &tmp_phys;
			break;
		}
		}

		*o_aspect_ratio = _simVars.sim.plane_domain_size[1] / _simVars.sim.plane_domain_size[0];
	}


	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				_simVars.timecontrol.current_simulation_time,
				_simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				_simVars.timecontrol.current_timestep_nr,
				_simVars.timecontrol.current_timestep_size,
				_simVars.diag.total_mass,
				_simVars.diag.total_energy,
				_simVars.diag.total_potential_enstrophy
			);

		return title_string;
	}


	void vis_pause()
	{
		_simVars.timecontrol.run_simulation_timesteps = !_simVars.timecontrol.run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			_simVars.misc.vis_id++;
			break;

		case 'V':
			_simVars.misc.vis_id--;
			break;
		}
	}


	bool instability_detected()
	{
		return !(	prog_h.toPhys().physical_reduce_boolean_all_finite() &&
					prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					prog_v.toPhys().physical_reduce_boolean_all_finite()
				);
	}
};


double compute_current_error(
		SimulationAdvection *simulationAdvection
)
{
	PlaneData_Spectral benchmark_h = simulationAdvection->get_advected_solution(_simVars.timecontrol.current_simulation_time);


	return (simulationAdvection->prog_h-benchmark_h).toPhys().physical_reduce_rms();
}


int main(
		int i_argc,
		char *i_argv[]
)
{
	_simVars.setupFromMainParameters(i_argc, i_argv);

	{
		// Quick and dirty solution: Just get the basic configuration
		sweet::ProgramArguments pa(false, false);
		pa.setup(i_argc, i_argv);

		_shackProgramSpecific.processProgramArguments(pa);
	}

	_simVars.outputConfig();
	_shackProgramSpecific.printClass();

	double u = _shackProgramSpecific.velocity_u;
	double v = _shackProgramSpecific.velocity_v;

	double total_speed;
	double turnaround_time;
	if (u == 0 && v == 0)
	{
		std::cerr << "Both velocity components are zero, EXIT" << std::endl;
		exit(1);
	}

	if (u != 0 && v == 0)
	{
		total_speed = u;
		turnaround_time = _simVars.sim.plane_domain_size[0]/u;
	}
	else if (u == 0 && v != 0)
	{
		total_speed = v;
		turnaround_time = _simVars.sim.plane_domain_size[1]/v;
	}
	else
	{
		total_speed = v;
		if (std::abs(_simVars.sim.plane_domain_size[1]/_simVars.sim.plane_domain_size[0]-v/u) > 0.000000001)
		{
			std::cerr << "ratio of domain sizes and speed have to be similar" << std::endl;
			exit(1);
		}

		total_speed = std::sqrt(u*u+v*v);
		double diagonal = std::sqrt(_simVars.sim.plane_domain_size[0]*_simVars.sim.plane_domain_size[0] + _simVars.sim.plane_domain_size[1]*_simVars.sim.plane_domain_size[1]);
		turnaround_time = diagonal/total_speed;
	}

	if (_simVars.misc.verbosity > 1)
	{
		std::cout << "Turnaround time: " << turnaround_time << std::endl;
		std::cout << "Total speed: " << total_speed << std::endl;
	}

#if SWEET_GUI
	if (_simVars.misc.gui_enabled)
	{
		SimulationAdvection *simulationAdvection = new SimulationAdvection;
		VisSweet<SimulationAdvection> visSweet(simulationAdvection);
		delete simulationAdvection;
		return 0;
	}
#endif

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	// allocate data storage for computed errors

	bool error_detected = false;

	if (_shackProgramSpecific.test_mode == 0)
	{
		std::ostringstream output_string_conv;

		double *computed_errors = new double[1024];
		double *conv_rate = new double[1024];

		std::size_t res_x = _simVars.disc.space_res_physical[0];
		std::size_t res_y = _simVars.disc.space_res_physical[1];

		std::size_t max_res = 128;

		if (res_x > max_res || res_y > max_res)
			max_res = std::max(res_x, res_y);

		for (	int res_iterator_id = 0;
				res_x <= max_res && res_y <= max_res;
				res_x *= 2, res_y *= 2, res_iterator_id++
		)
		{
			output_string_conv << std::endl;
			output_string_conv << res_x << "x" << res_y << "\t";

			std::cout << "*******************************************************************************" << std::endl;
			std::cout << "Testing convergence with resolution " << res_x << " x " << res_y << " and RK order " << _simVars.disc.timestepping_order << std::endl;
			std::cout << "*******************************************************************************" << std::endl;

			_simVars.disc.space_res_physical[0] = res_x;
			_simVars.disc.space_res_physical[1] = res_y;
			_simVars.disc.space_res_spectral[0] = 0;
			_simVars.disc.space_res_spectral[1] = 0;

			_simVars.timecontrol.current_simulation_time = 0;
			_simVars.timecontrol.current_timestep_nr = 0;
			_simVars.timecontrol.current_timestep_size *= 0.5;
			_simVars.timecontrol.setup_timestep_size  = _simVars.timecontrol.current_timestep_size;

			std::cout << " + current_timestep_size: " << _simVars.timecontrol.current_timestep_size << std::endl;

			_simVars.reset();

			planeDataConfigInstance.setupAutoSpectralSpace(_simVars.disc.space_res_physical, _simVars.misc.reuse_spectral_transformation_plans);

			SimulationAdvection *simulationAdvection = new SimulationAdvection(planeDataConfig);

			Stopwatch time;
			time.reset();

			while(true)
			{
				if (_simVars.misc.verbosity >= 10)
					std::cout << "time: " << _simVars.timecontrol.current_simulation_time << std::endl;

				simulationAdvection->run_timestep();

				if (simulationAdvection->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}

				bool print_output = false;
				if (turnaround_time <= _simVars.timecontrol.current_simulation_time)
					print_output = true;

				if (_simVars.timecontrol.max_simulation_time != -1)
					if (_simVars.timecontrol.current_simulation_time >= _simVars.timecontrol.max_simulation_time-_simVars.timecontrol.current_simulation_time*1e-10)
						print_output = true;

				if (print_output)
				{
					double &this_error = computed_errors[res_iterator_id];
					//double &this_conv_rate_space = conv_rate[res_iterator_id];

					/*
					 * TODO: Something's not right here.
					 * If we don't reduce the time step size, the errors are not converging.
					 */
					double error = compute_current_error(simulationAdvection);
					std::cout << "RMS error in height: " << error << std::endl;

//					double error_max = (simulationAdvection->prog_h-benchmark_h).reduce_maxAbs();
//					std::cout << "Max error in height: " << error_max << std::endl;

					this_error = error;

					double eps = 0.1;
					/*
					 * check convergence in space
					 */
					if (res_iterator_id > 0)
					{
						double &prev_error_space = computed_errors[(res_iterator_id-1)];

						double expected_conv_rate = std::pow(2.0, (double)(_simVars.disc.timestepping_order));
						double this_conv_rate_space = prev_error_space / this_error;

						std::cout << "          Norm2 convergence rate (space): " << this_conv_rate_space << ", expected: " << expected_conv_rate << std::endl;

						if (std::abs(this_conv_rate_space-expected_conv_rate) > eps*expected_conv_rate)
						{
							if (error < 10e-12)
							{
								std::cerr << "Warning: Ignoring this error, since it's below machine precision" << std::endl;
							}
							else
							{
								std::cerr << "Convergence rate threshold (" << eps*expected_conv_rate << ") exceeded" << std::endl;
								error_detected = true;
							}
						}

						output_string_conv << this_conv_rate_space << "\t";
					}
					break;
				}
			}	// while true

			time.stop();

			double seconds = time();

			std::cout << "Simulation time: " << seconds << " seconds" << std::endl;
			std::cout << "Time per time step: " << seconds/(double)_simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;

			delete simulationAdvection;

		}	// res

		delete [] computed_errors;
		delete [] conv_rate;

		std::cout << std::endl;
		std::cout << "Convergence rate in space (inc. resolution):";
		std::cout << output_string_conv.str() << std::endl;
	}
#if 1
	else if (_shackProgramSpecific.test_mode == 1)
	{
		std::ostringstream output_string_conv;

		double *computed_errors = new double[1024];
		double *conv_rate = new double[1024];

		planeDataConfigInstance.setupAutoSpectralSpace(_simVars.disc.space_res_physical, _simVars.misc.reuse_spectral_transformation_plans);


		for (	int cfl_iterator_id = 0;
				cfl_iterator_id < 7;
				_simVars.timecontrol.current_timestep_size *= 0.5, cfl_iterator_id++
		)
		{
			output_string_conv << std::endl;

			std::cout << "*********************************************************************************************************" << std::endl;
			std::cout << "Testing time convergence with time step size " << _simVars.timecontrol.current_timestep_size << " and RK order " << _simVars.disc.timestepping_order << std::endl;
			std::cout << "*********************************************************************************************************" << std::endl;

			SimulationAdvection simulationAdvection(planeDataConfig);
			simulationAdvection.reset();

			Stopwatch time(true);

			while(true)
			{
				if (_simVars.misc.verbosity >= 10)
					std::cout << "time: " << _simVars.timecontrol.current_simulation_time << std::endl;

				simulationAdvection.run_timestep();

				if (simulationAdvection.instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}

				bool print_output = false;
				if (turnaround_time <= _simVars.timecontrol.current_simulation_time)
					print_output = true;

				if (_simVars.timecontrol.max_simulation_time != -1)
					if (_simVars.timecontrol.current_simulation_time >= _simVars.timecontrol.max_simulation_time)
						print_output = true;

				if (print_output)
				{
					double &this_error = computed_errors[cfl_iterator_id];

					double error = compute_current_error(&simulationAdvection);
					std::cout << "Error in height: " << error << std::endl;

//					double error_max = (simulationAdvection->prog_h-benchmark_h).reduce_maxAbs();
//					std::cout << "Max error in height: " << error_max << std::endl;

					double cell_size_x = _simVars.sim.plane_domain_size[0]/(double)_simVars.disc.space_res_physical[0];
					double cell_size_y = _simVars.sim.plane_domain_size[1]/(double)_simVars.disc.space_res_physical[1];

					std::cout << "          dt = " << _simVars.timecontrol.current_timestep_size << "    dx = " << cell_size_x << " x " << cell_size_y << std::endl;

					this_error = error;

					double eps = 0.1;

					/*
					 * check convergence in time
					 */
					if (cfl_iterator_id > 0)
					{
						double &prev_error_space = computed_errors[(cfl_iterator_id-1)];

						double expected_conv_rate = std::pow(2.0, (double)(_simVars.disc.timestepping_order));
						double this_conv_rate_space = prev_error_space / this_error;

						std::cout << "          Norm2 convergence rate (time): " << this_conv_rate_space << ", expected: " << expected_conv_rate << std::endl;

						if (std::abs(this_conv_rate_space-expected_conv_rate) > eps*expected_conv_rate)
						{
							if (error < 10e-12)
							{
								std::cerr << "Warning: Ignoring this error, since it's below machine precision" << std::endl;
							}
							else
							{
								std::cerr << "Convergence rate threshold (" << eps*expected_conv_rate << ") exceeded" << std::endl;
								error_detected = true;
							}
						}

						output_string_conv << "r=" << this_conv_rate_space << "\t";
						output_string_conv << "dt=" << _simVars.timecontrol.current_timestep_size << "\t";
						output_string_conv << "dx=" << cell_size_x << "." << cell_size_x;
					}
					break;

				}

			}	// while true

			time.stop();

			double seconds = time();

			std::cout << "Simulation time: " << seconds << " seconds" << std::endl;
			std::cout << "Time per time step: " << seconds/(double)_simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;

		}	// res

		delete [] computed_errors;
		delete [] conv_rate;

		std::cout << std::endl;
		std::cout << "Convergence rate in time (inc. resolution):";
		std::cout << output_string_conv.str() << std::endl;
	}
#endif
	else
	{
		SWEETError("Not supported!");
	}

	if (error_detected)
	{
		std::cerr << "There was an error in the convergence tests" << std::endl;
		exit(1);
	}

	return 0;
}
