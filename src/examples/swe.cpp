
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include "sweet/Operators2D.hpp"

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>

SimulationParameters parameters;


class SimulationSWE
{
public:
	DataArray<2> prog_h, prog_u, prog_v;
	DataArray<2> eta;
	DataArray<2> tmp;

	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

public:
	SimulationSWE(
	)	:
		prog_h(parameters.res),
		prog_u(parameters.res),
		prog_v(parameters.res),

		eta(parameters.res),
		tmp(parameters.res),

		op(parameters.res, parameters.sim_domain_length, parameters.use_spectral_diffs)
	{
		reset();
	}


	void reset()
	{
		last_timestep_nr_update_diagnostics = -1;

		parameters.status_timestep_nr = 0;
		parameters.status_simulation_time = 0;


		prog_h.setAll(parameters.setup_h0);
		prog_u.setAll(0);
		prog_v.setAll(0);

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
				double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

				prog_h.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
				prog_u.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
				prog_v.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
			}
		}
	}


	void update_diagnostics()
	{

		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == parameters.status_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = parameters.status_timestep_nr;

		double normalization = (parameters.sim_domain_length[0]*parameters.sim_domain_length[1]) /
								((double)parameters.res[0]*(double)parameters.res[1]);

		// mass
		parameters.diagnostics_mass = prog_h.reduce_sum_quad() * normalization;

		// energy
		parameters.diagnostics_energy = 0.5*(
				prog_h*prog_h +
				prog_h*prog_u*prog_u +
				prog_h*prog_v*prog_v
			).reduce_sum_quad() * normalization;

		// potential enstropy
		eta = (op.diff_c_x(prog_v) - op.diff_c_y(prog_u) + parameters.sim_f) / prog_h;
		parameters.diagnostics_potential_entrophy = 0.5*(eta*eta*prog_h).reduce_sum_quad() * normalization;
	}


	void compute_upwinding_P_updates(
			const DataArray<2> &i_P,		///< prognostic variables (at T=tn)
			const DataArray<2> &i_u,		///< prognostic variables (at T=tn+dt)
			const DataArray<2> &i_v,		///< prognostic variables (at T=tn+dt)

			DataArray<2> &o_P_t				///< time updates (at T=tn+dt)
	)
	{
		std::cerr << "TODO: implement, is this really possible for non-staggered grid? (averaging of velocities required)" << std::endl;
		exit(-1);
		//             |                       |                       |
		// --v---------|-----------v-----------|-----------v-----------|
		//   h-1       u0          h0          u1          h1          u2
		//

		// same a above, but formulated in a finite-difference style
		o_P_t =
			(
				(
					// u is positive
					op.shift_right(i_P)*i_u.return_value_if_positive()	// inflow
					-i_P*op.shift_left(i_u.return_value_if_positive())					// outflow

					// u is negative
					+(i_P*i_u.return_value_if_negative())	// outflow
					-op.shift_left(i_P*i_u.return_value_if_negative())		// inflow
				)*(1.0/parameters.sim_cell_size[0])	// here we see a finite-difference-like formulation
				+
				(
					// v is positive
					op.shift_up(i_P)*i_v.return_value_if_positive()		// inflow
					-i_P*op.shift_down(i_v.return_value_if_positive())					// outflow

					// v is negative
					+(i_P*i_v.return_value_if_negative())	// outflow
					-op.shift_down(i_P*i_v.return_value_if_negative())	// inflow
				)*(1.0/parameters.sim_cell_size[1])
			);
	}





	void p_run_euler_timestep_update(
			const DataArray<2> &i_h,	///< prognostic variables
			const DataArray<2> &i_u,	///< prognostic variables
			const DataArray<2> &i_v,	///< prognostic variables

			DataArray<2> &o_h_t,	///< time updates
			DataArray<2> &o_u_t,	///< time updates
			DataArray<2> &o_v_t,	///< time updates

			double &o_dt,			///< time step restriction
			double i_fixed_dt = 0		///< if this value is not equal to 0, use this time step size instead of computing one
	)
	{
		/*
		 * non-conservative formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */
		o_u_t = -parameters.sim_g*op.diff_c_x(i_h) - i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u) + parameters.sim_f*i_v;
		o_v_t = -parameters.sim_g*op.diff_c_y(i_h) - i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v) - parameters.sim_f*i_u;

		if (parameters.sim_viscocity != 0)
		{
			o_v_t += (op.diff2_c_y(i_u) + op.diff2_c_y(i_v))*parameters.sim_viscocity;
			o_u_t += (op.diff2_c_x(i_u) + op.diff2_c_x(i_v))*parameters.sim_viscocity;
		}

		if (parameters.sim_hyper_viscocity != 0)
		{
			o_u_t += (op.diff2_c_x(op.diff2_c_x(i_u)) + op.diff2_c_x(op.diff2_c_x(i_v)))*parameters.sim_hyper_viscocity;
			o_v_t += (op.diff2_c_y(op.diff2_c_y(i_u)) + op.diff2_c_y(op.diff2_c_y(i_v)))*parameters.sim_hyper_viscocity;
		}


		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt != 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			/*
			 * If the timestep size parameter is negative, we use the absolute value of this one as the time step size
			 */
			if (i_fixed_dt < 0)
			{
				o_dt = -i_fixed_dt;
			}
			else
			{
				double limit_speed = std::max(parameters.sim_cell_size[0]/i_u.reduce_maxAbs(), parameters.sim_cell_size[1]/i_v.reduce_maxAbs());

				// limit by re
				double limit_visc = limit_speed;
		//        if (viscocity > 0)
		//           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

				// limit by gravitational acceleration
				double limit_gh = std::min(parameters.sim_cell_size[0], parameters.sim_cell_size[1])/std::sqrt(parameters.sim_g*i_h.reduce_maxAbs());

		//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
				o_dt = parameters.sim_CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
			}
		}

		// TODO: FIX
		o_h_t = -op.diff_c_y(i_v*i_h);
		o_h_t -= op.diff_c_x(i_u*i_h);
		return;

		if (!parameters.timestepping_leapfrog_like_update)
		{
			if (!parameters.timestepping_up_and_downwinding)
			{
				// standard update
				o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
				return;
			}
			else
			{
				// up/down winding
				compute_upwinding_P_updates(
						i_h,
						i_u,
						i_v,
						o_h_t
					);
				return;
			}
		}
		else
		{
			/*
			 * a kind of leapfrog:
			 *
			 * We use the hew v and u values to compute the update for p
			 *
			 * compute updated u and v values without using it
			 */
			if (!parameters.timestepping_up_and_downwinding)
			{
				// recompute U and V

				// update based on new u and v values
				o_h_t = -op.diff_c_x(
							i_h*(i_u+o_dt*o_u_t)
						)
						- op.diff_c_y(
							i_h*(i_v+o_dt*o_v_t)
						);
				return;
			}
			else
			{
				// update based on new u and v values
				compute_upwinding_P_updates(
						i_h,
						i_u+o_dt*o_u_t,
						i_v+o_dt*o_v_t,
						o_h_t
					);
				return;
			}
		}


	}



	void run_timestep()
	{
		double dt;
		timestepping.run_rk_timestep(
				this,
				&SimulationSWE::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
				prog_h, prog_u, prog_v,
				dt,
				parameters.timestepping_timestep_size,
				parameters.timestepping_runge_kutta_order
			);

		// provide information to parameters
		parameters.status_simulation_timestep_size = dt;
		parameters.status_simulation_time += dt;
		parameters.status_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}



public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (parameters.verbosity > 0)
		{
			update_diagnostics();

			if (parameters.status_timestep_nr == 0)
			{
				o_ostream << "T\tMASS\tENERGY\tPOT_ENSTROPHY";

				if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
					o_ostream << "\tABS_P_DT\tABS_U_DT\tABS_V_DT";

				o_ostream << std::endl;
			}

			o_ostream << parameters.status_simulation_time << "\t" << parameters.diagnostics_mass << "\t" << parameters.diagnostics_energy << "\t" << parameters.diagnostics_potential_entrophy;

			// this should be zero for the steady state test
			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				double test_val;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
					}

				test_val = (prog_h-tmp).reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << test_val;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// u space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
					}

				test_val = (prog_u-tmp).reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << test_val;

				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
					}

				test_val = (prog_v-tmp).reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << test_val;
			}

			o_ostream << std::endl;
		}
	}


public:
	bool should_quit()
	{
		if (parameters.max_timesteps_nr != -1 && parameters.max_timesteps_nr <= parameters.status_timestep_nr)
			return true;

		if (parameters.max_simulation_time != -1 && parameters.max_simulation_time <= parameters.status_simulation_time)
			return true;

		return false;
	}


	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}


	struct VisStuff
	{
		const DataArray<2>* data;
		const char *description;
	};

	VisStuff vis_arrays[4] =
	{
			{&prog_h,	"h"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
			{&eta,		"eta"}
	};

	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = parameters.sim_domain_length[1] / parameters.sim_domain_length[0];

#if 0

		if (parameters.vis_id == 0)
		{
			eta = op.diff_c_x;
//			eta = op.diff_c_y(prog_h);
			*o_dataArray = &eta;
		}
		else if (parameters.vis_id == 1)
		{
			eta = op.diff_c_x(prog_h)-op.diff_c_y(prog_h);
			*o_dataArray = &eta;
		}
		else if (parameters.vis_id == 2)
		{
			eta = op.diff_c_x(prog_h);
			*o_dataArray = &eta;
		}
		else if (parameters.vis_id == 3)
		{
			eta = op.diff_c_y(prog_h);
			*o_dataArray = &eta;
		}

#endif
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
		update_diagnostics();

		int id = parameters.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[2048];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				parameters.status_simulation_time,
				parameters.status_simulation_time/(60.0*60.0*24.0),
				parameters.status_timestep_nr,
				parameters.status_simulation_timestep_size,
				vis_arrays[id].description,
				parameters.diagnostics_mass, parameters.diagnostics_energy, parameters.diagnostics_potential_entrophy);

		return title_string;
	}


	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			parameters.vis_id++;
			break;

		case 'V':
			parameters.vis_id--;
			break;
		}
	}


	bool instability_detected()
	{
		return !(	prog_h.reduce_all_finite() &&
					prog_u.reduce_all_finite() &&
					prog_v.reduce_all_finite()
				);
	}
};



int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	parameters.setup(i_argc, i_argv);

	SimulationSWE *simulationSWE = new SimulationSWE;

	std::ostringstream buf;
	buf << std::setprecision(14);

#if SWEET_GUI
	VisSweet<SimulationSWE> visSweet(simulationSWE);
#else
	simulationSWE->reset();

	while(true)
	{
		if (parameters.verbosity > 1)
		{
			simulationSWE->timestep_output(buf);

			std::string output = buf.str();
			buf.str("");

			std::cout << output;

			if (parameters.verbosity > 2)
				std::cerr << output;
		}

		if (simulationSWE->should_quit())
			break;

		simulationSWE->run_timestep();

		if (simulationSWE->instability_detected())
		{
			std::cout << "INSTABILITY DETECTED" << std::endl;
			break;
		}
	}
#endif

	delete simulationSWE;

	return 1;
}
