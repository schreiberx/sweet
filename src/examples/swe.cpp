
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include "sweet/Operators2D.hpp"

#include <unistd.h>

SimulationParameters parameters;


class SimulationSWE
{
public:
	DataArray<2> h, u, v;
	DataArray<2> h_t, u_t, v_t;
	DataArray<2> eta;

	Operators2D op;

public:
	SimulationSWE(
	)	:
		h(parameters.res),
		u(parameters.res),
		v(parameters.res),
		h_t(parameters.res),
		u_t(parameters.res),
		v_t(parameters.res),

		eta(parameters.res),

		op(parameters.sim_cell_size, parameters.res)
	{
		// override gravitation
		parameters.sim_g = 1.0;

		reset();
	}


	void reset()
	{
		parameters.status_timestep_nr = 0;
		parameters.status_simulation_time = 0;


		h.data_setall(parameters.setup_h0);
		u.data_setall(0);
		v.data_setall(0);

		for (std::size_t j = 0; j < parameters.res[1]; j++)
		{
			for (std::size_t i = 0; i < parameters.res[0]; i++)
			{
				double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
				double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

				h.getDataRef(j,i) = SWEValidationBenchmarks::return_h(parameters, x, y);
				u.getDataRef(j,i) = SWEValidationBenchmarks::return_u(parameters, x, y);
				v.getDataRef(j,i) = SWEValidationBenchmarks::return_v(parameters, x, y);
			}
		}
	}


	bool run_timestep()
	{
		/*
		 * non-conservative formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y
		 *	v_t = -g * h_y - u * v_x - v * v_y
		 */
		u_t = -parameters.sim_g*op.diff_c_x(h) - u*op.diff_c_x(u) - v*op.diff_c_y(u);
		v_t = -parameters.sim_g*op.diff_c_y(h) - u*op.diff_c_x(v) - v*op.diff_c_y(v);

		// mass
		parameters.diagnostics_mass = h.reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// energy
		parameters.diagnostics_energy = 0.5*(
				h*h +
				h*u*u +
				h*v*v
			).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// potential enstropy
		// TODO: is this correct?
		eta = (op.diff_c_x(v) - op.diff_c_y(u)) / h;
		parameters.diagnostics_potential_entrophy = 0.5*(eta*eta*h).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);

		if (parameters.sim_viscocity != 0)
		{
			v_t += (op.diff2_c_y(u) + op.diff2_c_y(v))*parameters.sim_viscocity;
			u_t += (op.diff2_c_x(u) + op.diff2_c_x(v))*parameters.sim_viscocity;
		}

		if (parameters.sim_hyper_viscocity != 0)
		{
			u_t += (op.diff2_c_x(op.diff2_c_x(u)) + op.diff2_c_x(op.diff2_c_x(v)))*parameters.sim_hyper_viscocity;
			v_t += (op.diff2_c_y(op.diff2_c_y(u)) + op.diff2_c_y(op.diff2_c_y(v)))*parameters.sim_hyper_viscocity;
		}

		double limit_speed = std::max(parameters.sim_cell_size[0]/u.reduce_maxAbs(), parameters.sim_cell_size[1]/v.reduce_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
 //           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
		double limit_gh = std::min(parameters.sim_cell_size[0], parameters.sim_cell_size[1])/std::sqrt(parameters.sim_g*h.reduce_maxAbs());

//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
		double dt = parameters.sim_CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

		// provide information to parameters
		parameters.timestepping_timestep_size = dt;

		u += dt*u_t;
		v += dt*v_t;

		// update h based on updated velocities
		h_t = -op.diff_c_x(u*h) - op.diff_c_y(v*h);
		h += dt*h_t;

		parameters.status_simulation_time += dt;
		parameters.status_timestep_nr++;
	}


	bool should_quit()
	{
		return false;
	}


	/**
	 * postprocessing of frame: do time stepping
	 */
	bool vis_post_frame_processing(int i_num_iterations)
	{
		bool retval = true;
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				retval = run_timestep();

		return retval;
	}


	void vis_get_vis_data_array(
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		int id = parameters.vis_id % 7;

		switch(id)
		{
			case 0:	*o_dataArray = &h;	break;
			case 1: *o_dataArray = &u;	break;
			case 2: *o_dataArray = &v;	break;
			case 4: *o_dataArray = &eta;	break;
		}

		*o_aspect_ratio = parameters.sim_domain_length[1] / parameters.sim_domain_length[0];
	}

	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Time: %f, Timestep: %i, timestep size: %e, Mass: %e, Energy: %e, Potential Entrophy: %e", parameters.status_simulation_time, parameters.status_timestep_nr, parameters.timestepping_timestep_size, parameters.diagnostics_mass, parameters.diagnostics_energy, parameters.diagnostics_potential_entrophy);
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
};



int main(int i_argc, char *i_argv[])
{
	parameters.setup(i_argc, i_argv);

	SimulationSWE *simulationSWE = new SimulationSWE;

	VisSweet<SimulationSWE> visSweet(simulationSWE);

	delete simulationSWE;

	return 1;
}
