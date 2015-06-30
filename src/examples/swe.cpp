
#include <sweet/DataArray.hpp>
#include "VisSweet.hpp"
#include "Parameters.hpp"
#include "sweet/Operators2D.hpp"

#include <unistd.h>

Parameters parameters;


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

		op(parameters.cell_size, parameters.res)
	{

		if (0)
		{
			double diff_c_x_kernel[5][5] = {
					{0,0,0,0,0},{0,0,0,0,0},
					{1.0, -8.0, 0, 8.0, -1.0},
					{0,0,0,0,0},{0,0,0,0,0},
			};
			op.diff_c_x.setup_kernel(diff_c_x_kernel, 1.0/(12*parameters.cell_size[0]));

			double diff_c_y_kernel[5][5] = {
					{0,0, 1.0, 0,0},
					{0,0, -8.0, 0,0},
					{0,0,  0.0, 0,0},
					{0,0,  8.0, 0,0},
					{0,0, -1.0, 0,0},
			};
			op.diff_c_y.setup_kernel(diff_c_y_kernel, 1.0/(12*parameters.cell_size[1]));

			double diff2_c_x_kernel[5][5] = {
					{0,0,0,0,0},{0,0,0,0,0},
					{-1.0, 16.0, -30.0, 16.0, -1.0},
					{0,0,0,0,0},{0,0,0,0,0},
			};
			op.diff2_c_x.setup_kernel(diff2_c_x_kernel, 1.0/(12*parameters.cell_size[0]));

			double diff2_c_y_kernel[5][5] = {
					{0, 0,  -1.0, 0, 0},
					{0, 0,  16.0, 0, 0},
					{0, 0, -30.0, 0, 0},
					{0, 0,  16.0, 0, 0},
					{0, 0,  -1.0, 0, 0},
			};
			op.diff2_c_y.setup_kernel(diff2_c_y_kernel, 1.0/12*parameters.cell_size[1]);
		}

		reset();
	}


	void reset()
	{
		parameters.timestep_nr = 0;
		parameters.simulation_time = 0;

		h.data_setall(parameters.h0);

		double center_x = 0.7;
		double center_y = 0.6;

		if (parameters.setup_scenario == 0)
		{
			/*
			 * radial dam break
			 */
			double radius = 0.2;
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					if (radius*radius >= dx*dx+dy*dy)
						h.getDataRef(j,i) += 1.0;
				}
			}
		}

		if (parameters.setup_scenario == 1)
		{
			/*
			 * fun with Gaussian
			 */
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)/(double)parameters.res[0];
					double y = ((double)j+0.5)/(double)parameters.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					h.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
				}
			}
		}

		u.data_setall(0);
		v.data_setall(0);
	}


	void run_timestep()
	{
		/*
		 * non-conservative formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y
		 *	v_t = -g * h_y - u * v_x - v * v_y
		 */
		u_t = -parameters.g*op.diff_c_x(h) - u*op.diff_c_x(u) - v*op.diff_c_y(u);
		v_t = -parameters.g*op.diff_c_y(h) - u*op.diff_c_x(v) - v*op.diff_c_y(v);

		// mass
		parameters.mass = h.reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// energy
		parameters.energy = 0.5*(
				h*h +
				h*u*u +
				h*v*v
			).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// potential enstropy
		// TODO: is this correct?
		eta = (op.diff_c_x(v) - op.diff_c_y(u)) / h;
		parameters.potential_entrophy = 0.5*(eta*eta*h).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);

		if (parameters.viscocity != 0)
		{
			// TODO: is this correct?
			v_t += (op.diff2_c_y(u) + op.diff2_c_y(v))*parameters.viscocity;
			u_t += (op.diff2_c_x(u) + op.diff2_c_x(v))*parameters.viscocity;
		}

		if (parameters.hyper_viscocity != 0)
		{
			// TODO: is this correct?
			u_t += (op.diff2_c_x(op.diff2_c_x(u)) + op.diff2_c_x(op.diff2_c_x(v)))*parameters.hyper_viscocity;
			v_t += (op.diff2_c_y(op.diff2_c_y(u)) + op.diff2_c_y(op.diff2_c_y(v)))*parameters.hyper_viscocity;
		}

		double limit_speed = std::max(parameters.cell_size[0]/u.reduce_maxAbs(), parameters.cell_size[1]/v.reduce_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
 //           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
		double limit_gh = std::min(parameters.cell_size[0], parameters.cell_size[1])/std::sqrt(parameters.g*h.reduce_maxAbs());

//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
		double dt = parameters.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

		// provide information to parameters
		parameters.timestep_size = dt;


		u += dt*u_t;
		v += dt*v_t;

		h_t = -op.diff_c_x(u*h) - op.diff_c_y(v*h);
		h += dt*h_t;

		parameters.simulation_time += dt;
		parameters.timestep_nr++;
	}



	void vis_post_frame_processing(int i_num_iterations)
	{
		if (parameters.run_simulation)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}

	DataArray<2> &vis_get_vis_data_array()
	{
		return h;
	}

	const char* vis_get_status_string()
	{
		static char title_string[1024];
		sprintf(title_string, "Time: %f, Timestep: %i, timestep size: %e, Mass: %e, Energy: %e, Potential Entrophy: %e", parameters.simulation_time, parameters.timestep_nr, parameters.timestep_size, parameters.mass, parameters.energy, parameters.potential_entrophy);
		return title_string;
	}

	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
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
