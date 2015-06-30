
#include <sweet/DataArray.hpp>
#include "VisSweet.hpp"
#include "Parameters.hpp"
#include "sweet/Operators2D.hpp"
#include <unistd.h>

Parameters parameters;


class SimulationSWE
{
public:
	DataArray<2> h;
	DataArray<2> u;
	DataArray<2> v;

	DataArray<2> hu;
	DataArray<2> hv;

	DataArray<2> h_t;

	Operators2D op;

public:
	SimulationSWE(
	)	:
		h(parameters.res),
		u(parameters.res),
		v(parameters.res),

		hu(parameters.res),
		hv(parameters.res),

		h_t(parameters.res),

		op(parameters.cell_size, parameters.res)
	{
		reset();
	}


	void reset()
	{
		parameters.timestep_nr = 0;

		h.data_setall(parameters.h0);

		u.data_setall(parameters.bogus_var0);
		v.data_setall(parameters.bogus_var1);

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
	}


	void run_timestep()
	{
        double dt = parameters.CFL*std::min(parameters.cell_size[0]/u.reduce_maxAbs(), parameters.cell_size[1]/v.reduce_maxAbs());
        if (std::isinf(dt))
        	dt = parameters.CFL*parameters.cell_size[0]/0.000001;

        parameters.timestep_size = dt;

        // 0: staggered
        // 1: non-staggered
        // 2: up/downwinding

#define ADVECTION_METHOD	2

#if ADVECTION_METHOD == 0

        // staggered
		h -= dt*(
				op_d_f_x(op_avg_b_x(h)*u) +
				op_d_f_y(op_avg_b_y(h)*v)
			);
#endif

#if ADVECTION_METHOD == 1
		// non-staggered
		h = h - dt*(
				op_d_c_x(h*u) +
				op_d_c_y(h*v)
			);
#endif

#if ADVECTION_METHOD == 2

//		std::cout << u << std::endl;

		// up/downwinding on staggered grid (FV-like method)
		h += dt/(parameters.cell_size[0]*parameters.cell_size[1])*(
				parameters.cell_size[1]*(
					//             |                       |                       |
					// --v---------|-----------v-----------|-----------v-----------|
					//   h0        u0          h1          u1          h2          u2
					//

					// left edge
					op.shift_right(h*u.return_value_if_positive())		// inflow
					+(h*op.shift_right(u.return_value_if_negative()))	// outflow

					// right edge
					-h*u.return_value_if_positive()					// outflow
					-op.shift_left(h)*u.return_value_if_negative()	// inflow
				)
				+
				parameters.cell_size[0]*(
						// bottom edge
						op.shift_up(h*v.return_value_if_positive())	// inflow
						+(h*op.shift_up(v.return_value_if_negative()))	// outflow

						// top edge
						-h*v.return_value_if_positive()					// outflow
						-op.shift_down(h)*v.return_value_if_negative()	// inflow
				)
			);
#endif
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
		sprintf(title_string, "Timestep: %i, timestep size: %e, Mass: %e, Energy: %e, Potential Entrophy: %e", parameters.timestep_nr, parameters.timestep_size, parameters.mass, parameters.energy, parameters.potential_entrophy);
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
