
#include <sweet/DataArray.hpp>
#include "VisSweet.hpp"
#include "Parameters.hpp"
#include "sweet/Operators2D.hpp"
#include <unistd.h>

Parameters parameters;


class SimulationSWEStaggered
{
public:
	// prognostics
	DataArray<2> P, u, v;

	// diagnostics
	DataArray<2> H, U, V;
	DataArray<2> eta;

	DataArray<2> P_t, u_t, v_t;


	Operators2D op;

	/**
	 * See "The Dynamics of Finite-Difference Models of the Shallow-Water Equations", Robert Sadourny
	 *
	 * Prognostic:
	 *     V_t + \eta N x (P V) + grad(P + 0.5 V.V) = 0
	 *     P_t + div(P V) = 0
	 *
	 * Potential vorticity:
	 *     \eta = rot (V) / P
	 *
	 *   ______u______
	 *   |           |
	 *   |           |
	 *   v     P     v
	 *   |           |
	 *   |_____u_____|
	 */
public:
	SimulationSWEStaggered(
	)	:
		P(parameters.res),	// density/pressure
		u(parameters.res),	// velocity (x-direction)
		v(parameters.res),	// velocity (y-direction)

		H(parameters.res),	//
		U(parameters.res),	// mass flux (x-direction)
		V(parameters.res),	// mass flux (y-direction)
		eta(parameters.res),
		P_t(parameters.res),
		u_t(parameters.res),
		v_t(parameters.res),

		op(parameters.cell_size, parameters.res)
	{
		reset();
	}

	void reset()
	{
		parameters.timestep_nr = 0;
		parameters.simulation_time = 0;

		P.data_setall(parameters.h0);

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
						P.getDataRef(j,i) += 1.0;
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

					P.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
				}
			}
		}

		u.data_setall(0);
		v.data_setall(0);
	}


	void run_timestep()
	{
		U = op.avg_f_x(P)*u;
		V = op.avg_f_y(P)*v;

		H = P + 0.5*(op.avg_b_x(u*u) + op.avg_b_y(v*v));

		eta = (op.diff_f_x(v) - op.diff_f_y(u)) / op.avg_f_x(op.avg_f_y(P));

		// mass
		parameters.mass = P.reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// energy
		parameters.energy = 0.5*(
				P*P +
				P*op.avg_b_x(u*u) +
				P*op.avg_b_y(v*v)
			).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);
		// potential enstropy
		parameters.potential_entrophy = 0.5*(eta*eta*op.avg_f_x(op.avg_f_y(P))).reduce_sum() / (double)(parameters.res[0]*parameters.res[1]);

		u_t = op.avg_b_y(eta*op.avg_f_x(V)) - op.diff_f_x(H);
		v_t = op.avg_b_x(eta*op.avg_f_y(U)) - op.diff_f_y(H);
		P_t = -op.diff_b_x(U) - op.diff_b_y(V);

		if (parameters.viscocity != 0)
		{
			u_t += (op.diff2_c_x(u) + op.diff2_c_x(v))*parameters.viscocity;
			v_t += (op.diff2_c_y(u) + op.diff2_c_y(v))*parameters.viscocity;
		}

		if (parameters.hyper_viscocity != 0)
		{
			u_t += (op.diff2_c_x(op.diff2_c_x(u)) + op.diff2_c_x(op.diff2_c_x(v)))*parameters.hyper_viscocity;
			v_t += (op.diff2_c_y(op.diff2_c_y(u)) + op.diff2_c_y(op.diff2_c_y(v)))*parameters.hyper_viscocity;
		}

		double limit_speed = std::max(parameters.cell_size[0]/u.reduce_maxAbs(), parameters.cell_size[1]/v.reduce_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
 //           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
		double limit_gh = std::min(parameters.cell_size[0], parameters.cell_size[1])/std::sqrt(parameters.g*P.reduce_maxAbs());

//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
		double dt = parameters.CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

		// provide information to parameters
		parameters.timestep_size = dt;

//		std::cout << parameters.cell_size[0] << " " << parameters.cell_size[1] << std::endl;
#define DELAYED_P_UPDATE	1
#define USE_UP_AND_DOWNWINDING	0

#if DELAYED_P_UPDATE == 0

#if USE_UP_AND_DOWNWINDING == 0
		P += dt*P_t;
#else

		// up/downwinding on staggered grid (FV-like method)
		P += dt/(parameters.cell_size[0]*parameters.cell_size[1])*(
				parameters.cell_size[1]*(
					//             |                       |                       |
					// --v---------|-----------v-----------|-----------v-----------|
					//   h0        u0          h1          u1          h2          u2
					//

					// left edge
					op.shift_right(P*u.return_value_if_positive())		// inflow
					+(P*op.shift_right(u.return_value_if_negative()))	// outflow

					// right edge
					-P*u.return_value_if_positive()					// outflow
					-op.shift_left(P)*u.return_value_if_negative()	// inflow
				)
				+
				parameters.cell_size[0]*(
						// bottom edge
						op.shift_up(P*v.return_value_if_positive())	// inflow
						+(P*op.shift_up(v.return_value_if_negative()))	// outflow

						// top edge
						-P*v.return_value_if_positive()					// outflow
						-op.shift_down(P)*v.return_value_if_negative()	// inflow
				)
			);
#endif
#endif

		u += dt*u_t;
		v += dt*v_t;

#if USE_UP_AND_DOWNWINDING == 0
		// recompute U and V
		U = op.avg_f_x(P)*u;
		V = op.avg_f_y(P)*v;

		P_t = -op.diff_b_x(U) - op.diff_b_y(V);
		P += dt*P_t;
#else

		// up/downwinding on staggered grid (FV-like method)
		P += dt/(parameters.cell_size[0]*parameters.cell_size[1])*(
				parameters.cell_size[1]*(
					//             |                       |                       |
					// --v---------|-----------v-----------|-----------v-----------|
					//   h0        u0          h1          u1          h2          u2
					//

					// left edge
					op.shift_right(P*u.return_value_if_positive())		// inflow
					+(P*op.shift_right(u.return_value_if_negative()))	// outflow

					// right edge
					-P*u.return_value_if_positive()					// outflow
					-op.shift_left(P)*u.return_value_if_negative()	// inflow
				)
				+
				parameters.cell_size[0]*(
						// bottom edge
						op.shift_up(P*v.return_value_if_positive())	// inflow
						+(P*op.shift_up(v.return_value_if_negative()))	// outflow

						// top edge
						-P*v.return_value_if_positive()					// outflow
						-op.shift_down(P)*v.return_value_if_negative()	// inflow
				)
			);
#endif

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
		return P;
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

	SimulationSWEStaggered *simulationSWE = new SimulationSWEStaggered;

	VisSweet<SimulationSWEStaggered> visSweet(simulationSWE);

	delete simulationSWE;

	return 1;
}
