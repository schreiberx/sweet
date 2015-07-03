
#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
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
	DataArray<2> eta_u;
	DataArray<2> eta_v;
	DataArray<2> P_t, u_t, v_t;

	DataArray<2> f;

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
	 *   ____u0,1_____
	 *   |           |
	 *   |           |
	 * v0,0   P0,0   v1,0
	 *   |           |
	 *   |___u0,0____|
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
		eta_u(parameters.res),
		eta_v(parameters.res),
		P_t(parameters.res),
		u_t(parameters.res),
		v_t(parameters.res),

		f(parameters.res),

		op(parameters.cell_size, parameters.res)
	{
		reset();
	}

	void rot_coord(double angle, double &x, double &y)
	{
		angle *= 2.0*M_PI/360.0;
		double nx = std::cos(angle)*x - std::sin(angle)*y;
		double ny = std::sin(angle)*x + std::cos(angle)*y;
		x = nx;
		y = ny;
	}
	void rot_vector(double angle, double &x, double &y)
	{
		angle *= 2.0*M_PI/360.0;
		double nx = std::cos(angle)*x - std::sin(angle)*y;
		double ny = std::sin(angle)*x + std::cos(angle)*y;

		x = nx;
		y = ny;
	}

	void reset()
	{
		parameters.timestep_nr = 0;
		parameters.simulation_time = 0;

		P.data_setall(parameters.h0);
		u.data_setall(0);
		v.data_setall(0);

		double center_x = parameters.init_coord_x;
		double center_y = parameters.init_coord_y;

		if (parameters.setup_scenario == 0)
		{
			std::cout << "Setting up discontinuous radial dam break" << std::endl;

			/*
			 * radial dam break
			 */
			double radius = 100000;
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
					double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

					double dx = x-center_x;
					double dy = y-center_y;

					if (radius*radius >= dx*dx+dy*dy)
						P.getDataRef(j,i) += 1.0;
				}
			}
		}

		if (parameters.setup_scenario == 1)
		{
			std::cout << "Setting up Gaussian radial dam break" << std::endl;

			std::cout << center_x << std::endl;
			std::cout << center_y << std::endl;
			std::cout << std::endl;
			std::cout << parameters.sim_domain_length[0] << std::endl;
			std::cout << parameters.sim_domain_length[1] << std::endl;
			std::cout << std::endl;

			if (	std::abs(center_x/parameters.sim_domain_length[0]-0.5) > 0.001 ||
					std::abs(center_y/parameters.sim_domain_length[1]-0.5) > 0.001
			)
			{
				std::cerr << "ERROR: Gaussian radial dam break has to be centered. Otherwise, there will be discontinuities at the domain border" << std::endl;
				exit(1);
			}

			/*
			 * fun with Gaussian
			 */
			double radius = 1000000;
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double x = ((double)i+0.5)*parameters.sim_domain_length[0]/(double)parameters.res[0];
					double y = ((double)j+0.5)*parameters.sim_domain_length[1]/(double)parameters.res[1];

					double dx = x-center_x;
					double dy = y-center_y;

					dx /= radius;
					dy /= radius;

					P.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
				}
			}
		}

		if (parameters.setup_scenario == 2)
		{
			std::cout << "Setting up balanced steady state solution (ver1)" << std::endl;
			// see doc/balanced_steady_state_solution/*

			if (parameters.sim_f == 0)
			{
				std::cout << "Coriolis force required to setup balanced steady state solution!" << std::endl;
				exit(-1);
			}

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					{
						// P space
						double x = ((double)i+0.5)/(double)parameters.res[0];
						P.getDataRef(j,i) = std::sin(2.0*M_PI*x) + parameters.h0;
					}

					{
						// v space
						double x = ((double)i+0.5)/(double)parameters.res[0];

						v.getDataRef(j,i) = 2.0*M_PI*std::cos(2.0*M_PI*x)/parameters.sim_f;
					}
				}
			}
		}


		if (parameters.setup_scenario == 3)
		{
			std::cout << "Setting up balanced steady state solution (ver2)" << std::endl;
			// see doc/balanced_steady_state_solution/*

			if (parameters.sim_f == 0)
			{
				std::cout << "Coriolis force required to setup balanced steady state solution!" << std::endl;
				exit(-1);
			}

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					{
						// P space
						double y = ((double)j+0.5)/(double)parameters.res[1];
						P.getDataRef(j,i) = std::sin(2.0*M_PI*y)+parameters.h0;
					}
					{
						// u space
						double y = ((double)j+0.5)/(double)parameters.res[1];

						u.getDataRef(j,i) = -2.0*M_PI*std::cos(2.0*M_PI*y)/parameters.sim_f;
					}
				}
			}
		}


		if (parameters.setup_scenario == 4)
		{
			std::cout << "Setting up balanced steady state solution (ver3)" << std::endl;
			std::cout << "The rotations create instabilities in the linear terms, hence this does not work, yet" << std::endl;
			std::cout << "Furthermore, rotation doesn't make sense due to the missing periodic boundaries (discontinuities occur)" << std::endl;
			exit(-1);

			// see doc/balanced_steady_state_solution/*

			if (parameters.sim_f == 0)
			{
				std::cout << "Coriolis force required to setup balanced steady state solution!" << std::endl;
				exit(-1);
			}

			u.data_setall(0);
			v.data_setall(0);

			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					double angle = 45;

					{
						// P space
						double x = ((double)i+0.5)/(double)parameters.res[0];
						double y = ((double)j+0.5)/(double)parameters.res[1];
						rot_coord(angle, x, y);

						P.getDataRef(j,i) = std::sin(2.0*M_PI*x) + parameters.h0;
					}

					{
						// u space
						double x = ((double)i)/(double)parameters.res[0];
						double y = ((double)j+0.5)/(double)parameters.res[1];
						double xo = x;	double yo = y;
						rot_coord(angle, xo, yo);

						double vel_u = 0;
						double vel_v = 2.0*M_PI*std::cos(2.0*M_PI*x)/parameters.sim_f;
						double vel_x_o = vel_u+x;	double vel_y_o = vel_v+y;
						rot_coord(angle, vel_x_o, vel_y_o);

						double fin_vel_u = vel_x_o - xo;
						double fin_vel_v = vel_y_o - yo;

						u.getDataRef(j,i) = fin_vel_u;
					}
					{
						// v space
						double x = ((double)i+0.5)/(double)parameters.res[0];
						double y = ((double)j)/(double)parameters.res[1];
						double xo = x;	double yo = y;
						rot_coord(angle, xo, yo);

						double vel_u = 0;
						double vel_v = 2.0*M_PI*std::cos(2.0*M_PI*x)/parameters.sim_f;
						double vel_x_o = vel_u+x;	double vel_y_o = vel_v+y;
						rot_coord(angle, vel_x_o, vel_y_o);

						vel_u = vel_x_o - xo;
						vel_v = vel_y_o - yo;

						v.getDataRef(j,i) = vel_v;
//						v.getDataRef(j,i) = 2.0*M_PI*std::cos(2.0*M_PI*x)/parameters.sim_f;
					}
				}
			}
		}


		if (parameters.setup_scenario == 5)
		{
			for (std::size_t j = 0; j < parameters.res[1]; j++)
			{
				for (std::size_t i = 0; i < parameters.res[0]; i++)
				{
					{
						// beta plane
						double x = ((double)i)/(double)parameters.res[0];
						double y = ((double)j)/(double)parameters.res[1];
						f.getDataRef(j,i) = (y-0.5)*(y-0.5)*parameters.sim_f;
					}


					{
						double x = ((double)i+0.5)/(double)parameters.res[0];
						double y = ((double)j+0.5)/(double)parameters.res[1];

						double dx = x-center_x;
						double dy = y-center_y;

						P.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
					}
				}
			}
		}
	}



	bool run_timestep()
	{
		/*
		 * Note, that this grid does not follow the one in the paper
		 *
		 *             ^
		 *             |
		 *       ____v0,1_____
		 *       |           |
		 *       |           |
		 * <- u0,0  H/P0,0   u1,0 ->
		 *       |           |
		 * eta0,0|___v0,0____|
		 *             |
		 *             V
		 */

		U = op.avg_b_x(P)*u;
		V = op.avg_b_y(P)*v;

		H = P + 0.5*(op.avg_f_x(u*u) + op.avg_f_y(v*v));

		if (parameters.setup_scenario != 5)
		{
			eta_u = (op.diff_b_x(v) - op.diff_b_y(u) + parameters.sim_f/parameters.sim_domain_length[0]) / op.avg_b_x(op.avg_b_y(P));
			eta_v = (op.diff_b_x(v) - op.diff_b_y(u) + parameters.sim_f/parameters.sim_domain_length[1]) / op.avg_b_x(op.avg_b_y(P));
		}
		else
		{
			eta_u = (op.diff_b_x(v) - op.diff_b_y(u) + f) / op.avg_b_x(op.avg_b_y(P));
			eta_v = eta_u;
		}

		double normalization = (parameters.sim_domain_length[0]*parameters.sim_domain_length[1]) / ((double)parameters.res[0]*(double)parameters.res[1]);
		// mass
		parameters.mass = P.reduce_sum() * normalization;

		// energy
		parameters.energy = 0.5*(
				P*P +
				P*op.avg_f_x(u*u) +
				P*op.avg_f_y(v*v)
			).reduce_sum() * normalization;

		// potential enstropy
		parameters.potential_entrophy = 0.5*(eta_u*eta_u*op.avg_b_x(op.avg_b_y(P))).reduce_sum() * normalization;

		u_t = op.avg_f_y(eta_u*op.avg_b_x(V)) - op.diff_b_x(H);
		v_t = -op.avg_f_x(eta_v*op.avg_b_y(U)) - op.diff_b_y(H);

//		std::cout << u_t << std::endl;


		if (parameters.verbosity > 0)
		{
			if (parameters.timestep_nr == 0)
			{
				std::cout << "MASS\tENERGY\tPOT_ENSTROPHY";

				if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
					std::cout << "\tABS_U_DT";

				std::cout << std::endl;
			}

			std::cout << parameters.mass << "\t" << parameters.energy << "\t" << parameters.potential_entrophy;

			// this should be zero for the steady state test
			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				double test_val = u_t.reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				std::cout << "\t" << test_val;

				test_val = v_t.reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				std::cout << "\t" << test_val;

				test_val = P_t.reduce_sumAbs() / (double)(parameters.res[0]*parameters.res[1]);
				std::cout << "\t" << test_val;
			}

			std::cout << std::endl;
		}

		if (parameters.sim_viscocity != 0)
		{
			u_t += (op.diff2_c_x(u) + op.diff2_c_x(v))*parameters.sim_viscocity;
			v_t += (op.diff2_c_y(u) + op.diff2_c_y(v))*parameters.sim_viscocity;
		}

		if (parameters.sim_hyper_viscocity != 0)
		{
			u_t += (op.diff2_c_x(op.diff2_c_x(u)) + op.diff2_c_x(op.diff2_c_x(v)))*parameters.sim_hyper_viscocity;
			v_t += (op.diff2_c_y(op.diff2_c_y(u)) + op.diff2_c_y(op.diff2_c_y(v)))*parameters.sim_hyper_viscocity;
		}

		double limit_speed = std::max(parameters.cell_size[0]/u.reduce_maxAbs(), parameters.cell_size[1]/v.reduce_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
//           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
//		double limit_gh = limit_speed;
		double limit_gh = std::min(parameters.cell_size[0], parameters.cell_size[1])/std::sqrt(parameters.g*P.reduce_maxAbs());

//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
		double dt = parameters.sim_CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

		// provide information to parameters
		parameters.timestep_size = dt;


		bool leapfrog_like_update = true;
		bool up_and_downwinding = false;


		// standard explicit euler on velocities
		u += dt*u_t;
		v += dt*v_t;

		if (!leapfrog_like_update && !up_and_downwinding)
		{
			P_t = -op.diff_f_x(U) - op.diff_f_y(V);
			P += dt*P_t;
		}
		else
		{
			if (!up_and_downwinding)
			{
				if (leapfrog_like_update)
				{
					// recompute U and V
					U = op.avg_b_x(P)*u;
					V = op.avg_b_y(P)*v;
				}

				P_t = -op.diff_f_x(U) - op.diff_f_y(V);
				P += dt*P_t;
			}
			else
			{

				//             |                       |                       |
				// --v---------|-----------v-----------|-----------v-----------|
				//   h-1       u0          h0          u1          h1          u2
				//

				// same a above, but formulated in a finite-difference style
				P += dt
					*(
						(
							// u is positive
							op.shift_right(P)*u.return_value_if_positive()	// inflow
							-P*op.shift_left(u.return_value_if_positive())					// outflow

							// u is negative
							+(P*u.return_value_if_negative())	// outflow
							-op.shift_left(P*u.return_value_if_negative())		// inflow
						)*(1.0/parameters.cell_size[0])	// here we see a finite-difference-like formulation
						+
						(
							// v is positive
							op.shift_up(P)*v.return_value_if_positive()		// inflow
							-P*op.shift_down(v.return_value_if_positive())					// outflow

							// v is negative
							+(P*v.return_value_if_negative())	// outflow
							-op.shift_down(P*v.return_value_if_negative())	// inflow
						)*(1.0/parameters.cell_size[1])
					);
			}
		}

		parameters.simulation_time += dt;
		parameters.timestep_nr++;


		if (parameters.max_timesteps != -1 && parameters.max_timesteps <= parameters.timestep_nr)
			return false;

		if (parameters.max_simulation_time != -1 && parameters.max_simulation_time <= parameters.simulation_time)
			return false;

		return true;
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
		case 0:	*o_dataArray = &P;		break;
		case 1: *o_dataArray = &u;		break;
		case 2: *o_dataArray = &v;		break;
		case 3: *o_dataArray = &H;		break;
		case 4: *o_dataArray = &eta_u;	break;
		case 5: *o_dataArray = &U;		break;
		case 6: *o_dataArray = &V;		break;
		}

		*o_aspect_ratio = parameters.sim_domain_length[1] / parameters.sim_domain_length[0];
	}

	const char* vis_get_status_string()
	{
		const char *vis_id_strings[] =
		{
				"P", "u", "v", "eta", "H", "U", "V"
		};
		int vis_id = parameters.vis_id % 7;

		static char title_string[1024];
		sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %e, Vis: %s, Mass: %e, Energy: %e, Potential Entrophy: %e",
				parameters.simulation_time,
				parameters.simulation_time/(60.0*60.0*24.0),
				parameters.timestep_nr,
				parameters.timestep_size,
				vis_id_strings[vis_id],
				parameters.mass, parameters.energy, parameters.potential_entrophy);
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
			if (parameters.vis_id >= 7)
				parameters.vis_id = 6;
			break;

		case 'V':
			parameters.vis_id--;
			if (parameters.vis_id < 0)
				parameters.vis_id = 0;
			break;
		}
	}
};




int main(int i_argc, char *i_argv[])
{
	parameters.setup(i_argc, i_argv);

	SimulationSWEStaggered *simulationSWE = new SimulationSWEStaggered;

#if SWEET_GUI
	VisSweet<SimulationSWEStaggered> visSweet(simulationSWE);
#else
	simulationSWE->run();
#endif

	delete simulationSWE;

	return 1;
}
