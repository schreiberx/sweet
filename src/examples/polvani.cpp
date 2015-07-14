
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif


#include <sweet/DataArray.hpp>
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationParameters.hpp>
#include <sweet/TimesteppingRK.hpp>
#include <sweet/SWEValidationBenchmarks.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/Stopwatch.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>


SimulationParameters parameters;



DataArray<2> *init_h = nullptr;
DataArray<2> *init_u = nullptr;
DataArray<2> *init_v = nullptr;


class SimulationSWE
{
public:
	DataArray<2> prog_h, prog_u, prog_v;
	DataArray<2> eta;
	DataArray<2> tmp;

	Operators2D op;

	TimesteppingRK timestepping;

	int last_timestep_nr_update_diagnostics = -1;

	double benchmark_diff_h;
	double benchmark_diff_u;
	double benchmark_diff_v;

public:
	SimulationSWE()	:
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

		benchmark_diff_h = 0;
		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		parameters.status_timestep_nr = 0;
		parameters.status_simulation_time = 0;

#if 0
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

#else

		if (init_h != nullptr)
		{
			prog_h = *init_h;
			prog_u = *init_u;
			prog_v = *init_v;
		}
#endif
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
			o_u_t += (op.diff2_c_x(i_u) + op.diff2_c_x(i_v))*parameters.sim_viscocity;
			o_v_t += (op.diff2_c_y(i_u) + op.diff2_c_y(i_v))*parameters.sim_viscocity;
		}

		if (parameters.sim_hyper_viscocity != 0)
		{
			o_u_t += (op.diff2_c_x(op.diff2_c_x(i_u)) + op.diff2_c_x(op.diff2_c_x(i_v)))*parameters.sim_hyper_viscocity;
			o_v_t += (op.diff2_c_y(op.diff2_c_y(i_u)) + op.diff2_c_y(op.diff2_c_y(i_v)))*parameters.sim_hyper_viscocity;
		}


		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt > 0)
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
				double limit_visc = std::numeric_limits<double>::infinity();
		//        if (viscocity > 0)
		//           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

				// limit by gravitational acceleration
				double limit_gh = std::min(parameters.sim_cell_size[0], parameters.sim_cell_size[1])/std::sqrt(parameters.sim_g*i_h.reduce_maxAbs());

				if (parameters.verbosity > 2)
					std::cerr << "limit_speed: " << limit_speed << ", limit_visc: " << limit_visc << ", limit_gh: " << limit_gh << std::endl;

				o_dt = parameters.sim_CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);
			}
		}

		if (!parameters.timestepping_leapfrog_like_update)
		{
			if (!parameters.timestepping_up_and_downwinding)
			{
				// standard update
				o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
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

				if (parameters.sim_viscocity != 0)
					o_h_t += (op.diff2_c_x(i_h) + op.diff2_c_y(i_h))*parameters.sim_viscocity;

				if (parameters.sim_hyper_viscocity != 0)
					o_h_t += (op.diff2_c_x(op.diff2_c_x(i_h)) + op.diff2_c_y(op.diff2_c_y(i_h)))*parameters.sim_hyper_viscocity;
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
			}
		}

		if (parameters.sim_potential_viscocity != 0)
			o_h_t += (op.diff2_c_x(i_h) + op.diff2_c_y(i_h))*parameters.sim_potential_viscocity;

		if (parameters.sim_potential_hyper_viscocity != 0)
			o_h_t += (op.diff2_c_x(op.diff2_c_x(i_h)) + op.diff2_c_y(op.diff2_c_y(i_h)))*parameters.sim_potential_hyper_viscocity;

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
					o_ostream << "\tDIFF_P\tDIFF_U\tDIFF_V";

				o_ostream << std::endl;
			}

			o_ostream << parameters.status_simulation_time << "\t" << parameters.diagnostics_mass << "\t" << parameters.diagnostics_energy << "\t" << parameters.diagnostics_potential_entrophy;

			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// h
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_h(parameters, x, y));
					}

				benchmark_diff_h = (prog_h-tmp).reduce_sumAbs_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_h;

				// set data to something to overcome assertion error
				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// u space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_u(parameters, x, y));
					}

				benchmark_diff_u = (prog_u-tmp).reduce_sumAbs_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_u;

				for (std::size_t j = 0; j < parameters.res[1]; j++)
					for (std::size_t i = 0; i < parameters.res[0]; i++)
					{
						// v space
						double x = (((double)i+0.5)/(double)parameters.res[0])*parameters.sim_domain_length[0];
						double y = (((double)j+0.5)/(double)parameters.res[1])*parameters.sim_domain_length[1];

						tmp.set(j, i, SWEValidationBenchmarks::return_v(parameters, x, y));
					}

				benchmark_diff_v = (prog_v-tmp).reduce_sumAbs_quad() / (double)(parameters.res[0]*parameters.res[1]);
				o_ostream << "\t" << benchmark_diff_v;
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
	void vis_post_frame_processing(
			int i_num_iterations
	)
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
		DataArray<2> &o = (DataArray<2> &)*vis_arrays[id].data;
		o.setAllSpec(0, 0);
		o.setSpec(0, 0, 1, 0);

		if (id == 0)
		{
			o.setSpec(0, 0, 1, 0);
		}
		else if (id == 1)
		{
			o.setSpec(0, 0, 1, 1);
		}
		else if (id == 2)
		{
			o.setSpec(0, 1, 1, 0);
		}
		else if (id == 3)
		{
			o.setSpec(0, 1, 1, 1);
//			o.setSpec(0, o.resolution_spec[0]-1, 1, 0);
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



void compute_polvani_initialization(
		DataArray<2> &h,
		DataArray<2> &u,
		DataArray<2> &v
)
{

	{
		double m = 25;
		double k0 = 14;

		double R = 0.01;
		double F = 0.04;

		if (parameters.bogus_var0 > 0)
			R = parameters.bogus_var0;

		if (parameters.bogus_var1 > 0)
			F = parameters.bogus_var1;

		std::cout << "Using R of " << R << std::endl;
		std::cout << "Using F of " << F << std::endl;

		double B = (R*R)/(F*F);

		// initialize seed
	//	srand(0x15051982);
		srand(time(NULL));

		Operators2D op(parameters.res, parameters.sim_domain_length, parameters.use_spectral_diffs);

		/*
		 * see Polvani et. al: Coherent structures of shallow-water turbulence
		 */

		/**
		 * Helpers
		 */
		DataArray<2> laplace_op = op.diff2_c_x+op.diff2_c_y;

		/**
		 * STEP 1) Initialize PSI
		 */

		DataArray<2> energy_init(parameters.res);

		for (std::size_t j = 0; j < energy_init.resolution_spec[1]; j++)
		{
			for (std::size_t i = 0; i < energy_init.resolution_spec[0]; i++)
			{
				std::size_t ka = i;
				std::size_t kb = (j < energy_init.resolution_spec[1]/2 ? j : energy_init.resolution_spec[1]-j);
				assert(kb >= 0 && kb <= energy_init.resolution_spec[1]/2);

				double k = std::sqrt((double)(ka*ka)+(double)(kb*kb));

				// compute energy spectrum
				double energy = pow(k, m*0.5)/pow(k+k0, m);
	//			double energy = exp(-pow(k-k0,2));

				assert(energy >= 0);

				energy_init.setSpec(j, i, energy, 0);
			}
		}

		/**
		 * See paper:
		 * The energy should be set, so that the rms speed of the initial velocity is 1:
		 *
		 * non-dimensional RMS speed = 1
		 *
		 * => rms speed = sqrt(1/N  \sum( sqrt(u^2+v^2) )) = 1
		 *
		 * rms(u) = rms(v) = 1/sqrt(2)
		 *
		 * => energy: |E| = 0.5*sqrt(1/N * 1)
		 */

		// compute rms
		double scale_energy = 1.0/energy_init.reduce_rms();;

		// normalize to get desired velocity
		scale_energy *= 0.5;

		energy_init *= scale_energy;

		std::cout << "ENERGY rms: " << energy_init.reduce_rms() << std::endl;

		DataArray<2> psi(parameters.res);

		for (std::size_t j = 0; j < psi.resolution_spec[1]; j++)
		{
			for (std::size_t i = 0; i < psi.resolution_spec[0]; i++)
			{
				std::size_t ka = i;
				std::size_t kb = (j < psi.resolution_spec[1]/2 ? j : psi.resolution_spec[1]-j);
				assert(kb >= 0 && kb <= energy_init.resolution_spec[1]/2);

				double k = std::sqrt((double)(ka*ka)+(double)(kb*kb));

				// compute energy spectrum
				double energy = energy_init.getSpec_Re(j, i);

				if (k == 0)
				{
					psi.setSpec(j, i, 0, 0);
					continue;
				}

				double psi_abs = std::sqrt(energy*2.0/(k*k));

				// compute random number \in [0;1[
				double r = (double)rand()/((double)RAND_MAX+1);

				// generate random phase shift
				double psi_re = cos(2.0*M_PIl*r);
				double psi_im = sin(2.0*M_PIl*r);

				psi_re *= psi_abs;
				psi_im *= psi_abs;

				psi.setSpec(j, i, psi_re, psi_im);
			}
		}

		std::cout << "CART Max: " << psi.reduce_maxAbs() << std::endl;
		std::cout << "SPEC Max: " << psi.reduce_spec_maxAbs() << std::endl;
		std::cout << "SPEC Centroid: " << psi.reduce_getFrequencyCentroid() << std::endl;

		/**
		 * STEP 2) Compute h
		 */
		// equation 2.5b
		h = op.laplace(psi)+(2.0*R)*op.arakawa_jacobian(op.diff_c_x(psi), op.diff_c_y(psi));
		// (D*D)^-1
		h = h.spec_div_element_wise(laplace_op);

		// Add 1.0 for non-dimensional average height (not sure, if this is correct)
//		h += 1.0;


		/**
		 * STEP 3) Solve with iterations for \xi
		 */
		DataArray<2> xi(parameters.res);
		xi.setAllSpec(0, 0);

		DataArray<2> prev_xi(parameters.res);
		prev_xi.setAllSpec(0, 0);

		double prev_inf_norm = -1;
		double eps = 1e-9;


		int i = 0;
		for (; i < 100; i++)
		{
			// cache this value
			DataArray<2> laplace_psi = op.laplace(psi);

			/**
			 * we first solve psi_t for equation (2.5a)
			 */

			DataArray<2> psi_t =
				(
					op.arakawa_jacobian(psi, laplace_psi)
					+ (1.0/R)*op.laplace(xi)
					+ op.diff_c_x(laplace_psi*op.diff_c_x(xi))
					+ op.diff_c_y(laplace_psi*op.diff_c_y(xi))
				).spec_div_element_wise(laplace_op);


			////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////
			DataArray<2> psi_t_rhs =
				(
					op.arakawa_jacobian(psi, laplace_psi)
					+ (1.0/R)*op.laplace(xi)
					+ op.diff_c_x(laplace_psi*op.diff_c_x(xi))
					+ op.diff_c_y(laplace_psi*op.diff_c_y(xi))
				);

			DataArray<2> psi_t_lhs = op.laplace(psi_t);

			double check = (psi_t_lhs-psi_t_rhs).reduce_maxAbs();
			std::cout << check << std::endl;
			////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////

			/**
			 * Solve iterative equation (after equation 2.6)
			 */
			DataArray<2> rhs(parameters.res);
			rhs.setAllSpec(0, 0);

			// RHS, 1st line
			rhs += -op.arakawa_jacobian(psi, op.laplace(xi));
			rhs += op.laplace(op.arakawa_jacobian(psi, h));

			// RHS, 2nd line
			// special jacobian, see Footnote 2 on page 186
			DataArray<2> spec_jacob =
					op.diff2_c_x(psi_t)*op.diff2_c_y(psi)
					+ op.diff2_c_x(psi)*op.diff2_c_y(psi_t)
					- 2.0*op.diff_c_x(op.diff_c_y(psi))*op.diff_c_x(op.diff_c_y(psi_t))
				;
			rhs += 2.0*R*spec_jacob;
			rhs += -op.diff_c_x(laplace_psi*op.diff_c_x(xi)) - op.diff_c_y(laplace_psi*op.diff_c_y(xi));

			// RHS, 3rd line
			rhs += op.laplace(
					op.diff_c_x(h*op.diff_c_x(xi)) +
					op.diff_c_y(h*op.diff_c_y(xi))
				);

			// LHS, 1st line
			DataArray<2> lhs(parameters.res);
			lhs.setAllSpec(0, 0);

			// IMPORTANT! Apply last operator element-wise
			lhs += (1.0/R)*(1.0-B*laplace_op).operator()(laplace_op);

			xi = rhs.spec_div_element_wise(lhs);

			double new_inf_norm;
			// convergence check
			if (0)
			{
				// compute difference on both sides of equations in xi

				DataArray<2> psi_t =
					(
						op.arakawa_jacobian(psi, laplace_psi)
						+ (1.0/R)*op.laplace(xi)
						+ op.diff_c_x(laplace_psi*op.diff_c_x(xi))
						+ op.diff_c_y(laplace_psi*op.diff_c_y(xi))
					).spec_div_element_wise(laplace_op);

				DataArray<2> spec_jacob =
						op.diff2_c_x(psi_t)*op.diff2_c_y(psi)
						+ op.diff2_c_x(psi)*op.diff2_c_y(psi_t)
						- 2.0*op.diff_c_x(op.diff_c_y(psi))*op.diff_c_x(op.diff_c_y(psi_t))
					;

				new_inf_norm =
						(
							(1.0/R)*(op.laplace(xi)-B*op.laplace(op.laplace(xi)))	/// LHS
							-
							(-1.0)*op.arakawa_jacobian(psi, op.laplace(xi))
							+ op.laplace(op.arakawa_jacobian(psi, h))
							+ 2.0*R*op.arakawa_jacobian(op.diff_c_x(psi), op.diff_c_y(psi))
								- op.diff_c_x(op.laplace(psi)*op.diff_c_x(xi))
								- op.diff_c_y(op.laplace(psi)*op.diff_c_y(xi))
							+ op.laplace(
									op.diff_c_x(h*op.diff2_c_x(xi))
									+ op.diff_c_y(h*op.diff2_c_y(xi))
								)
						).reduce_maxAbs();//*normalize_norm;
			}
			else
			{
				new_inf_norm = (prev_xi-xi).reduce_maxAbs();
			}

			double normalization = xi.reduce_sumAbs()/(double)(parameters.res[0]*parameters.res[1]);

			double convergence = std::abs(prev_inf_norm - new_inf_norm)/normalization;
			std::cout << "INF norm / convergence: " << new_inf_norm << " / " << convergence << std::endl;

			if (convergence < eps)
			{
				std::cout << "CONVERGED with norm: " << new_inf_norm << std::endl;
				break;
			}
			prev_inf_norm = new_inf_norm;
			prev_xi = xi;
		};

		if (i == 100)
		{
			std::cerr << "NO CONVERGENCE!" << std::endl;
			exit(1);
		}

		double epsilon = R * (std::max(1.0,R)/std::max(1.0,B));

		u = op.diff_c_y(psi) + epsilon*op.diff_c_x(xi);
		v = -op.diff_c_x(psi) + epsilon*op.diff_c_y(xi);

		double rms_u = u.reduce_rms();
		u /= rms_u*sqrt(2.0);
	//	rms_u = u.reduce_rms();

		double rms_v = v.reduce_rms();
		v /= rms_v*sqrt(2.0);
	//	rms_v = v.reduce_rms();

		std::cout << "Target RMS u, v: " << 1.0/sqrt(2.0) << ", " << 1.0/sqrt(2.0) << std::endl;
		std::cout << "Current RMS u, v: " << rms_u << ", " << rms_v << std::endl;

		/**
		 * Check results
		 */
		{
//			double normalize_norm = 1.0/(double)(parameters.res[0]*parameters.res[1]);

			DataArray<2> laplace_psi = op.laplace(psi);
			double psi_l1_norm =
					(	op.laplace(psi)
						-
						op.arakawa_jacobian(psi, op.laplace(psi))
						+(1.0/R)*op.laplace(xi)
						+op.diff_c_x(op.laplace(psi)*op.diff_c_x(xi))
						+op.diff_c_y(op.laplace(psi)*op.diff_c_y(xi))
					).reduce_maxAbs();//*normalize_norm;
			std::cout << "PSI l1 norm: " << psi_l1_norm << std::endl;

			double h_l1_norm =
					(
						op.laplace(h)
						-
						op.laplace(psi)
						+2.0*R*op.arakawa_jacobian(op.diff_c_x(psi), op.diff_c_y(psi))
					).reduce_maxAbs();//*normalize_norm;
			std::cout << "H l1 norm: " << h_l1_norm << std::endl;


			DataArray<2> psi_t =
				(
					op.arakawa_jacobian(psi, laplace_psi)
					+ (1.0/R)*op.laplace(xi)
					+ op.diff_c_x(laplace_psi*op.diff_c_x(xi))
					+ op.diff_c_y(laplace_psi*op.diff_c_y(xi))
				).spec_div_element_wise(laplace_op);

			DataArray<2> spec_jacob =
					op.diff2_c_x(psi_t)*op.diff2_c_y(psi)
					+ op.diff2_c_x(psi)*op.diff2_c_y(psi_t)
					- 2.0*op.diff_c_x(op.diff_c_y(psi))*op.diff_c_x(op.diff_c_y(psi_t))
				;

			double xi_l1_norm =
					(
						(1.0/R)*(op.laplace(xi)-B*op.laplace(op.laplace(xi)))	/// LHS
						-
						(-1.0)*op.arakawa_jacobian(psi, op.laplace(xi))
						+ op.laplace(op.arakawa_jacobian(psi, h))
						+ 2.0*R*op.arakawa_jacobian(op.diff_c_x(psi), op.diff_c_y(psi))
							- op.diff_c_x(op.laplace(psi)*op.diff_c_x(xi))
							- op.diff_c_y(op.laplace(psi)*op.diff_c_y(xi))
						+ op.laplace(
								op.diff_c_x(h*op.diff2_c_x(xi))
								+ op.diff_c_y(h*op.diff2_c_y(xi))
							)
					).reduce_maxAbs();//*normalize_norm;
			std::cout << "xi l1 norm: " << xi_l1_norm << std::endl;
		}

		DataArray<2> energy = 0.5*(u*u+v*v);
		std::cout << "ENERGY Centroid: " << energy.reduce_getFrequencyCentroid() << std::endl;
		std::cout << std::endl;
		std::cout << "R: " << R << std::endl;
		std::cout << "F: " << F << std::endl;
		std::cout << "B: " << B << std::endl;
		std::cout << "XIrms/PSIrms: " << xi.reduce_rms()/psi.reduce_rms() << std::endl;
	}

}


int main(int i_argc, char *i_argv[])
{
	std::cout << std::setprecision(14);
	std::cerr << std::setprecision(14);

	parameters.setup(i_argc, i_argv);

	DataArray<2> local_h(parameters.res);
	DataArray<2> local_u(parameters.res);
	DataArray<2> local_v(parameters.res);

	if (1)
	{
		compute_polvani_initialization(local_h, local_u, local_v);

		::init_h = &local_h;
		::init_u = &local_u;
		::init_v = &local_v;
	}


	if (1)
	{
		std::cout << std::setprecision(14);
		std::cerr << std::setprecision(14);

		SimulationSWE *simulationSWE = new SimulationSWE;

		std::ostringstream buf;
		buf << std::setprecision(14);

	#if SWEET_GUI
		VisSweet<SimulationSWE> visSweet(simulationSWE);
	#else
		simulationSWE->reset();

		Stopwatch time;
		time.reset();

		double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

		if (parameters.verbosity > 1)
		{
			simulationSWE->update_diagnostics();
			diagnostics_energy_start = parameters.diagnostics_energy;
			diagnostics_mass_start = parameters.diagnostics_mass;
			diagnostics_potential_entrophy_start = parameters.diagnostics_potential_entrophy;
		}

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

		time.stop();

		double seconds = time();

		std::cout << "Simulation time: " << seconds << " seconds" << std::endl;
		std::cout << "Time per time step: " << seconds/(double)parameters.status_timestep_nr << " sec/ts" << std::endl;



		if (parameters.verbosity > 1)
		{
			std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((parameters.diagnostics_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
			std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((parameters.diagnostics_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((parameters.diagnostics_potential_entrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start) << std::endl;

			if (parameters.setup_scenario == 2 || parameters.setup_scenario == 3 || parameters.setup_scenario == 4)
			{
				std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark_diff_h << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark_diff_u << std::endl;
				std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark_diff_v << std::endl;
			}
		}
	#endif

		delete simulationSWE;
	}

	return 1;
}
