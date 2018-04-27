/*
 * Burgers_Plane_TS_l_irk_mms.cpp
 *
 *  Created on: 16 April 2018
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_l_irk_mms.hpp"

bool Burgers_Plane_TS_l_irk_mms::table_created = false;
int Burgers_Plane_TS_l_irk_mms::table_size = 0;
double** Burgers_Plane_TS_l_irk_mms::table = nullptr;

/*
 * Implementation of an implicit Runge-Kutta method for orders 1 and 2 with the method of manufactured solutions.
 * Solutions are either two sinusoidal waves or a smoothened saw-tooth function (Only working in combination with
 * l_irk_n_sl_mms.
 * Order 1: u^{n+1} = u^n + dt*f(u^{n+1})
 * Order 2: u^{n+1} = u^n + dt*k_1
 * 			k_1 = f(u^n+dt/2 + k_1)
 */
void Burgers_Plane_TS_l_irk_mms::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_l_irk_mms: Only constant time step size allowed");

	if (op.diff_c_y(io_u).reduce_maxAbs()>1e-11)
		FatalError("Burgers_Plane_TS_l_irk_mms: Implemented in 1D only!");


	// setup dummy data
	PlaneData tmp(io_u.planeDataConfig);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_all(0,0);
#endif
	tmp.physical_set_all(0);

	// Setting explicit right hand side and operator of the left hand side
	PlaneData rhs_u = io_u;
//	PlaneData rhs_v = io_v;
	int itime = (int)(simVars.timecontrol.current_simulation_time/simVars.timecontrol.current_timestep_size+0.5)+(int)(i_fixed_dt/simVars.timecontrol.current_timestep_size+0.5);
	if (timestepping_order == 2)
		itime = (int)(2*simVars.timecontrol.current_simulation_time/simVars.timecontrol.current_timestep_size+0.5)+(int)(2*i_fixed_dt/simVars.timecontrol.current_timestep_size+0.5);
	if (second_time)
		itime +=(int)(2*i_fixed_dt/simVars.timecontrol.current_timestep_size+0.5);
	for (std::size_t jj=0; jj<tmp.planeDataConfig->physical_res[1]; jj++)
		for (std::size_t ii=0; ii<tmp.planeDataConfig->physical_res[0]; ii++)
			tmp.p_physical_set(jj,ii,table[itime][ii]);

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{
		PlaneData lhs = io_u;

		if (timestepping_order == 1)
		{
			rhs_u = rhs_u + i_fixed_dt*tmp;
			lhs = ((-i_fixed_dt)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);

			io_u = rhs_u.spectral_div_element_wise(lhs);
//			io_v = rhs_v.spectral_div_element_wise(lhs);
		}
		else if (timestepping_order == 2)
		{
			rhs_u = simVars.sim.viscosity*(op.diff2_c_x(rhs_u) + op.diff2_c_y(rhs_u)) + tmp;
//			rhs_v = simVars.sim.viscosity*(op.diff2_c_x(rhs_v) + op.diff2_c_y(rhs_v));
			lhs = ((-0.5*i_fixed_dt)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);

			PlaneData k1_u = rhs_u.spectral_div_element_wise(lhs);
//			PlaneData k1_v = rhs_v.spectral_div_element_wise(lhs);

			io_u = io_u + i_fixed_dt*k1_u;
//			io_v = io_v + i_fixed_dt*k1_v;

			second_time = !second_time;
		}
		else
			FatalError("This timestepping order is not available with l_irk");

	} else { //Jacobi
		FatalError("NOT available");
	}

}


void Burgers_Plane_TS_l_irk_mms::return_initial(PlaneData& init)
{
	for (std::size_t jj=0; jj<init.planeDataConfig->physical_res[1]; jj++)
		for (std::size_t ii=0; ii<init.planeDataConfig->physical_res[0]; ii++)
			init.p_physical_set(jj,ii,table[0][jj]);
}


void Burgers_Plane_TS_l_irk_mms::setup_look_up_table(double start, double end, double step_size)
{
	if (!table_created)
	{
		// Calculate look-up table
		double tp = 2*M_PIl;
		double x = 0.0;
		double time = 0.0;
		double freq = simVars.sim.f0;
		double eps = 0.1;
		// Calculate number of iterations by rounding the result of max_simulation_time / timestep_size
		int iiMax = (int)((end-start)/step_size+0.5);
		if (timestepping_order == 2)
			iiMax = iiMax*2;
		table_size = iiMax + 1;

		table = new double*[table_size];
		for (int ii=0; ii<table_size; ii++)
		{
			table[ii] = new double[simVars.disc.res_physical[0]];
			time = ii*step_size;
			if ((timestepping_order == 2) && (ii > 0))
				time = time/2-step_size/4;
			for (int jj=0; jj<simVars.disc.res_physical[0]; jj++)
			{
				x = (((double)jj+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];

				if (simVars.setup.benchmark_scenario_id == 8)
				{
					table[ii][jj] = tp * std::sin(tp*x) * std::cos(tp*time) + tp * std::sin(tp*freq*x) * std::cos(tp*freq*time)
						+ (std::sin(tp*x) * std::sin(tp*time) + 1/freq * std::sin(tp*freq*x) * std::sin(tp*freq*time))
						* (tp * std::cos(tp*x) * std::sin(tp*time) + tp * std::cos(tp*freq*x) * std::sin(tp*freq*time))
						- simVars.sim.viscosity * (-tp*tp * std::sin(tp*x) * std::sin(tp*time)
						- tp*tp*freq * std::sin(tp*freq*x) * std::sin(tp*freq*time));
				}
				else if (simVars.setup.benchmark_scenario_id == 12)
				{
					double tmpvar = 0.0;
					double A1 = 0.0;
					double A2 = 0.0;
					double argument = 0.0;
					double AA = 0.0;
					for (int kk=1; kk<freq; kk++)
					{
						argument = tp*kk*x + M_PIl*kk*(1-time);
						AA = eps/sinh(M_PIl*kk*eps/2);
						tmpvar += (simVars.sim.viscosity*std::sin(argument)*4*M_PIl*kk - std::cos(argument))*kk*AA;
						A1 += std::cos(argument)*kk*AA*M_PIl;
						A2 += std::sin(argument)*0.5*AA;
					}
					tmpvar *= 0.5*M_PIl;
					tmpvar += A1*A2;
					table[ii][jj] = tmpvar;
				}
				else
					FatalError("This scenario is not supported");
			}
		}
		table_created = true;
	}
}

/*
 * Setup
 */
void Burgers_Plane_TS_l_irk_mms::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_l_irk_mms::Burgers_Plane_TS_l_irk_mms(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	second_time = false;
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_l_irk_mms::~Burgers_Plane_TS_l_irk_mms()
{
	if (table != nullptr)
	{
		for (int ii=0; ii<table_size; ii++)
			delete[] table[ii];
		delete[] table;
		table = nullptr;
	}
}

