/*
 * Burgers_Plane_TS_ln_imex_mms.cpp
 *
 *  Created on: 17 April 2018
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_imex_mms.hpp"

bool Burgers_Plane_TS_ln_imex_mms::table_created = false;
int Burgers_Plane_TS_ln_imex_mms::table_size = 0;
double** Burgers_Plane_TS_ln_imex_mms::table = nullptr;

/*
 * Implementation of the IMEX method applied to Burgers' equation with the method of manufacturing solutions according
 * to the algorithm provided in
 * Ascher et al. (1997) Implicit-Explicit Runge-Kutta Methods for Time-Dependent Partial Differential Equations
 * for a timestep from u^n to u^{n+1}. The manufactured solutions are either two sinusoidal waves or a smoothened saw-tooth
 * function.
 *
 * First order implementation is IMEX(1,1,1)
 * 	(The IMEX(1,2,1) implementation is available as-well, see below)
 * Second order implementation is IMEX(1,2,2)
 */
void Burgers_Plane_TS_ln_imex_mms::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.misc.verbosity > 2)
		std::cout << "Burgers_Plane::run_timestep_imex_mms()" << std::endl;

	// Use local variables for calculation
	PlaneData u=io_u;
	PlaneData v=io_v;
	double dt = i_fixed_dt;

	// Setup dummy data
	PlaneData tmp(io_u.planeDataConfig);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
	tmp.spectral_set_all(0,0);
#endif
	tmp.physical_set_all(0);

	int itime = (int)(simVars.timecontrol.current_simulation_time/simVars.timecontrol.current_timestep_size+0.5)+(int)(dt/simVars.timecontrol.current_timestep_size+0.5);
	for (std::size_t jj=0; jj<tmp.planeDataConfig->physical_res[1]; jj++)
		for (std::size_t ii=0; ii<tmp.planeDataConfig->physical_res[0]; ii++)
			tmp.p_physical_set(jj,ii,table[itime][ii]);

	// Initialize variables for right-hand side
	PlaneData rhs_u = u;
	PlaneData rhs_v = v;

	// Calculate \hat{k}_1 = f(u^n)
	PlaneData explK_u = -(u*op.diff_c_x(u)+v*op.diff_c_y(u));
	PlaneData explK_v = -(u*op.diff_c_x(v)+v*op.diff_c_y(v));

	// Calculate k_1
	if (timestepping_order == 1)
	{
		rhs_u = simVars.sim.viscosity*(op.diff2_c_x(u+dt*explK_u)+op.diff2_c_y(u+dt*explK_u))+tmp;
		rhs_v = simVars.sim.viscosity*(op.diff2_c_x(v+dt*explK_v)+op.diff2_c_y(v+dt*explK_v));
	}
	else if (timestepping_order == 2)
	{
		rhs_u = simVars.sim.viscosity*(op.diff2_c_x(u+dt/2*explK_u)+op.diff2_c_y(u+dt/2*explK_u))+tmp;
		rhs_v = simVars.sim.viscosity*(op.diff2_c_x(v+dt/2*explK_v)+op.diff2_c_y(v+dt/2*explK_v));
	}
	else
		FatalError("The chosen timestepping-order is not possible with IMEX");

	if (simVars.disc.use_spectral_basis_diffs) //spectral
	{

		PlaneData lhs = u;
		if (timestepping_order == 1)
		{
			lhs = ((-dt)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		else
		{
			lhs = ((-dt*0.5)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
		}
		PlaneData implK_u = rhs_u.spectral_div_element_wise(lhs);
		PlaneData implK_v = rhs_v.spectral_div_element_wise(lhs);

		// Calculate \hat{k}_2
		if (timestepping_order == 2)
		{
			explK_u = - ((u+dt/2*explK_u+dt/2*implK_u)*op.diff_c_x(u+dt/2*explK_u+dt/2*implK_u)
					+ (v+dt/2*explK_v+dt/2*implK_v)*op.diff_c_y(u+dt/2*explK_u+dt/2*implK_u));
			explK_v = - ((u+dt/2*explK_u+dt/2*implK_u)*op.diff_c_x(v+dt/2*explK_v+dt/2*implK_v)
					+ (v+dt/2*explK_v+dt/2*implK_v)*op.diff_c_y(v+dt/2*explK_v+dt/2*implK_v));
		}
#if 0 //used for Ascher(1,2,1) from Ascher et al. (1997)
		else
		{
			explK_u = - ((u+dt*explK_u+dt*implK_u)*op.diff_c_x(u+dt*explK_u+dt*implK_u)
					+ (v+dt*explK_v+dt*implK_v)*op.diff_c_y(u+dt*explK_u+dt*implK_u));
			explK_v = - ((u+dt*explK_u+dt*implK_u)*op.diff_c_x(v+dt*explK_v+dt*implK_v)
					+ (v+dt*explK_v+dt*implK_v)*op.diff_c_y(v+dt*explK_v+dt*implK_v));
		}
#endif

		// Calculate u^{n+1} = u^n + dt*k_1 + dt*\hat{k}_2
		io_u = u + dt*explK_u + dt*implK_u;
		io_v = v + dt*explK_v + dt*implK_v;

	} else { //Jacobi
		FatalError("NOT available");
	}
}


void Burgers_Plane_TS_ln_imex_mms::return_initial(PlaneData& init)
{
	for (std::size_t jj=0; jj<init.planeDataConfig->physical_res[1]; jj++)
		for (std::size_t ii=0; ii<init.planeDataConfig->physical_res[0]; ii++)
			init.p_physical_set(jj,ii,table[0][jj]);
}

void Burgers_Plane_TS_ln_imex_mms::setup_look_up_table(double start, double end, double step_size)
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
		table_size = iiMax + 1;

		table = new double*[table_size];
		for (int ii=0; ii<table_size; ii++)
		{
			table[ii] = new double[simVars.disc.res_physical[0]];
			time = ii*step_size;
			if ((timestepping_order == 2) && (ii > 0))
				time = time-step_size/2;
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
void Burgers_Plane_TS_ln_imex_mms::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Burgers_Plane_TS_ln_imex_mms::Burgers_Plane_TS_ln_imex_mms(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Burgers_Plane_TS_ln_imex_mms::~Burgers_Plane_TS_ln_imex_mms()
{
	if (table != nullptr)
	{
		for (int ii=0; ii<table_size; ii++)
			delete[] table[ii];
		delete[] table;
		table = nullptr;
	}
}

