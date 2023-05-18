/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 *
 */

#include "SWE_Plane_Mori_Zwanzig_TS_n_erk.hpp"

#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


bool SWE_Plane_Mori_Zwanzig_TS_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	PDESWEPlaneMoriZwanzigTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_P(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);

}


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_bilinear(
		const sweet::PlaneData_Spectral &i_h_A,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_A,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_A,	///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_B,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_B,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_B,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u_A*ops->diff_c_x(i_u_B) - i_v_A*ops->diff_c_y(i_u_B);
	o_v_t = -i_u_A*ops->diff_c_x(i_v_B) - i_v_A*ops->diff_c_y(i_v_B);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u_B*ops->diff_c_x(i_h_A) + i_v_B*ops->diff_c_y(i_h_A));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u_B*i_h_A) - ops->diff_c_y(i_v_B*i_h_A);

}



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_Q(
		const sweet::PlaneData_Spectral &i_h_SQ,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_SQ,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_SQ,	///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_FQ,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_FQ,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_FQ,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_SQ_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_SQ_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_SQ_t,	///< time updates

		sweet::PlaneData_Spectral &o_h_FQ_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_FQ_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_FQ_t,	///< time updates

		double i_timestamp
)
{

		sweet::PlaneData_Spectral h_N_SF(i_h_SQ.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SF(i_u_SQ.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SF(i_v_SQ.planeDataConfig);
		sweet::PlaneData_Spectral h_N_SS(i_h_SQ.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SS(i_u_SQ.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SS(i_v_SQ.planeDataConfig);
		sweet::PlaneData_Spectral h_N_FF(i_h_SQ.planeDataConfig);
		sweet::PlaneData_Spectral u_N_FF(i_u_SQ.planeDataConfig);
		sweet::PlaneData_Spectral v_N_FF(i_v_SQ.planeDataConfig);
		sweet::PlaneData_Spectral h_N_FF_2(i_h_SQ.planeDataConfig);
		sweet::PlaneData_Spectral u_N_FF_2(i_u_SQ.planeDataConfig);
		sweet::PlaneData_Spectral v_N_FF_2(i_v_SQ.planeDataConfig);

		// compute nonlinear terms
		this->euler_timestep_update_bilinear(
							i_h_SQ, i_u_SQ, i_v_SQ,
							i_h_FQ, i_u_FQ, i_v_FQ,
							h_N_SF, u_N_SF, v_N_SF,
							i_timestamp
					);

		this->euler_timestep_update_bilinear(
							i_h_SQ, i_u_SQ, i_v_SQ,
							i_h_SQ, i_u_SQ, i_v_SQ,
							h_N_SS, u_N_SS, v_N_SS,
							i_timestamp
					);

		this->euler_timestep_update_bilinear(
							i_h_FQ, i_u_FQ, i_v_FQ,
							i_h_FQ, i_u_FQ, i_v_FQ,
							h_N_FF, u_N_FF, v_N_FF,
							i_timestamp
					);

		h_N_FF_2 = h_N_FF;
		u_N_FF_2 = u_N_FF;
		v_N_FF_2 = v_N_FF;


		// project
		this->projection.project_S(h_N_SF, u_N_SF, v_N_SF);
		this->projection.project_S(h_N_FF, u_N_FF, v_N_FF);
		this->projection.project_F(h_N_SS, u_N_SS, v_N_SS);
		this->projection.project_F(h_N_FF_2, u_N_FF_2, v_N_FF_2);

		// time update
		o_h_SQ_t = h_N_SF + h_N_FF;
		o_u_SQ_t = u_N_SF + u_N_FF;
		o_v_SQ_t = v_N_SF + v_N_FF;

		o_h_FQ_t = h_N_SS + h_N_FF_2;
		o_u_FQ_t = u_N_SS + u_N_FF_2;
		o_v_FQ_t = v_N_SS + v_N_FF_2;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_SF(
		const sweet::PlaneData_Spectral &i_h_S,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_S,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_S,	///< prognostic variables

		const sweet::PlaneData_Spectral &i_h_F,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u_F,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v_F,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_S_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_S_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_S_t,	///< time updates

		sweet::PlaneData_Spectral &o_h_F_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_F_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_F_t,	///< time updates

		double i_timestamp
)
{

		sweet::PlaneData_Spectral h_N_SF(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SF(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SF(i_v_S.planeDataConfig);
		sweet::PlaneData_Spectral h_N_SS(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SS(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SS(i_v_S.planeDataConfig);
		sweet::PlaneData_Spectral h_N_FF(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_FF(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_FF(i_v_S.planeDataConfig);
		sweet::PlaneData_Spectral h_N_SF_2(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SF_2(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SF_2(i_v_S.planeDataConfig);
		sweet::PlaneData_Spectral h_N_SS_2(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_SS_2(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_SS_2(i_v_S.planeDataConfig);
		sweet::PlaneData_Spectral h_N_FF_2(i_h_S.planeDataConfig);
		sweet::PlaneData_Spectral u_N_FF_2(i_u_S.planeDataConfig);
		sweet::PlaneData_Spectral v_N_FF_2(i_v_S.planeDataConfig);

		// compute nonlinear terms
		this->euler_timestep_update_bilinear(
							i_h_S, i_u_S, i_v_S,
							i_h_S, i_u_S, i_v_S,
							h_N_SS, u_N_SS, v_N_SS,
							i_timestamp
					);

		this->euler_timestep_update_bilinear(
							i_h_S, i_u_S, i_v_S,
							i_h_F, i_u_F, i_v_F,
							h_N_SF, u_N_SF, v_N_SF,
							i_timestamp
					);

		this->euler_timestep_update_bilinear(
							i_h_F, i_u_F, i_v_F,
							i_h_F, i_u_F, i_v_F,
							h_N_FF, u_N_FF, v_N_FF,
							i_timestamp
					);

		h_N_SS_2 = h_N_SS;
		u_N_SS_2 = u_N_SS;
		v_N_SS_2 = v_N_SS;
		h_N_SF_2 = h_N_SF;
		u_N_SF_2 = u_N_SF;
		v_N_SF_2 = v_N_SF;
		h_N_FF_2 = h_N_FF;
		u_N_FF_2 = u_N_FF;
		v_N_FF_2 = v_N_FF;


		// project
		this->projection.project_S(h_N_SS, u_N_SS, v_N_SS);
		this->projection.project_S(h_N_SF, u_N_SF, v_N_SF);
		this->projection.project_S(h_N_FF, u_N_FF, v_N_FF);
		this->projection.project_F(h_N_SS_2, u_N_SS_2, v_N_SS_2);
		this->projection.project_F(h_N_SF_2, u_N_SF_2, v_N_SF_2);
		this->projection.project_F(h_N_FF_2, u_N_FF_2, v_N_FF_2);

		// time update
		o_h_S_t = h_N_SS + h_N_SF + h_N_FF;
		o_u_S_t = u_N_SS + u_N_SF + u_N_FF;
		o_v_S_t = v_N_SS + v_N_SF + v_N_FF;

		o_h_F_t = h_N_SS_2 + h_N_SF_2 + h_N_FF_2;
		o_u_F_t = u_N_SS_2 + u_N_SF_2 + u_N_FF_2;
		o_v_F_t = v_N_SS_2 + v_N_SF_2 + v_N_FF_2;
}

void SWE_Plane_Mori_Zwanzig_TS_n_erk::runTimestep_P(
		sweet::PlaneData_Spectral &io_h_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_n_erk: Only constant time step size allowed");

	///////////////////////////
	// solve equation for SP //
	///////////////////////////

	sweet::PlaneData_Spectral h_N(io_h_SP.planeDataConfig);
	sweet::PlaneData_Spectral u_N(io_u_SP.planeDataConfig);
	sweet::PlaneData_Spectral v_N(io_v_SP.planeDataConfig);

	h_N = io_h_SP;
	u_N = io_u_SP;
	v_N = io_v_SP;

	timestepping_rk_P.runTimestep(
			this,
			&SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_P,	///< pointer to function to compute euler time step updates
			h_N, u_N, v_N,
			i_dt,
			timestepping_order_nonlinear_P,
			i_simulation_timestamp
		);

	this->projection.project_S(h_N, u_N, v_N);

	io_h_SP += h_N;
	io_u_SP += u_N;
	io_v_SP += v_N;
}


void SWE_Plane_Mori_Zwanzig_TS_n_erk::runTimestep_Q(
		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables


		double i_dt,
		double i_simulation_timestamp
)
{
	this->runTimestep_Q_SF(
					&SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_Q,

					io_h_pert_SQ,
					io_u_SQ,
					io_v_SQ,

					io_h_pert_FQ,
					io_u_FQ,
					io_v_FQ,

					i_dt,
					i_simulation_timestamp
		);
}

void SWE_Plane_Mori_Zwanzig_TS_n_erk::runTimestep_SF(
		sweet::PlaneData_Spectral &io_h_pert_S,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_S,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_S,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_F,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_F,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_F,	///< prognostic variables


		double i_dt,
		double i_simulation_timestamp
)
{
	this->runTimestep_Q_SF(
					&SWE_Plane_Mori_Zwanzig_TS_n_erk::euler_timestep_update_nonlinear_Q,

					io_h_pert_S,
					io_u_S,
					io_v_S,

					io_h_pert_F,
					io_u_F,
					io_v_F,

					i_dt,
					i_simulation_timestamp
		);
}


void SWE_Plane_Mori_Zwanzig_TS_n_erk::runTimestep_Q_SF(
		void (SWE_Plane_Mori_Zwanzig_TS_n_erk::*i_compute_euler_timestep_update)(
				const sweet::PlaneData_Spectral &i_h_A,
				const sweet::PlaneData_Spectral &i_u_A,
				const sweet::PlaneData_Spectral &i_v_A,

				const sweet::PlaneData_Spectral &i_h_B,
				const sweet::PlaneData_Spectral &i_u_B,
				const sweet::PlaneData_Spectral &i_v_B,

				sweet::PlaneData_Spectral &o_h_A,
				sweet::PlaneData_Spectral &o_u_A,
				sweet::PlaneData_Spectral &o_v_A,

				sweet::PlaneData_Spectral &o_h_B,
				sweet::PlaneData_Spectral &o_u_B,
				sweet::PlaneData_Spectral &o_v_B,

				double i_dt
		),

		sweet::PlaneData_Spectral &io_h_pert_A,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_A,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_A,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_B,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_B,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_B,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	/////////////////////////////////////////////////
	// solve equation for SQ and FQ                //
	// implement RK for system U^Q = [U^Q_S,U^Q_F] //
	/////////////////////////////////////////////////
	if (timestepping_order_nonlinear_Q == 1)
	{

		//this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A,
								io_u_A,
								io_v_A,
								io_h_pert_B,
								io_u_B,
								io_v_B,
								*RK_h_A_t[0],
								*RK_u_A_t[0],
								*RK_v_A_t[0],
								*RK_h_B_t[0],
								*RK_u_B_t[0],
								*RK_v_B_t[0],
								i_simulation_timestamp
				);

		io_h_pert_A += i_dt**RK_h_A_t[0];
		io_u_A      += i_dt**RK_u_A_t[0];
		io_v_A      += i_dt**RK_v_A_t[0];
		io_h_pert_B += i_dt**RK_h_B_t[0];
		io_u_B      += i_dt**RK_u_B_t[0];
		io_v_B      += i_dt**RK_v_B_t[0];

	}
	else if (timestepping_order_nonlinear_Q == 2)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/2 | 1/2
		 * --------------
		 *     | 0   1    b
		 */
		double a2[1] = {0.5};
		double b[2] = {0.0, 1.0};
		double c[1] = {0.5};

		// STAGE 1
		//this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A,
								io_u_A,
								io_v_A,
								io_h_pert_B,
								io_u_B,
								io_v_B,
								*RK_h_A_t[0],
								*RK_u_A_t[0],
								*RK_v_A_t[0],
								*RK_h_B_t[0],
								*RK_u_B_t[0],
								*RK_v_B_t[0],
								i_simulation_timestamp
				);

		// STAGE 2
		///this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A + ( i_dt * a2[0] * (*RK_h_A_t[0]) ),
								io_u_A      + ( i_dt * a2[0] * (*RK_u_A_t[0]) ),
								io_v_A      + ( i_dt * a2[0] * (*RK_v_A_t[0]) ),
								io_h_pert_B + ( i_dt * a2[0] * (*RK_h_B_t[0]) ),
								io_u_B      + ( i_dt * a2[0] * (*RK_u_B_t[0]) ),
								io_v_B      + ( i_dt * a2[0] * (*RK_v_B_t[0]) ),
								*RK_h_A_t[1],
								*RK_u_A_t[1],
								*RK_v_A_t[1],
								*RK_h_B_t[1],
								*RK_u_B_t[1],
								*RK_v_B_t[1],
								i_simulation_timestamp + c[0] * i_dt
				);

		io_h_pert_A += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_h_A_t[1]) );
		io_u_A      += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_u_A_t[1]) );
		io_v_A      += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_v_A_t[1]) );
		io_h_pert_B += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_h_B_t[1]) );
		io_u_B      += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_u_B_t[1]) );
		io_v_B      += i_dt * (/* b[0]*(*RK_h_t[0]) +*/ b[1] * (*RK_v_B_t[1]) );

	}
	else if (timestepping_order_nonlinear_Q == 3)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/3 | 1/3
		 * 2/3 | 0    2/3
		 * --------------
		 *     | 1/4  0   3/4
		 */
		double a2[1] = {1.0/3.0};
		double a3[2] = {0.0, 2.0/3.0};
		double b[3] = {1.0/4.0, 0.0, 3.0/4.0};
		double c[2] = {1.0/3.0, 2.0/3.0};

		// STAGE 1
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A,
								io_u_A,
								io_v_A,
								io_h_pert_B,
								io_u_B,
								io_v_B,
								*RK_h_A_t[0],
								*RK_u_A_t[0],
								*RK_v_A_t[0],
								*RK_h_B_t[0],
								*RK_u_B_t[0],
								*RK_v_B_t[0],
								i_simulation_timestamp
				);

		// STAGE 2
		///this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A + ( i_dt * a2[0] * (*RK_h_A_t[0]) ),
								io_u_A      + ( i_dt * a2[0] * (*RK_u_A_t[0]) ),
								io_v_A      + ( i_dt * a2[0] * (*RK_v_A_t[0]) ),
								io_h_pert_B + ( i_dt * a2[0] * (*RK_h_B_t[0]) ),
								io_u_B      + ( i_dt * a2[0] * (*RK_u_B_t[0]) ),
								io_v_B      + ( i_dt * a2[0] * (*RK_v_B_t[0]) ),
								*RK_h_A_t[1],
								*RK_u_A_t[1],
								*RK_v_A_t[1],
								*RK_h_B_t[1],
								*RK_u_B_t[1],
								*RK_v_B_t[1],
								i_simulation_timestamp + c[0] * i_dt
				);

		// STAGE 3
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A +  i_dt*( a3[0]*(*RK_h_A_t[0]) + a3[1]*(*RK_h_A_t[1]) ),
								io_u_A      +  i_dt*( a3[0]*(*RK_u_A_t[0]) + a3[1]*(*RK_u_A_t[1]) ),
								io_v_A      +  i_dt*( a3[0]*(*RK_v_A_t[0]) + a3[1]*(*RK_v_A_t[1]) ),
								io_h_pert_B +  i_dt*( a3[0]*(*RK_h_B_t[0]) + a3[1]*(*RK_h_B_t[1]) ),
								io_u_B      +  i_dt*( a3[0]*(*RK_u_B_t[0]) + a3[1]*(*RK_u_B_t[1]) ),
								io_v_B      +  i_dt*( a3[0]*(*RK_v_B_t[0]) + a3[1]*(*RK_v_B_t[1]) ),
								*RK_h_A_t[2],
								*RK_u_A_t[2],
								*RK_v_A_t[2],
								*RK_h_B_t[2],
								*RK_u_B_t[2],
								*RK_v_B_t[2],
								i_simulation_timestamp + c[1] * i_dt
				);

		io_h_pert_A += i_dt * ( (b[0]*(*RK_h_A_t[0])) + (b[1]*(*RK_h_A_t[1]))  + (b[2]*(*RK_h_A_t[2])) );
		io_u_A      += i_dt * ( (b[0]*(*RK_u_A_t[0])) + (b[1]*(*RK_u_A_t[1]))  + (b[2]*(*RK_u_A_t[2])) );
		io_v_A      += i_dt * ( (b[0]*(*RK_v_A_t[0])) + (b[1]*(*RK_v_A_t[1]))  + (b[2]*(*RK_v_A_t[2])) );
		io_h_pert_B += i_dt * ( (b[0]*(*RK_h_B_t[0])) + (b[1]*(*RK_h_B_t[1]))  + (b[2]*(*RK_h_B_t[2])) );
		io_u_B      += i_dt * ( (b[0]*(*RK_u_B_t[0])) + (b[1]*(*RK_u_B_t[1]))  + (b[2]*(*RK_u_B_t[2])) );
		io_v_B      += i_dt * ( (b[0]*(*RK_v_B_t[0])) + (b[1]*(*RK_v_B_t[1]))  + (b[2]*(*RK_v_B_t[2])) );

	}
	else if (timestepping_order_nonlinear_Q == 4)
	{
		// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
		// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
		/*
		 * c     a
		 * 0   |
		 * 1/2 | 1/2
		 * 1/2 | 0    1/2
		 * 1   | 0    0    1
		 * --------------
		 *     | 1/6  1/3  1/3  1/6
		 */
		double a2[1] = {0.5};
		double a3[2] = {0.0, 0.5};
		double a4[3] = {0.0, 0.0, 1.0};
		double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
		double c[3] = {0.5, 0.5, 1.0};

		// STAGE 1
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A,
								io_u_A,
								io_v_A,
								io_h_pert_B,
								io_u_B,
								io_v_B,
								*RK_h_A_t[0],
								*RK_u_A_t[0],
								*RK_v_A_t[0],
								*RK_h_B_t[0],
								*RK_u_B_t[0],
								*RK_v_B_t[0],
								i_simulation_timestamp
				);

		// STAGE 2
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A + ( i_dt * a2[0] * (*RK_h_A_t[0]) ),
								io_u_A      + ( i_dt * a2[0] * (*RK_u_A_t[0]) ),
								io_v_A      + ( i_dt * a2[0] * (*RK_v_A_t[0]) ),
								io_h_pert_B + ( i_dt * a2[0] * (*RK_h_B_t[0]) ),
								io_u_B      + ( i_dt * a2[0] * (*RK_u_B_t[0]) ),
								io_v_B      + ( i_dt * a2[0] * (*RK_v_B_t[0]) ),
								*RK_h_A_t[1],
								*RK_u_A_t[1],
								*RK_v_A_t[1],
								*RK_h_B_t[1],
								*RK_u_B_t[1],
								*RK_v_B_t[1],
								i_simulation_timestamp + c[0] * i_dt
				);

		// STAGE 3
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A +  i_dt*( /*a3[0]*(*RK_h_A_t[0]) +*/ a3[1]*(*RK_h_A_t[1]) ),
								io_u_A      +  i_dt*( /*a3[0]*(*RK_u_A_t[0]) +*/ a3[1]*(*RK_u_A_t[1]) ),
								io_v_A      +  i_dt*( /*a3[0]*(*RK_v_A_t[0]) +*/ a3[1]*(*RK_v_A_t[1]) ),
								io_h_pert_B +  i_dt*( /*a3[0]*(*RK_h_B_t[0]) +*/ a3[1]*(*RK_h_B_t[1]) ),
								io_u_B      +  i_dt*( /*a3[0]*(*RK_u_B_t[0]) +*/ a3[1]*(*RK_u_B_t[1]) ),
								io_v_B      +  i_dt*( /*a3[0]*(*RK_v_B_t[0]) +*/ a3[1]*(*RK_v_B_t[1]) ),
								*RK_h_A_t[2],
								*RK_u_A_t[2],
								*RK_v_A_t[2],
								*RK_h_B_t[2],
								*RK_u_B_t[2],
								*RK_v_B_t[2],
								i_simulation_timestamp + c[1] * i_dt
				);

		// STAGE 4
		////this->euler_timestep_update_nonlinear_Q(
		(this->*i_compute_euler_timestep_update)(
								io_h_pert_A +  i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_A_t[2]) ),
								io_u_A      +  i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_A_t[2]) ),
								io_v_A      +  i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_A_t[2]) ),
								io_h_pert_B +  i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_B_t[2]) ),
								io_u_B      +  i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_B_t[2]) ),
								io_v_B      +  i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_B_t[2]) ),
								*RK_h_A_t[3],
								*RK_u_A_t[3],
								*RK_v_A_t[3],
								*RK_h_B_t[3],
								*RK_u_B_t[3],
								*RK_v_B_t[3],
								i_simulation_timestamp + c[2] * i_dt
				);

		io_h_pert_A += i_dt * ( (b[0]*(*RK_h_A_t[0])) + (b[1]*(*RK_h_A_t[1]))  + (b[2]*(*RK_h_A_t[2]))  + (b[3]*(*RK_h_A_t[3])));
		io_u_A      += i_dt * ( (b[0]*(*RK_u_A_t[0])) + (b[1]*(*RK_u_A_t[1]))  + (b[2]*(*RK_u_A_t[2]))  + (b[3]*(*RK_u_A_t[3])));
		io_v_A      += i_dt * ( (b[0]*(*RK_v_A_t[0])) + (b[1]*(*RK_v_A_t[1]))  + (b[2]*(*RK_v_A_t[2]))  + (b[3]*(*RK_v_A_t[3])));
		io_h_pert_B += i_dt * ( (b[0]*(*RK_h_B_t[0])) + (b[1]*(*RK_h_B_t[1]))  + (b[2]*(*RK_h_B_t[2]))  + (b[3]*(*RK_h_B_t[3])));
		io_u_B      += i_dt * ( (b[0]*(*RK_u_B_t[0])) + (b[1]*(*RK_u_B_t[1]))  + (b[2]*(*RK_u_B_t[2]))  + (b[3]*(*RK_u_B_t[3])));
		io_v_B      += i_dt * ( (b[0]*(*RK_v_B_t[0])) + (b[1]*(*RK_v_B_t[1]))  + (b[2]*(*RK_v_B_t[2]))  + (b[3]*(*RK_v_B_t[3])));

	}
	else
		SWEETError("SWE_Plane_TS_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");

}




void SWE_Plane_Mori_Zwanzig_TS_n_erk::runTimestep(
		sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

		sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables


		double i_dt,
		double i_simulation_timestamp
)
{

	SWEETError("Don't use this function! Use instead run_Timestep_[P,Q,SF]");

	////if (i_dt <= 0)
	////	SWEETError("SWE_Plane_TS_n_erk: Only constant time step size allowed");

	////this->runTimestep_P(
	////			io_h_pert_SP,
	////			io_u_SP,
	////			io_v_SP,

	////			i_dt,
	////			i_simulation_timestamp
	////	);

	////this->runTimestep_Q(
	////			io_h_pert_SQ,
	////			io_u_SQ,
	////			io_v_SQ,

	////			io_h_pert_FQ,
	////			io_u_FQ,
	////			io_v_FQ,

	////			i_dt,
	////			i_simulation_timestamp
	////	);

}




	void SWE_Plane_Mori_Zwanzig_TS_n_erk::setupBuffers(
			const sweet::PlaneData_Config *i_planeDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
	)
	{
		if (RK_h_A_t != nullptr)	///< already allocated?
			return;

		int N = i_rk_order;

		if (N <= 0 || N > 4)
			SWEETError("Invalid order for RK time stepping (Please set --timestepping-order and/or --timestepping-order2)");

		RK_h_A_t = new sweet::PlaneData_Spectral*[N];
		RK_u_A_t = new sweet::PlaneData_Spectral*[N];
		RK_v_A_t = new sweet::PlaneData_Spectral*[N];
		RK_h_B_t = new sweet::PlaneData_Spectral*[N];
		RK_u_B_t = new sweet::PlaneData_Spectral*[N];
		RK_v_B_t = new sweet::PlaneData_Spectral*[N];

		for (int i = 0; i < N; i++)
		{
			RK_h_A_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
			RK_u_A_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
			RK_v_A_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
			RK_h_B_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
			RK_u_B_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
			RK_v_B_t[i] = new sweet::PlaneData_Spectral(i_planeDataConfig);
		}
	}


/*
 * Setup
 */
bool SWE_Plane_Mori_Zwanzig_TS_n_erk::setup(
	sweet::PlaneOperators *io_ops,
	const sweet::PlaneData_Config *i_planeDataConfig,
	std::string i_equation
)
{

	PDESWEPlaneMoriZwanzigTS_BaseInterface::setup(io_ops, i_equation);

	ops = io_ops;

	equation = i_equation;

	RK_h_A_t = nullptr;
	RK_u_A_t = nullptr;
	RK_v_A_t = nullptr;
	RK_h_B_t = nullptr;
	RK_u_B_t = nullptr;
	RK_v_B_t = nullptr;


	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	timestepping_order_nonlinear_P = shackPDESWETimeDisc->timestepping_order_P;
	timestepping_order_nonlinear_Q = shackPDESWETimeDisc->timestepping_order_Q;

	timestepping_rk_P.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear_P);

	this->setupBuffers(i_planeDataConfig, timestepping_order_nonlinear_Q);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for n_erk");

	return true;
}

SWE_Plane_Mori_Zwanzig_TS_n_erk::~SWE_Plane_Mori_Zwanzig_TS_n_erk()
{
	int N = this->timestepping_order_nonlinear_Q;

	if (RK_h_A_t != nullptr)
	{
		for (int i = 0; i < N; i++)
		{
			delete RK_h_A_t[i];
			delete RK_u_A_t[i];
			delete RK_v_A_t[i];
			delete RK_h_B_t[i];
			delete RK_u_B_t[i];
			delete RK_v_B_t[i];
		}

		delete [] RK_h_A_t;
		delete [] RK_u_A_t;
		delete [] RK_v_A_t;
		delete [] RK_h_B_t;
		delete [] RK_u_B_t;
		delete [] RK_v_B_t;

		RK_h_A_t = nullptr;
		RK_u_A_t = nullptr;
		RK_v_A_t = nullptr;
		RK_h_B_t = nullptr;
		RK_u_B_t = nullptr;
		RK_v_B_t = nullptr;
	}
}

