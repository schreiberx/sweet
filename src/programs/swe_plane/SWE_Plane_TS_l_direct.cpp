/*
 * SWE_Plane_TS_l_direct.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../swe_plane/SWE_Plane_TS_l_direct.hpp"

#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>
#include <sweet/plane/PlaneStaggering.hpp>


void SWE_Plane_TS_l_direct::setup(
		const std::string &i_function_name
)
{
	rexiFunctions.setup(i_function_name);
}



void SWE_Plane_TS_l_direct::run_timestep(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.disc.space_grid_use_c_staggering)
		run_timestep_cgrid(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
	else
		run_timestep_agrid(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
}




/**
 * Computation of analytical solution on staggered grid
 */
void SWE_Plane_TS_l_direct::run_timestep_cgrid(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,		///< prognostic variables
		PlaneData &io_v,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	// For output, variables need to be on unstaggered A-grid
	PlaneData t_u(io_h_pert.planeDataConfig);
	PlaneData t_v(io_h_pert.planeDataConfig);

	if (!simVars.disc.space_grid_use_c_staggering)
		FatalError("Expected staggering");

	planeDataGridMapping.mapCtoA_u(io_u, t_u);
	planeDataGridMapping.mapCtoA_v(io_v, t_v);

	simVars.disc.space_grid_use_c_staggering = false;

	run_timestep_agrid(
			io_h_pert, t_u, t_v,
			i_dt, i_simulation_timestamp
	);

	simVars.disc.space_grid_use_c_staggering = true;

	planeDataGridMapping.mapAtoC_u(t_u, io_u);
	planeDataGridMapping.mapAtoC_v(t_v, io_v);
}



void SWE_Plane_TS_l_direct::run_timestep_agrid(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	run_timestep_agrid_planedata(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
#else
	run_timestep_agrid_planedatacomplex(io_h_pert, io_u, io_v, i_dt, i_simulation_timestamp);
#endif
}


#if SWEET_USE_PLANE_SPECTRAL_SPACE

/**
 * This method computes the analytical solution based on the given initial values.
 *
 * See Embid/Madja/1996, Terry/Beth/2014, page 16
 * and
 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
 * for the dimension full formulation.
 */
void SWE_Plane_TS_l_direct::run_timestep_agrid_planedata(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("Staggering not supported");

	//if (i_dt < 0)
	//	FatalError("SWE_Plane_TS_l_direct: Only constant time step size allowed (please set --dt )");


	typedef std::complex<T> complex;
	complex I(0.0, 1.0);

	T dt = i_dt;


	/*
	 * This implementation works directly on PlaneData
	 */
	T s0 = simVars.sim.plane_domain_size[0];
	T s1 = simVars.sim.plane_domain_size[1];

	io_h_pert.request_data_spectral();
	io_u.request_data_spectral();
	io_v.request_data_spectral();

	T f = simVars.sim.plane_rotating_f0;
	T h = simVars.sim.h0;
	T g = simVars.sim.gravitation;

	T sqrt_h = rexiFunctions.l_sqrt(h);
	T sqrt_g = rexiFunctions.l_sqrt(g);


#if SWEET_THREADING_SPACE
	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
#endif
	for (std::size_t ik1 = 0; ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < io_h_pert.planeDataConfig->spectral_data_size[0]; ik0++)
		{
			T k1;
			if (ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]/2)
				k1 = (T)ik1;
			else
				k1 = (T)((int)ik1-(int)io_h_pert.planeDataConfig->spectral_data_size[1]);

			T k0 = (T)ik0;

			complex U[3];
			U[0] = io_h_pert.spectral_get(ik1, ik0);
			U[1] = io_u.spectral_get(ik1, ik0);
			U[2] = io_v.spectral_get(ik1, ik0);

			complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
			complex c = -k1*I;

			b = b*rexiFunctions.pi2/s0;
			c = c*rexiFunctions.pi2/s1;

			/*
			 * Matrix with Eigenvectors (column-wise)
			 */
			complex v[3][3];

			/*
			 * Eigenvalues
			 */
			complex lambda[3];

			if (simVars.sim.plane_rotating_f0 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,0%7D,%7Bg*c,0,0%7D%7D
				 */
				if (k0 == 0 && k1 == 0)
				{
					v[0][0] = 1;
					v[1][0] = 0;
					v[2][0] = 0;

					v[0][1] = 0;
					v[1][1] = 1;
					v[2][1] = 0;

					v[0][2] = 0;
					v[1][2] = 0;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = 0;
					lambda[2] = 0;
				}
				else if (k0 == 0)
				{
					v[0][0] = 0;
					v[1][0] = 1;
					v[2][0] = 0;

					v[0][1] = -sqrt_h/sqrt_g;
					v[1][1] = 0;
					v[2][1] = 1;

					v[0][2] = sqrt_h/sqrt_g;
					v[1][2] = 0;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -c*sqrt_g*sqrt_h;
					lambda[2] = c*sqrt_g*sqrt_h;;
				}
				else if (k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,0%7D,%7Bg*c*0,0,0%7D%7D
					 */

					v[0][0] = 0;
					v[1][0] = 0;
					v[2][0] = 1;

					v[0][1] = -sqrt_h/sqrt_g;
					v[1][1] = 1;
					v[2][1] = 0;

					v[0][2] = sqrt_h/sqrt_g;
					v[1][2] = 1;
					v[2][2] = 0;

					lambda[0] = 0;
					lambda[1] = -b*sqrt_g*sqrt_h;
					lambda[2] = b*sqrt_g*sqrt_h;
				}
				else
				{
					v[0][0] = 0;
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
					v[1][1] = b/c;
					v[2][1] = 1.0;

					v[0][2] = (sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
					v[1][2] = b/c;
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
				}
			}
			else
			{
				if (k0 == 0 && k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
					 */
					v[0][0] = 0;
					v[1][0] = -I;
					v[2][0] = 1;

					v[0][1] = 0;
					v[1][1] = I;
					v[2][1] = 1;

					v[0][2] = 1;
					v[1][2] = 0;
					v[2][2] = 0;

					lambda[0] = I*f;
					lambda[1] = -I*f;
					lambda[2] = 0;
				}
				else if (k0 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b*0,h*c%7D,%7Bg*b*0,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */
					v[0][0] = f/(c*g);
					v[1][0] = 1;
					v[2][0] = 0;

					v[0][1] = -(c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[1][1] =  -f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[2][1] = 1;

					v[0][2] = (c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[1][2] = f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
				}
				else if (k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,f%7D,%7Bg*c*0,-f,0%7D%7D
					 */
					v[0][0] = -f/(b*g);
					v[1][0] = 0;
					v[2][0] = 1;

					v[0][1] = -(b*h)/f;
					v[1][1] = rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
					v[2][1] = 1;

					v[0][2] = -(b*h)/f;
					v[1][2] = -rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
				}
				else
				{
					/*
					 * Compute EV's of
					 * Linear operator
					 *
					 * [ 0  hb  hc ]
					 * [ gb  0   f ]
					 * [ gc -f   0 ]
					 *
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */

					v[0][0] = -f/(b*g);
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
					lambda[2] =  rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
				}
			}

			/*
			 * Invert Eigenvalue matrix
			 */
			complex v_inv[3][3];

			v_inv[0][0] =  (v[1][1]*v[2][2] - v[1][2]*v[2][1]);
			v_inv[0][1] = -(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
			v_inv[0][2] =  (v[0][1]*v[1][2] - v[0][2]*v[1][1]);

			v_inv[1][0] = -(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
			v_inv[1][1] =  (v[0][0]*v[2][2] - v[0][2]*v[2][0]);
			v_inv[1][2] = -(v[0][0]*v[1][2] - v[0][2]*v[1][0]);

			v_inv[2][0] =  (v[1][0]*v[2][1] - v[1][1]*v[2][0]);
			v_inv[2][1] = -(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
			v_inv[2][2] =  (v[0][0]*v[1][1] - v[0][1]*v[1][0]);

			complex s = v[0][0]*v_inv[0][0] + v[0][1]*v_inv[1][0] + v[0][2]*v_inv[2][0];

			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					v_inv[j][i] /= s;

			complex UEV[3] = {0.0, 0.0, 0.0};
			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					UEV[k] += v_inv[k][j] * U[j];


			for (int k = 0; k < 3; k++)
			{
				std::complex<T> &lam = lambda[k];

				std::complex<T> K = rexiFunctions.eval(lam*dt);

				UEV[k] = K*UEV[k];
			}

			for (int k = 0; k < 3; k++)
				U[k] = 0.0;

			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					U[k] += v[k][j] * UEV[j];


#if SWEET_QUADMATH
			std::complex<double> tmp0(U[0].real(), U[0].imag());
			io_h_pert.p_spectral_set(ik1, ik0, tmp0);

			std::complex<double> tmp1(U[1].real(), U[1].imag());
			io_u.p_spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[2].real(), U[2].imag());
			io_v.p_spectral_set(ik1, ik0, tmp2);
#else
			io_h_pert.p_spectral_set(ik1, ik0, U[0]);
			io_u.p_spectral_set(ik1, ik0, U[1]);
			io_v.p_spectral_set(ik1, ik0, U[2]);
#endif
		}
	}

	io_h_pert.spectral_zeroAliasingModes();
	io_u.spectral_zeroAliasingModes();
	io_v.spectral_zeroAliasingModes();
}

#endif


void SWE_Plane_TS_l_direct::run_timestep_agrid_planedatacomplex(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt < 0)
		FatalError("SWE_Plane_TS_l_direct: Only constant time step size allowed");

	typedef std::complex<T> complex;

	complex I(0.0, 1.0);


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	PlaneDataComplex i_h_pert = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_h_pert);
	PlaneDataComplex i_u = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_u);
	PlaneDataComplex i_v = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_v);
#else
	PlaneDataComplex i_h_pert = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_h_pert);
	PlaneDataComplex i_u = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_u);
	PlaneDataComplex i_v = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_v);
#endif

	PlaneDataComplex o_h_pert(io_h_pert.planeDataConfig);
	PlaneDataComplex o_u(io_h_pert.planeDataConfig);
	PlaneDataComplex o_v(io_h_pert.planeDataConfig);

	T dt = i_dt;

	T s0 = simVars.sim.plane_domain_size[0];
	T s1 = simVars.sim.plane_domain_size[1];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	o_h_pert.spectral_space_data_valid = true;
	o_h_pert.physical_space_data_valid = false;

	o_u.spectral_space_data_valid = true;
	o_u.physical_space_data_valid = false;

	o_v.spectral_space_data_valid = true;
	o_v.physical_space_data_valid = false;
#endif

	T f = simVars.sim.plane_rotating_f0;
	T h = simVars.sim.h0;
	T g = simVars.sim.gravitation;

	T sqrt_h = rexiFunctions.l_sqrt(h);
	T sqrt_g = rexiFunctions.l_sqrt(g);

	SWEET_THREADING_SPACE_PARALLEL_FOR_SIMD_COLLAPSE2
	for (std::size_t ik1 = 0; ik1 < i_h_pert.planeDataConfig->spectral_complex_data_size[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < i_h_pert.planeDataConfig->spectral_complex_data_size[0]; ik0++)
		{
			T k1;
			if (ik1 < i_h_pert.planeDataConfig->spectral_complex_data_size[1]/2)
				k1 = (T)ik1;
			else
				k1 = -(T)((int)ik1-(int)i_h_pert.planeDataConfig->spectral_complex_data_size[1]);

			T k0;
			if (ik0 < i_h_pert.planeDataConfig->spectral_complex_data_size[0]/2)
				k0 = (T)ik0;
			else
				k0 = (T)((int)ik0-(int)i_h_pert.planeDataConfig->spectral_complex_data_size[0]);

			complex U[3];
			U[0] = i_h_pert.p_spectral_get(ik1, ik0);
			U[1] = i_u.p_spectral_get(ik1, ik0);
			U[2] = i_v.p_spectral_get(ik1, ik0);

			complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
			complex c = -k1*I;

			b = b*rexiFunctions.pi2/s0;
			c = c*rexiFunctions.pi2/s1;

			/*
			 * Matrix with Eigenvectors (column-wise)
			 */
			complex v[3][3];

			/*
			 * Eigenvalues
			 */
			complex lambda[3];

			if (simVars.sim.plane_rotating_f0 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,0%7D,%7Bg*c,0,0%7D%7D
				 */
				if (k0 == 0 && k1 == 0)
				{
					v[0][0] = 1;
					v[1][0] = 0;
					v[2][0] = 0;

					v[0][1] = 0;
					v[1][1] = 1;
					v[2][1] = 0;

					v[0][2] = 0;
					v[1][2] = 0;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = 0;
					lambda[2] = 0;
				}
				else if (k0 == 0)
				{
					v[0][0] = 0;
					v[1][0] = 1;
					v[2][0] = 0;

					v[0][1] = -sqrt_h/sqrt_g;
					v[1][1] = 0;
					v[2][1] = 1;

					v[0][2] = sqrt_h/sqrt_g;
					v[1][2] = 0;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -c*sqrt_g*sqrt_h;
					lambda[2] = c*sqrt_g*sqrt_h;;
				}
				else if (k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,0%7D,%7Bg*c*0,0,0%7D%7D
					 */

					v[0][0] = 0;
					v[1][0] = 0;
					v[2][0] = 1;

					v[0][1] = -sqrt_h/sqrt_g;
					v[1][1] = 1;
					v[2][1] = 0;

					v[0][2] = sqrt_h/sqrt_g;
					v[1][2] = 1;
					v[2][2] = 0;

					lambda[0] = 0;
					lambda[1] = -b*sqrt_g*sqrt_h;
					lambda[2] = b*sqrt_g*sqrt_h;
				}
				else
				{
					v[0][0] = 0;
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
					v[1][1] = b/c;
					v[2][1] = 1.0;

					v[0][2] = (sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
					v[1][2] = b/c;
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
				}
			}
			else
			{
				if (k0 == 0 && k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
					 */
					v[0][0] = 0;
					v[1][0] = -I;
					v[2][0] = 1;

					v[0][1] = 0;
					v[1][1] = I;
					v[2][1] = 1;

					v[0][2] = 1;
					v[1][2] = 0;
					v[2][2] = 0;

					lambda[0] = I*f;
					lambda[1] = -I*f;
					lambda[2] = 0;
				}
				else if (k0 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b*0,h*c%7D,%7Bg*b*0,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */
					v[0][0] = f/(c*g);
					v[1][0] = 1;
					v[2][0] = 0;

					v[0][1] = -(c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[1][1] =  -f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[2][1] = 1;

					v[0][2] = (c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[1][2] = f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
				}
				else if (k1 == 0)
				{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,f%7D,%7Bg*c*0,-f,0%7D%7D
					 */
					v[0][0] = -f/(b*g);
					v[1][0] = 0;
					v[2][0] = 1;

					v[0][1] = -(b*h)/f;
					v[1][1] = rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
					v[2][1] = 1;

					v[0][2] = -(b*h)/f;
					v[1][2] = -rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
				}
				else
				{
					/*
					 * Compute EV's of
					 * Linear operator
					 *
					 * [ 0  hb  hc ]
					 * [ gb  0   f ]
					 * [ gc -f   0 ]
					 *
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */

					v[0][0] = -f/(b*g);
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
					lambda[2] =  rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
				}
			}

			/*
			 * Invert Eigenvalue matrix
			 */
			complex v_inv[3][3];

			v_inv[0][0] =  (v[1][1]*v[2][2] - v[1][2]*v[2][1]);
			v_inv[0][1] = -(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
			v_inv[0][2] =  (v[0][1]*v[1][2] - v[0][2]*v[1][1]);

			v_inv[1][0] = -(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
			v_inv[1][1] =  (v[0][0]*v[2][2] - v[0][2]*v[2][0]);
			v_inv[1][2] = -(v[0][0]*v[1][2] - v[0][2]*v[1][0]);

			v_inv[2][0] =  (v[1][0]*v[2][1] - v[1][1]*v[2][0]);
			v_inv[2][1] = -(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
			v_inv[2][2] =  (v[0][0]*v[1][1] - v[0][1]*v[1][0]);

			complex s = v[0][0]*v_inv[0][0] + v[0][1]*v_inv[1][0] + v[0][2]*v_inv[2][0];

			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					v_inv[j][i] /= s;

			complex UEV[3] = {0.0, 0.0, 0.0};
			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					UEV[k] += v_inv[k][j] * U[j];


			for (int k = 0; k < 3; k++)
			{
				std::complex<T> &lam = lambda[k];

				std::complex<T> K = rexiFunctions.eval(lam*dt);

				UEV[k] = K*UEV[k];
			}

			for (int k = 0; k < 3; k++)
				U[k] = 0.0;

			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					U[k] += v[k][j] * UEV[j];

			/*
			 * Make sure that symmetry is correct
			 */
			if (ik0 > i_h_pert.planeDataConfig->spectral_complex_data_size[0]/2)
			{
				for (int i = 0; i < 3; i++)
					U[i].imag(-U[i].imag());
			}

			/*
			 * Convert to EV space
			 */
#if SWEET_QUADMATH
			std::complex<double> tmp0(U[0].real(), U[0].imag());
			o_h_pert.p_spectral_set(ik1, ik0, tmp0);

			std::complex<double> tmp1(U[1].real(), U[1].imag());
			o_u.p_spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[2].real(), U[2].imag());
			o_v.p_spectral_set(ik1, ik0, tmp2);
#else
			o_h_pert.p_spectral_set(ik1, ik0, U[0]);
			o_u.p_spectral_set(ik1, ik0, U[1]);
			o_v.p_spectral_set(ik1, ik0, U[2]);
#endif
		}
	}

#if SWEET_DEBUG
	o_h_pert.test_realphysical();
	o_u.test_realphysical();
	o_v.test_realphysical();
#endif

	o_h_pert.spectral_zeroAliasingModes();
	o_u.spectral_zeroAliasingModes();
	o_v.spectral_zeroAliasingModes();

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	io_h_pert = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_h_pert);
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_v);
#else
	io_h_pert = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_h_pert);
	io_u = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_v);
#endif
}


SWE_Plane_TS_l_direct::SWE_Plane_TS_l_direct(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	if (simVars.disc.space_grid_use_c_staggering)
		planeDataGridMapping.setup(i_simVars, op.planeDataConfig);
}



SWE_Plane_TS_l_direct::~SWE_Plane_TS_l_direct()
{
}

