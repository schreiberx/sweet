/*
 * SWE_Plane_TS_l_direct.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "SWE_Plane_TS_l_direct.hpp"
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/Staggering.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>



void SWE_Plane_TS_l_direct::run_timestep(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (simVars.disc.use_staggering)
		run_timestep_cgrid(io_h_pert, io_u, io_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);
	else
		run_timestep_agrid(io_h_pert, io_u, io_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);
}




/**
 * Computation of analytical solution on staggered grid
 */
void SWE_Plane_TS_l_direct::run_timestep_cgrid(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,		///< prognostic variables
		PlaneData &io_v,		///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	// For output, variables need to be on unstaggered A-grid
	PlaneData t_u(io_h_pert.planeDataConfig);
	PlaneData t_v(io_h_pert.planeDataConfig);

	if (!simVars.disc.use_staggering)
		FatalError("Expected staggering");

	planeDataGridMapping.mapCtoA_u(io_u, t_u);
	planeDataGridMapping.mapCtoA_v(io_v, t_v);

	simVars.disc.use_staggering = false;

	run_timestep_agrid(
			io_h_pert, t_u, t_v,
			o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time
	);

	simVars.disc.use_staggering = true;

	planeDataGridMapping.mapAtoC_u(t_u, io_u);
	planeDataGridMapping.mapAtoC_v(t_v, io_v);
}



void SWE_Plane_TS_l_direct::run_timestep_agrid(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	run_timestep_agrid_planedata(io_h_pert, io_u, io_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);
#else
	run_timestep_agrid_planedatacomplex(io_h_pert, io_u, io_v, o_dt, i_fixed_dt, i_simulation_timestamp, i_max_simulation_time);
#endif
}


#if SWEET_USE_PLANE_SPECTRAL_SPACE

/**
 * This method computes the analytical solution based on the given initial values.
 *
 * See Embid/Madja/1996, Terry/Beth/2014, page 16
 * and
 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
 * for the dimensionful formulation.
 */
void SWE_Plane_TS_l_direct::run_timestep_agrid_planedata(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (simVars.disc.use_staggering)
		FatalError("Staggering not supported");

	if (i_fixed_dt < 0)
		FatalError("SWE_Plane_TS_l_direct: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;

	typedef std::complex<double> complex;

	complex I(0.0, 1.0);


	/*
	 * This implementation works directly on PlaneData
	 */
	double s0 = simVars.sim.domain_size[0];
	double s1 = simVars.sim.domain_size[1];

	io_h_pert.request_data_spectral();
	io_u.request_data_spectral();
	io_v.request_data_spectral();

	double f = simVars.sim.f0;
	double h = simVars.sim.h0;
	double g = simVars.sim.gravitation;

	double sqrt_h = std::sqrt(h);
	double sqrt_g = std::sqrt(g);


	for (std::size_t ik1 = 0; ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]; ik1++)
	{
		double k1;
		if (ik1 < io_h_pert.planeDataConfig->spectral_data_size[1]/2)
			k1 = (double)ik1;
		else
			k1 = (double)((int)ik1-(int)io_h_pert.planeDataConfig->spectral_data_size[1]);

		for (std::size_t ik0 = 0; ik0 < io_h_pert.planeDataConfig->spectral_data_size[0]; ik0++)
		{
			double k0 = (double)ik0;

			complex U[3];
			U[0] = io_h_pert.spectral_get(ik1, ik0);
			U[1] = io_u.spectral_get(ik1, ik0);
			U[2] = io_v.spectral_get(ik1, ik0);

			complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
			complex c = -k1*I;

			b = b*2.0*M_PI/s0;
			c = c*2.0*M_PI/s1;

			/*
			 * Matrix with Eigenvectors (column-wise)
			 */
			complex v[3][3];

			/*
			 * Eigenvalues
			 */
			complex lambda[3];

			if (simVars.sim.f0 == 0)
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

					v[0][1] = -(sqrt_h*std::sqrt(b*b + c*c))/(c*sqrt_g);
					v[1][1] = b/c;
					v[2][1] = 1.0;

					v[0][2] = (sqrt_h*std::sqrt(b*b + c*c))/(c*sqrt_g);
					v[1][2] = b/c;
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -std::sqrt(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = std::sqrt(b*b + c*c)*sqrt_h*sqrt_g;
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

					v[0][1] = -(c*h)/std::sqrt(-f*f + c*c*g*h);
					v[1][1] =  -f/std::sqrt(-f*f + c*c*g*h);
					v[2][1] = 1;

					v[0][2] = (c*h)/std::sqrt(-f*f + c*c*g*h);
					v[1][2] = f/std::sqrt(-f*f + c*c*g*h);
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -std::sqrt(c*c*g*h-f*f);
					lambda[2] = std::sqrt(c*c*g*h-f*f);
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
					v[1][1] = std::sqrt(-f*f + b*b*g*h)/f;
					v[2][1] = 1;

					v[0][2] = -(b*h)/f;
					v[1][2] = -std::sqrt(-f*f + b*b*g*h)/f;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -std::sqrt(b*b*g*h-f*f);
					lambda[2] = std::sqrt(b*b*g*h-f*f);
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

					v[0][1] = -(c*f*h + b*h*std::sqrt(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*std::sqrt(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -std::sqrt(b*b*g*h + c*c*g*h - f*f);
					lambda[2] =  std::sqrt(b*b*g*h + c*c*g*h - f*f);
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
				UEV[k] = std::exp(i_fixed_dt*lambda[k])*UEV[k];

			for (int k = 0; k < 3; k++)
				U[k] = 0.0;

			for (int k = 0; k < 3; k++)
				for (int j = 0; j < 3; j++)
					U[k] += v[k][j] * UEV[j];

			io_h_pert.p_spectral_set(ik1, ik0, U[0]);
			io_u.p_spectral_set(ik1, ik0, U[1]);
			io_v.p_spectral_set(ik1, ik0, U[2]);
		}
	}

	io_h_pert.spectral_zeroAliasingModes();
	io_u.spectral_zeroAliasingModes();
	io_v.spectral_zeroAliasingModes();

	o_dt = i_fixed_dt;
}

#endif


void SWE_Plane_TS_l_direct::run_timestep_agrid_planedatacomplex(
		PlaneData &io_h_pert,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt < 0)
		FatalError("SWE_Plane_TS_l_direct: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time - i_simulation_timestamp;

	typedef std::complex<double> complex;

	complex I(0.0, 1.0);


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
#warning "WARNING: Not doing this in spectral space leads to a loss of the highest mode"
	/*
	 * WARNING: This leads to a loss of precision due to the highest mode which cannot be
	 * tracked due to the Nyquist theorem
	 */
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

	double s0 = simVars.sim.domain_size[0];
	double s1 = simVars.sim.domain_size[1];

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	o_h_pert.spectral_space_data_valid = true;
	o_h_pert.physical_space_data_valid = false;

	o_u.spectral_space_data_valid = true;
	o_u.physical_space_data_valid = false;

	o_v.spectral_space_data_valid = true;
	o_v.physical_space_data_valid = false;
#endif

	double f = simVars.sim.f0;
	double h = simVars.sim.h0;
	double g = simVars.sim.gravitation;

	double sqrt_h = std::sqrt(h);
	double sqrt_g = std::sqrt(g);

	for (std::size_t ik1 = 0; ik1 < i_h_pert.planeDataConfig->spectral_complex_data_size[1]; ik1++)
	{
		double k1;
		if (ik1 < i_h_pert.planeDataConfig->spectral_complex_data_size[1]/2)
			k1 = (double)ik1;
		else
			k1 = -(double)((int)ik1-(int)i_h_pert.planeDataConfig->spectral_complex_data_size[1]);

		for (std::size_t ik0 = 0; ik0 < i_h_pert.planeDataConfig->spectral_complex_data_size[0]; ik0++)
		{
			double k0;
			if (ik0 < i_h_pert.planeDataConfig->spectral_complex_data_size[0]/2)
				k0 = (double)ik0;
			else
				k0 = (double)((int)ik0-(int)i_h_pert.planeDataConfig->spectral_complex_data_size[0]);

			complex U[3];
			U[0] = i_h_pert.p_spectral_get(ik1, ik0);
			U[1] = i_u.p_spectral_get(ik1, ik0);
			U[2] = i_v.p_spectral_get(ik1, ik0);

			complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
			complex c = -k1*I;

			b = b*2.0*M_PI/s0;
			c = c*2.0*M_PI/s1;

			/*
			 * Matrix with Eigenvectors (column-wise)
			 */
			complex v[3][3];

			/*
			 * Eigenvalues
			 */
			complex lambda[3];

			if (simVars.sim.f0 == 0)
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

					v[0][1] = -(sqrt_h*std::sqrt(b*b + c*c))/(c*sqrt_g);
					v[1][1] = b/c;
					v[2][1] = 1.0;

					v[0][2] = (sqrt_h*std::sqrt(b*b + c*c))/(c*sqrt_g);
					v[1][2] = b/c;
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -std::sqrt(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = std::sqrt(b*b + c*c)*sqrt_h*sqrt_g;
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

					v[0][1] = -(c*h)/std::sqrt(-f*f + c*c*g*h);
					v[1][1] =  -f/std::sqrt(-f*f + c*c*g*h);
					v[2][1] = 1;

					v[0][2] = (c*h)/std::sqrt(-f*f + c*c*g*h);
					v[1][2] = f/std::sqrt(-f*f + c*c*g*h);
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -std::sqrt(c*c*g*h-f*f);
					lambda[2] = std::sqrt(c*c*g*h-f*f);
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
					v[1][1] = std::sqrt(-f*f + b*b*g*h)/f;
					v[2][1] = 1;

					v[0][2] = -(b*h)/f;
					v[1][2] = -std::sqrt(-f*f + b*b*g*h)/f;
					v[2][2] = 1;

					lambda[0] = 0;
					lambda[1] = -std::sqrt(b*b*g*h-f*f);
					lambda[2] = std::sqrt(b*b*g*h-f*f);
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

					v[0][1] = -(c*f*h + b*h*std::sqrt(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*std::sqrt(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*std::sqrt(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					lambda[0] = 0.0;
					lambda[1] = -std::sqrt(b*b*g*h + c*c*g*h - f*f);
					lambda[2] =  std::sqrt(b*b*g*h + c*c*g*h - f*f);
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
				UEV[k] = std::exp(i_fixed_dt*lambda[k])*UEV[k];

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
			o_h_pert.p_spectral_set(ik1, ik0, U[0]);
			o_u.p_spectral_set(ik1, ik0, U[1]);
			o_v.p_spectral_set(ik1, ik0, U[2]);
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
#warning "WARNING: Not doing this in spectral space leads to a loss of the highest mode"
	/*
	 * WARNING: This leads to a loss of precision due to the highest mode which cannot be
	 * tracked due to the Nyquist theorem
	 */
	io_h_pert = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_h_pert);
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_v);
#else
	io_h_pert = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_h_pert);
	io_u = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_v);
#endif
	o_dt = i_fixed_dt;
}


SWE_Plane_TS_l_direct::SWE_Plane_TS_l_direct(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	if (simVars.disc.use_staggering)
		planeDataGridMapping.setup(i_simVars, op.planeDataConfig);
}



SWE_Plane_TS_l_direct::~SWE_Plane_TS_l_direct()
{
}

