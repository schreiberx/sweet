/*
 * SWE_Plane_TS_l_direct.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "Burgers_Plane_TS_l_direct.hpp"
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/Staggering.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>



void Burgers_Plane_TS_l_direct:: setup()
{
}



void Burgers_Plane_TS_l_direct::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_prev_u,	///< prognostic variables
		PlaneData &io_prev_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.disc.use_staggering)
		run_timestep_cgrid(io_u, io_v, i_dt, i_simulation_timestamp);
	else
		run_timestep_agrid(io_u, io_v, i_dt, i_simulation_timestamp);
}




/**
 * Computation of analytical solution on staggered grid
 */
void Burgers_Plane_TS_l_direct::run_timestep_cgrid(
		PlaneData &io_u,		///< prognostic variables
		PlaneData &io_v,		///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	// For output, variables need to be on unstaggered A-grid
	PlaneData t_u(io_u.planeDataConfig);
	PlaneData t_v(io_u.planeDataConfig);

	if (!simVars.disc.use_staggering)
		FatalError("Expected staggering");

	planeDataGridMapping.mapCtoA_u(io_u, t_u);
	planeDataGridMapping.mapCtoA_v(io_v, t_v);

	simVars.disc.use_staggering = false;

	run_timestep_agrid(
			t_u, t_v,
			i_dt, i_simulation_timestamp
	);

	simVars.disc.use_staggering = true;

	planeDataGridMapping.mapAtoC_u(t_u, io_u);
	planeDataGridMapping.mapAtoC_v(t_v, io_v);
}



void Burgers_Plane_TS_l_direct::run_timestep_agrid(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	run_timestep_agrid_planedata(io_u, io_v, i_dt, i_simulation_timestamp);
#else
	run_timestep_agrid_planedatacomplex(io_u, io_v, i_dt, i_simulation_timestamp);
#endif
}


#if SWEET_USE_PLANE_SPECTRAL_SPACE

void Burgers_Plane_TS_l_direct::run_timestep_agrid_planedata(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (simVars.disc.use_staggering)
		FatalError("Staggering not supported");

	if (i_dt < 0)
		FatalError("Burgers_Plane_TS_l_direct: Only constant time step size allowed");


	typedef std::complex<T> complex;

	T dt = i_dt;

	/*
	 * This implementation works directly on PlaneData
	 */

	io_u.request_data_spectral();
	io_v.request_data_spectral();

	for (std::size_t ik1 = 0; ik1 < io_u.planeDataConfig->spectral_data_size[1]; ik1++)
	{
		T k1;
		if (ik1 < io_u.planeDataConfig->spectral_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = (T)((int)ik1-(int)io_u.planeDataConfig->spectral_data_size[1]);

		for (std::size_t ik0 = 0; ik0 < io_u.planeDataConfig->spectral_data_size[0]; ik0++)
		{
			T k0 = (T)ik0;

			complex U[2];
			U[0] = io_u.spectral_get(ik1, ik0);
			U[1] = io_v.spectral_get(ik1, ik0);

			/*
			 * Eigenvalues
			 */
			complex lambda;

         T a = pi2()*pi2()/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0]);
         T b = pi2()*pi2()/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1]);
			lambda = -simVars.sim.viscosity*(k0*k0*a+k1*k1*b);

			for (int k = 0; k < 2; k++)
			{
				std::complex<T> &lam = lambda;

				std::complex<T> K;

				K = dt*lam;
				K = l_expcplx(K);

				U[k] = K*U[k];
			}

#if BURGERS_PLANE_TS_L_DIRECT_QUADPRECISION
			std::complex<double> tmp1(U[0].real(), U[0].imag());
			io_u.p_spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[1].real(), U[1].imag());
			io_v.p_spectral_set(ik1, ik0, tmp2);
#else
			io_u.p_spectral_set(ik1, ik0, U[0]);
			io_v.p_spectral_set(ik1, ik0, U[1]);
#endif
		}
	}

	io_u.spectral_zeroAliasingModes();
	io_v.spectral_zeroAliasingModes();
}

#endif


void Burgers_Plane_TS_l_direct::run_timestep_agrid_planedatacomplex(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt < 0)
		FatalError("Burgers_Plane_TS_l_direct: Only constant time step size allowed");

	typedef std::complex<T> complex;


#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	PlaneDataComplex i_u = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_u);
	PlaneDataComplex i_v = Convert_PlaneData_To_PlaneDataComplex::physical_convert(io_v);
#else
	PlaneDataComplex i_u = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_u);
	PlaneDataComplex i_v = Convert_PlaneData_To_PlaneDataComplex::spectral_convert(io_v);
#endif

	PlaneDataComplex o_u(io_u.planeDataConfig);
	PlaneDataComplex o_v(io_u.planeDataConfig);

	T dt = i_dt;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	o_u.spectral_space_data_valid = true;
	o_u.physical_space_data_valid = false;

	o_v.spectral_space_data_valid = true;
	o_v.physical_space_data_valid = false;
#endif

	for (std::size_t ik1 = 0; ik1 < i_u.planeDataConfig->spectral_complex_data_size[1]; ik1++)
	{
		T k1;
		if (ik1 < i_u.planeDataConfig->spectral_complex_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = -(T)((int)ik1-(int)i_u.planeDataConfig->spectral_complex_data_size[1]);

		for (std::size_t ik0 = 0; ik0 < i_u.planeDataConfig->spectral_complex_data_size[0]; ik0++)
		{
			T k0;
			if (ik0 < i_u.planeDataConfig->spectral_complex_data_size[0]/2)
				k0 = (T)ik0;
			else
				k0 = (T)((int)ik0-(int)i_u.planeDataConfig->spectral_complex_data_size[0]);

			complex U[2];
			U[0] = i_u.p_spectral_get(ik1, ik0);
			U[1] = i_v.p_spectral_get(ik1, ik0);

			/*
			 * Eigenvalues
			 */
			complex lambda;
         T a = pi2()*pi2()/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0]);
         T b = pi2()*pi2()/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1]);
			lambda = -simVars.sim.viscosity*(k0*k0*a+k1*k1*b);

			for (int k = 0; k < 2; k++)
			{
				std::complex<T> &lam = lambda;

				std::complex<T> K;

				K = dt*lam;
				K = l_expcplx(K);

				U[k] = K*U[k];
			}


#if BURGERS_PLANE_TS_L_DIRECT_QUADPRECISION
			std::complex<double> tmp1(U[0].real(), U[0].imag());
			o_u.p_spectral_set(ik1, ik0, tmp1);

			std::complex<double> tmp2(U[1].real(), U[1].imag());
			o_v.p_spectral_set(ik1, ik0, tmp2);
#else
			o_u.p_spectral_set(ik1, ik0, U[0]);
			o_v.p_spectral_set(ik1, ik0, U[1]);
#endif
		}
	}

#if SWEET_DEBUG
	o_u.test_realphysical();
	o_v.test_realphysical();
#endif

	o_u.spectral_zeroAliasingModes();
	o_v.spectral_zeroAliasingModes();

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	io_u = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::physical_convert(o_v);
#else
	io_u = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_u);
	io_v = Convert_PlaneDataComplex_To_PlaneData::spectral_convert_physical_real_only(o_v);
#endif
}


Burgers_Plane_TS_l_direct::Burgers_Plane_TS_l_direct(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	if (simVars.disc.use_staggering)
		planeDataGridMapping.setup(i_simVars, op.planeDataConfig);
}



Burgers_Plane_TS_l_direct::~Burgers_Plane_TS_l_direct()
{
}

