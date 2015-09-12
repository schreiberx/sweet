/*
 * rexi_swe.hpp
 *
 *  Created on: 24 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#include "RexiSWE.hpp"

#ifndef SWEET_REXI_PARALLEL_SUM
#	define SWEET_REXI_PARALLEL_SUM 1
#endif

/**
 * Compute the REXI sum massively parallel *without* a parallelization with parfor in space
 */
#if SWEET_REXI_PARALLEL_SUM
#	include <omp.h>
#endif



RexiSWE::RexiSWE()	:
	perThreadVars(nullptr)
{
#if !SWEET_USE_LIBFFT
	std::cerr << "Spectral space required for solvers, use compile optino --libfft=enable" << std::endl;
	exit(-1);
#endif

#if SWEET_REXI_PARALLEL_SUM
	#pragma omp parallel
	#pragma omp master
		num_threads = omp_get_num_threads();
#else
	num_threads = 1;
#endif
}

RexiSWE::~RexiSWE()
{
	if (perThreadVars != nullptr)
		delete [] perThreadVars;
}



/**
 * setup the REXI
 */
void RexiSWE::setup(
		double i_tau,			///< time step size
		double i_h,				///< sampling size
		int i_M,				///< number of sampling points
		int i_L,				///< number of sampling points for Gaussian approx
		double i_f,				///< Coriolis force
		std::size_t *i_resolution,		///< resolution of domain
		const double *i_domain_size,	///< size of domain
		bool i_rexi_half,		///< use half-pole reduction
		bool i_use_finite_differences
)
{
	M = i_M;
	h = i_h;
	tau = i_tau;
	f = i_f;

	rexi.setup(h, M, i_L, i_rexi_half);



	std::size_t N = rexi.alpha.size();
	block_size = N/num_threads;
	if (block_size*num_threads != N)
		block_size++;

	if (perThreadVars != nullptr)
		delete [] perThreadVars;

	perThreadVars = new PerThreadVars[num_threads];

	/**
	 * We split the setup from the utilization here.
	 *
	 * This is necessary, since it has to be assured that
	 * the FFTW plans are initialized before using them.
	 */
#if SWEET_REXI_PARALLEL_SUM
#	pragma omp parallel for schedule(static)
#endif
	for (int i = 0; i < num_threads; i++)
	{
		perThreadVars[i].op_diff_c_x.setup(i_resolution);
		perThreadVars[i].op_diff_c_y.setup(i_resolution);
		perThreadVars[i].op_diff2_c_x.setup(i_resolution);
		perThreadVars[i].op_diff2_c_y.setup(i_resolution);
		perThreadVars[i].eta0.setup(i_resolution);
		perThreadVars[i].u0.setup(i_resolution);
		perThreadVars[i].v0.setup(i_resolution);
		perThreadVars[i].h_sum.setup(i_resolution);
		perThreadVars[i].u_sum.setup(i_resolution);
		perThreadVars[i].v_sum.setup(i_resolution);
	}


#if SWEET_REXI_PARALLEL_SUM
#	pragma omp parallel for schedule(static)
#endif
	for (int i = 0; i < num_threads; i++)
	{
		// initialize all values to account for first touch policy
		perThreadVars[i].op_diff_c_x.setAll(0, 0);
		perThreadVars[i].op_diff_c_x.op_setup_diff_x(i_domain_size, i_use_finite_differences);

		perThreadVars[i].op_diff_c_y.setAll(0, 0);
		perThreadVars[i].op_diff_c_y.op_setup_diff_y(i_domain_size, i_use_finite_differences);

		perThreadVars[i].op_diff2_c_x.setAll(0, 0);
		perThreadVars[i].op_diff2_c_x.op_setup_diff2_x(i_domain_size, i_use_finite_differences);

		perThreadVars[i].op_diff2_c_y.setAll(0, 0);
		perThreadVars[i].op_diff2_c_y.op_setup_diff2_y(i_domain_size, i_use_finite_differences);

		perThreadVars[i].eta0.setAll(0, 0);
		perThreadVars[i].u0.setAll(0, 0);
		perThreadVars[i].v0.setAll(0, 0);

		perThreadVars[i].h_sum.setAll(0, 0);
		perThreadVars[i].u_sum.setAll(0, 0);
		perThreadVars[i].v_sum.setAll(0, 0);
	}
}



/**
 * Solve the REXI of \f$ U(t) = exp(L*t) \f$
 *
 * See
 * 		doc/rexi/understanding_rexi.pdf
 * for further information
 */
void RexiSWE::run_timestep(
	DataArray<2> &io_h,
	DataArray<2> &io_u,
	DataArray<2> &io_v,

	Operators2D &op,
	const SimulationVariables &i_parameters
)
{
	typedef std::complex<double> complex;


	std::size_t N = rexi.alpha.size();

#if SWEET_REXI_PARALLEL_SUM
#	pragma omp parallel for schedule(static) default(none) shared(i_parameters, io_h, io_u, io_v, N)
#endif
	for (int i = 0; i < num_threads; i++)
	{
		double eta_bar = i_parameters.setup.h0;
		double g = i_parameters.sim.g;

		Complex2DArrayFFT &op_diff_c_x = perThreadVars[i].op_diff_c_x;
		Complex2DArrayFFT &op_diff_c_y = perThreadVars[i].op_diff_c_y;
		Complex2DArrayFFT &op_diff2_c_x = perThreadVars[i].op_diff2_c_x;
		Complex2DArrayFFT &op_diff2_c_y = perThreadVars[i].op_diff2_c_y;

		Complex2DArrayFFT &eta0 = perThreadVars[i].eta0;
		Complex2DArrayFFT &u0 = perThreadVars[i].u0;
		Complex2DArrayFFT &v0 = perThreadVars[i].v0;

		Complex2DArrayFFT &h_sum = perThreadVars[i].h_sum;
		Complex2DArrayFFT &u_sum = perThreadVars[i].u_sum;
		Complex2DArrayFFT &v_sum = perThreadVars[i].v_sum;

		/*
		 * INITIALIZATION - THIS IS THE NON-PARALLELIZABLE PART!
		 */
		h_sum.setAll(0, 0);
		u_sum.setAll(0, 0);
		v_sum.setAll(0, 0);

		eta0.loadRealFromDataArray(io_h);
		u0.loadRealFromDataArray(io_u);
		v0.loadRealFromDataArray(io_v);

		// convert to spectral space
		// scale with inverse of tau
		eta0 = eta0.toSpec()*(1.0/tau);
		u0 = u0.toSpec()*(1.0/tau);
		v0 = v0.toSpec()*(1.0/tau);

#if SWEET_REXI_PARALLEL_SUM
		int thread_id = omp_get_thread_num();

		std::size_t start = block_size*thread_id;
		std::size_t end = std::min(N, start+block_size);

#else

		std::size_t start = 0;
		std::size_t end = N;

#endif

		/*
		 * DO SUM IN PARALLEL
		 */
		for (std::size_t n = start; n < end; n++)
		{
			// load alpha (a) and scale by inverse of tau
			// we flip the sign to account for the -L used in exp(\tau (-L))
			complex alpha = -rexi.alpha[n]/tau;
			complex beta = -rexi.beta_re[n];

			// load kappa (k)
			complex kappa = alpha*alpha + f*f;

			// compute
			// 		kappa - g * eta_bar * D2
			// NOTE!!! We add kappa in Cartesian space, hence add this value to all frequency components to account for scaling all frequencies!!!
			// This is *NOT* straightforward and different to adding a constant for computations.
			// We account for this by seeing the LHS as a set of operators which have to be joint later by a sum.
			Complex2DArrayFFT lhs = (-g*eta_bar*(op_diff2_c_x + op_diff2_c_y)).addScalar_Cart(kappa);

			Complex2DArrayFFT rhs =
					(kappa/alpha) * eta0
					- eta_bar*(op_diff_c_x(u0) + op_diff_c_y(v0))
					- (f*eta_bar/alpha) * (op_diff_c_x(v0) - op_diff_c_y(u0))
				;

			Complex2DArrayFFT eta = rhs.spec_div_element_wise(lhs);

			Complex2DArrayFFT uh = u0 - g*op_diff_c_x(eta);
			Complex2DArrayFFT vh = v0 - g*op_diff_c_y(eta);

			Complex2DArrayFFT u1 = alpha/kappa * uh     + f/kappa * vh;
			Complex2DArrayFFT v1 =    -f/kappa * uh + alpha/kappa * vh;

			h_sum += eta.toCart()*beta;
			u_sum += u1.toCart()*beta;
			v_sum += v1.toCart()*beta;
		}
	}

#if SWEET_REXI_PARALLEL_SUM
	io_h.set_all(0);
	io_u.set_all(0);
	io_v.set_all(0);
	for (int n = 0; n < num_threads; n++)
	{
		#pragma omp parallel for schedule(static)
		for (std::size_t i = 0; i < io_h.array_data_cartesian_length; i++)
		{
			io_h.array_data_cartesian_space[i] += perThreadVars[n].h_sum.data[i<<1];
			io_u.array_data_cartesian_space[i] += perThreadVars[n].u_sum.data[i<<1];
			io_v.array_data_cartesian_space[i] += perThreadVars[n].v_sum.data[i<<1];
		}
	}
#else

	io_h = perThreadVars[0].h_sum.getRealWithDataArray();
	io_u = perThreadVars[0].u_sum.getRealWithDataArray();
	io_v = perThreadVars[0].v_sum.getRealWithDataArray();

#endif
}



inline std::complex<double> conj(const std::complex<double> &v)
{
	return std::complex<double>(v.real(), -v.imag());
}



/**
 * This method computes the analytical solution based on the given initial values.
 *
 * See Embid/Madja/1996, Terry/Beth/2014, page 16
 * and
 * 		doc/swe_solution_for_L/sympy_L_spec_decomposition.py
 * for the dimensionful formulation.
 *
 * Don't use this function to frequently, since it always computes
 * the required coefficients on-the-fly which is expensive.
 */
void RexiSWE::run_timestep_direct_solution(
		DataArray<2> &io_h,
		DataArray<2> &io_u,
		DataArray<2> &io_v,

		double i_timestep_size,

		Operators2D &op,
		const SimulationVariables &i_parameters
)
{
	typedef std::complex<double> complex;

	double eta_bar = i_parameters.setup.h0;
	double g = i_parameters.sim.g;
	double f = i_parameters.sim.f0;
	complex I(0.0,1.0);

	Complex2DArrayFFT i_h(io_h.resolution);
	Complex2DArrayFFT i_u(io_h.resolution);
	Complex2DArrayFFT i_v(io_h.resolution);

	Complex2DArrayFFT o_h(io_h.resolution);
	Complex2DArrayFFT o_u(io_h.resolution);
	Complex2DArrayFFT o_v(io_h.resolution);

	i_h.loadRealFromDataArray(io_h);
	i_h = i_h.toSpec();

	i_u.loadRealFromDataArray(io_u);
	i_u = i_u.toSpec();

	i_v.loadRealFromDataArray(io_v);
	i_v = i_v.toSpec();

	double s0 = i_parameters.sim.domain_size[0];
	double s1 = i_parameters.sim.domain_size[1];

	for (std::size_t ik1 = 0; ik1 < i_h.resolution[1]; ik1++)
	{
		for (std::size_t ik0 = 0; ik0 < i_h.resolution[0]; ik0++)
		{
			if (ik0 == i_h.resolution[0]/2 || ik1 == i_h.resolution[1]/2)
			{
				o_h.set(ik1, ik0, 0, 0);
				o_u.set(ik1, ik0, 0, 0);
				o_v.set(ik1, ik0, 0, 0);
			}

			complex U_hat[3];
			U_hat[0] = i_h.get(ik1, ik0);
			U_hat[1] = i_u.get(ik1, ik0);
			U_hat[2] = i_v.get(ik1, ik0);

			double k0, k1;
			if (ik0 < i_h.resolution[0]/2)
				k0 = (double)ik0;
			else
				k0 = (double)((int)ik0-(int)i_h.resolution[0]);

			if (ik1 < i_h.resolution[1]/2)
				k1 = (double)ik1;
			else
				k1 = (double)((int)ik1-(int)i_h.resolution[1]);

			/*
			 * dimensionful formulation
			 * see doc/swe_solution_for_L
			 */

			double H0 = eta_bar;

			//////////////////////////////////////
			// GENERATED CODE START
			//////////////////////////////////////
			complex eigenvalues[3];
			complex eigenvectors[3][3];

			if (k0 == 0 && k1 == 0)
			{
//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1);

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*f;
				eigenvalues[2] = f;

				eigenvectors[0][0] = 1.00000000000000;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = 0.0;
				eigenvectors[1][1] = -1.0*I;
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 0.0;
				eigenvectors[2][1] = I;
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k0 == 0)
			{
//					complex wg = std::sqrt((complex)s0*s0*(f*f*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k1*k1 - 1.0*f*f*s1*s1);

				eigenvectors[0][0] = (1.0/2.0)*I*1.0/M_PI*f*1.0/g*1.0/k1*s1;
				eigenvectors[0][1] = 1.00000000000000;
				eigenvectors[0][2] = 0.0;
				eigenvectors[1][0] = -2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = -1.0*I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*M_PI*H0*k1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = I*f*s1/std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else if (k1 == 0)
			{
//					complex wg = std::sqrt((complex)s1*s1*(f*f*s0*s0 + 4.0*M_PI*M_PI*g*g*k0*k0));

				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvalues[2] = -1.0*I*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = 0.0;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[1][1] = -1.0*I*1.0/f*1.0/s0*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0 + f*f*s0*s0);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = 2.0*I*M_PI*H0*1.0/f*k0*1.0/s0;
				eigenvectors[2][1] = 1.0/f*1.0/s0*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0 - 1.0*f*f*s0*s0);
				eigenvectors[2][2] = 1.00000000000000;
			}
			else
			{
//					complex K2 = M_PI*M_PI*k0*k0 + M_PI*M_PI*k1*k1;
				complex w = std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);

//					complex wg = std::sqrt((complex)f*f*s0*s0*s1*s1 + 4.0*M_PI*M_PI*g*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*g*g*k1*k1*s0*s0);
				eigenvalues[0] = 0.0;
				eigenvalues[1] = -1.0*1.0/s0*1.0/s1*std::sqrt((complex)4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 + 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 + f*f*s0*s0*s1*s1);
				eigenvalues[2] = -1.0*I*1.0/s0*1.0/s1*std::sqrt((complex)-4.0*M_PI*M_PI*H0*g*k0*k0*s1*s1 - 4.0*M_PI*M_PI*H0*g*k1*k1*s0*s0 - 1.0*f*f*s0*s0*s1*s1);

				eigenvectors[0][0] = -1.0/2.0*I*1.0/M_PI*f*1.0/g*1.0/k0*s0;
				eigenvectors[0][1] = -1.0*1.0/k0*k1*s0*1.0/s1;
				eigenvectors[0][2] = 1.00000000000000;
				eigenvectors[1][0] = 2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 + f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[1][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 - 1.0*I*f*w);
				eigenvectors[1][2] = 1.00000000000000;
				eigenvectors[2][0] = -2.0*M_PI*H0*1.0/s0*1.0/w*(I*k0*s1*s1*(4.0*I*M_PI*M_PI*H0*g*k0*k1 - 1.0*f*w) - 1.0*k1*s0*s0*(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1))*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1);
				eigenvectors[2][1] = 1.0/s0*s1*1.0/(4.0*M_PI*M_PI*H0*g*k1*k1 + f*f*s1*s1)*(4.0*M_PI*M_PI*H0*g*k0*k1 + I*f*w);
				eigenvectors[2][2] = 1.00000000000000;
			}




			//////////////////////////////////////
			// GENERATED CODE END
			//////////////////////////////////////


			if (f == 0)
			{
				/*
				 * override if f == 0, see ./sympy_L_spec_decomposition.py executed with LNr=4
				 */
				if (k0 != 0 || k1 != 0)
				{
					double K2 = K2;

					eigenvalues[0] = 0.0;
					eigenvalues[1] = -2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);
					eigenvalues[2] = 2.0*M_PI*sqrt(H0)*sqrt((double)g)*sqrt(k0*k0 + k1*k1);

					eigenvectors[0][0] = 0.0;
					eigenvectors[0][1] = -1.0*k1/sqrt(k0*k0 + k1*k1);
					eigenvectors[0][2] = k0/sqrt(k0*k0 + k1*k1);
					eigenvectors[1][0] = -1.0*sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[1][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][0] = sqrt(H0)*sqrt(k0*k0 + k1*k1)/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][1] = sqrt((double)g)*k0/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
					eigenvectors[2][2] = sqrt((double)g)*k1/sqrt(H0*(k0*k0 + k1*k1) + g*k0*k0 + g*k1*k1);
				}
				else
				{

					eigenvalues[0] = 0.0;
					eigenvalues[1] = 0.0;
					eigenvalues[2] = 0.0;

					eigenvectors[0][0] = 1.00000000000000;
					eigenvectors[0][1] = 0.0;
					eigenvectors[0][2] = 0.0;
					eigenvectors[1][0] = 0.0;
					eigenvectors[1][1] = 1.00000000000000;
					eigenvectors[1][2] = 0.0;
					eigenvectors[2][0] = 0.0;
					eigenvectors[2][1] = 0.0;
					eigenvectors[2][2] = 1.00000000000000;
				}
			}


			/*
			 * Compute inverse of Eigenvectors.
			 * This generalizes to the case that the Eigenvectors are not orthonormal.
			 */
			complex eigenvectors_inv[3][3];

			eigenvectors_inv[0][0] =  (eigenvectors[1][1]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][1] = -(eigenvectors[0][1]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][1]);
			eigenvectors_inv[0][2] =  (eigenvectors[0][1]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][1]);

			eigenvectors_inv[1][0] = -(eigenvectors[1][0]*eigenvectors[2][2] - eigenvectors[1][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][1] =  (eigenvectors[0][0]*eigenvectors[2][2] - eigenvectors[0][2]*eigenvectors[2][0]);
			eigenvectors_inv[1][2] = -(eigenvectors[0][0]*eigenvectors[1][2] - eigenvectors[0][2]*eigenvectors[1][0]);

			eigenvectors_inv[2][0] =  (eigenvectors[1][0]*eigenvectors[2][1] - eigenvectors[1][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][1] = -(eigenvectors[0][0]*eigenvectors[2][1] - eigenvectors[0][1]*eigenvectors[2][0]);
			eigenvectors_inv[2][2] =  (eigenvectors[0][0]*eigenvectors[1][1] - eigenvectors[0][1]*eigenvectors[1][0]);

			complex s = eigenvectors[0][0]*eigenvectors_inv[0][0] + eigenvectors[0][1]*eigenvectors_inv[1][0] + eigenvectors[0][2]*eigenvectors_inv[2][0];

			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					eigenvectors_inv[j][i] /= s;


			// check
			for (int j = 0; j < 3; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					if (
							std::isnan(eigenvectors[j][i].real()) || std::isinf(eigenvectors[j][i].real())	||
							std::isnan(eigenvectors[j][i].imag()) || std::isinf(eigenvectors[j][i].imag())
					)
					{
						std::cerr << "Invalid number in Eigenvector " << j << " detected: " << eigenvectors[j][0] << ", " << eigenvectors[j][1] << ", " << eigenvectors[j][2] << std::endl;
					}

					if (
							std::isnan(eigenvectors_inv[j][i].real()) || std::isinf(eigenvectors_inv[j][i].real())	||
							std::isnan(eigenvectors_inv[j][i].imag()) || std::isinf(eigenvectors_inv[j][i].imag())
					)
					{
						std::cerr << "Invalid number in inverse of Eigenvector " << j << " detected: " << eigenvectors_inv[j][0] << ", " << eigenvectors_inv[j][1] << ", " << eigenvectors_inv[j][2] << std::endl;
					}
				}
			}

			/*
			 * Solve based on previously computed data.
			 * Note, that this data can be also precomputed and reused every time.
			 */
			complex UEV0_sp[3];
			for (int k = 0; k < 3; k++)
			{
				UEV0_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					UEV0_sp[k] += eigenvectors_inv[j][k] * U_hat[j];
			}

			complex omega[3];
			omega[0] = std::exp(-I*eigenvalues[0]*i_timestep_size);
			omega[1] = std::exp(-I*eigenvalues[1]*i_timestep_size);
			omega[2] = std::exp(-I*eigenvalues[2]*i_timestep_size);

			complex U_hat_sp[3];
			for (int k = 0; k < 3; k++)
			{
				U_hat_sp[k] = {0, 0};
				for (int j = 0; j < 3; j++)
					U_hat_sp[k] += eigenvectors[j][k] * omega[j] * UEV0_sp[j];
			}

			o_h.set(ik1, ik0, U_hat_sp[0]);
			o_u.set(ik1, ik0, U_hat_sp[1]);
			o_v.set(ik1, ik0, U_hat_sp[2]);
		}
	}

	io_h = o_h.toCart().getRealWithDataArray();
	io_u = o_u.toCart().getRealWithDataArray();
	io_v = o_v.toCart().getRealWithDataArray();
}
