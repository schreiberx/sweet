
//#if !SWEET_USE_SPECTRAL_SPACE
//	#error "Spectral space not activated"
//#endif

#if SWEET_GUI
#	error	"GUI not supported"
#endif



#include <sweet/DataArray.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/Operators2D.hpp>

#include <math.h>
#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

SimulationVariables simVars;

#if SWEET_DEBUG_MODE
#include <fenv.h>
static void __attribute__ ((constructor))
trapfpe ()
{
	/* Enable some exceptions.  At startup all exceptions are masked.  */

	feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

int main(int i_argc, char *i_argv[])
{
#if SWEET_DEBUG_MODE
	trapfpe();
#endif

	// override flag
	SimulationVariables simVars;
	simVars.disc.use_spectral_basis_diffs = true;

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
		return -1;

	if (simVars.disc.use_spectral_basis_diffs)
		std::cout << "Using spectral diffs" << std::endl;
	else
		std::cout << "Using kernel-based diffs" << std::endl;


	double prev_error_diff_x = 0;
	double prev_error_diff_y = 0;

	double prev_error_diff2_x = 0;
	double prev_error_diff2_y = 0;

	double prev_error_lap = 0;
	double prev_error_bilap = 0;

	/*
	 * iterate over resolutions, starting by res[0] given e.g. by program parameter -n
	 */
	std::size_t res_x = simVars.disc.res[0];
	std::size_t res_y = simVars.disc.res[1];

	//std::size_t max_res = 2048;
	std::size_t max_res = 128;

	if (res_x > max_res || res_y > max_res)
		max_res = std::max(res_x, res_y);

	for (; res_x <= max_res && res_y <= max_res; res_x *= 2, res_y *= 2)
	{

		double tolerance_increase = sqrt(res_x) + sqrt(res_y);

		/*
		 * error tolerance for machine accuracy
		 *
		 * We assume 1e-12 for double precision
		 */
		double eps = 1e-9*tolerance_increase;

		/*
		 * error tolerance for convergence
		 *
		 * Here, we are very patronizing due to flickering convergence for coarse solutions which
		 * are not really representable in the Fouerier space where the discretization errors
		 * are dominating.
		 */
		double eps_convergence = 1e-4;


		std::cout << "*************************************************************" << std::endl;
		std::cout << "Testing operators with resolution " << res_x << " x " << res_y << std::endl;
		std::cout << "*************************************************************" << std::endl;
		std::size_t res[2] = {res_x, res_y};

		simVars.disc.res[0] = res[0];
		simVars.disc.res[1] = res[1];
		simVars.reset();

		/*
		 * keep h in the outer regions to allocate it only once and avoid reinitialization of FFTW
		 */
		DataArray<2> h(res);


		{
			std::cout << "**********************************************" << std::endl;
			std::cout << "> Resolution (" << res_x << "x" << res_y << ")" << std::endl;
			std::cout << "> Domain size (" << simVars.sim.domain_size[0] << "x" << simVars.sim.domain_size[1] << ")" << std::endl;
			std::cout << "**********************************************" << std::endl;
			std::cout << "error tol = " << eps << std::endl;
			std::cout << "**********************************************" << std::endl;


			/**
			 * 2nd order differential operators on high mode functions (nyquist frequency)
			 *
			 * note, that the function on which the 2nd diff operator is computed on has
			 * to be scaled up be a factor of domain_size^2, since e.g.
			 *
			 * diff(sin(2 pi x / size), x, x) = 4.0 pi^2 sin(2 pi x / size) / size^2
			 */
			{
				std::cout << std::endl;
				std::cout << " Testing differentiation for different frequencies (including Nyquist)" << std::endl;

				DataArray<2> h_diff2_x(res);
				DataArray<2> h_diff2_y(res);
				DataArray<2> h_diff_x(res);
				DataArray<2> h_diff_y(res);
				DataArray<2> h_bilaplace(res);

				Operators2D op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

				double freq_x = 0;
				double freq_y = 0;

				//Nyquist freq
				std::size_t nyq=simVars.disc.res[0]/2;

				//Vary frequencies
				for (std::size_t k = 0; k <= 4; k++)
				{
					if(k==0) //Fix a given frequency
					{
						freq_x = 5;
						freq_y = 5;
					}
					else
					{    // Vary with k
						freq_x = ((double)k* (double) nyq)/ 4.0;
						freq_y = ((double)k* (double) nyq)/ 4.0;
					}

					//std::cout << freq_x << std::endl;
					for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					{
						for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
						{
							double x = (((double)i+0.5)/(double)simVars.disc.res[0])*simVars.sim.domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.res[1])*simVars.sim.domain_size[1];
							h.set(j, i,	sin(2.0*freq_x*M_PIl*x)*sin(2.0*freq_y*M_PIl*y)	);
							h_diff_x.set(j, i,
									2.0*freq_x*M_PIl*cos(2.0*freq_x*M_PIl*x)*sin(2.0*freq_y*M_PIl*y)/(simVars.sim.domain_size[0])
							);

							h_diff_y.set(j, i,
									2.0*freq_y*M_PIl*sin(2.0*freq_x*M_PIl*x)*cos(2.0*freq_y*M_PIl*y)/(simVars.sim.domain_size[1])
							);
							h_diff2_x.set(j, i,
									4.0*freq_x*freq_x*M_PIl*M_PIl*(-1.0)*sin(2.0*freq_x*M_PIl*x)*sin(2.0*freq_y*M_PIl*y)/(simVars.sim.domain_size[0]*simVars.sim.domain_size[0])
							);

							h_diff2_y.set(j, i,
									4.0*freq_y*freq_y*M_PIl*M_PIl*(-1.0)*sin(2.0*freq_x*M_PIl*x)*sin(2.0*freq_y*M_PIl*y)/(simVars.sim.domain_size[1]*simVars.sim.domain_size[1])
							);
						}
					}
					//This assumes freq_x = freq_y
					h_bilaplace=8.0*freq_x*freq_x*M_PIl*M_PIl*8.0*freq_x*freq_x*M_PIl*M_PIl*h;

					double err_x = (op.diff_c_x(h)-h_diff_x).reduce_maxAbs()/(2.0*freq_x*M_PIl); // .reduce_norm2_quad();//*normalization*(simVars.sim.domain_size[0])/(2.0*M_PIl);
					double err_y = (op.diff_c_y(h)-h_diff_y).reduce_maxAbs()/(2.0*freq_y*M_PIl); //.reduce_norm2_quad();//*normalization*(simVars.sim.domain_size[1])/(2.0*M_PIl);
					// diff2 normalization = 4.0 pi^2 / L^2
					double err2_x = (op.diff2_c_x(h)-h_diff2_x).reduce_maxAbs()/(2.0*freq_x*M_PIl)/(2.0*freq_x*M_PIl); //.reduce_norm2_quad();//*normalization*(simVars.sim.domain_size[0]*simVars.sim.domain_size[0])/(4.0*M_PIl*M_PIl);
					double err2_y = (op.diff2_c_y(h)-h_diff2_y).reduce_maxAbs()/(2.0*freq_y*M_PIl)/(2.0*freq_y*M_PIl); //.reduce_norm2_quad();//*normalization*(simVars.sim.domain_size[1]*simVars.sim.domain_size[1])/(4.0*M_PIl*M_PIl);

					double err_laplace = (op.laplace(h)-h_diff2_x-h_diff2_y).reduce_maxAbs()/(2.0*freq_y*M_PIl)/(2.0*freq_y*M_PIl);

					double err_bilaplace = (op.laplace(op.laplace(h))-h_bilaplace).reduce_maxAbs()/(8.0*freq_x*freq_x*M_PIl*M_PIl*8.0*freq_x*freq_x*M_PIl*M_PIl);

					if (simVars.disc.use_spectral_basis_diffs)
					{
						std::cout << "frequency = " << freq_x << " of " << simVars.disc.res[0]/2 << std::endl;
						std::cout << "error diff x = " << err_x << std::endl;
						std::cout << "error diff y = " << err_y << std::endl;
						std::cout << "error diff2 x = " << err2_x << std::endl;
						std::cout << "error diff2 y = " << err2_y << std::endl;
						std::cout << "error laplace = " << err_laplace << std::endl;
						std::cout << "error bilaplace = " << err_bilaplace << std::endl;

						if ( std::max({err_x, err_y, err2_x, err2_y, err_laplace})  > eps)
						{
							std::cerr << "SPEC: Error threshold for diff operators too high for spectral differentiation!" << std::endl;
							exit(-1);
						}
						if ( err_bilaplace  > eps*10000) //there is more error associated because of the magnitude of round off errors (forth order operator)?
						{
							std::cerr << "SPEC: Error threshold for diff operator bilaplacian too high for spectral differentiation!" << std::endl;
							exit(-1);
						}

					}
					else
					{
						if(k==0)
						{
							double conv_x = prev_error_diff_x/err_x;
							double conv_y = prev_error_diff_y/err_y;
							double conv2_x = prev_error_diff2_x/err2_x;
							double conv2_y = prev_error_diff2_y/err2_y;
							double conv_lap = prev_error_lap/err_laplace;
							double conv_bilap = prev_error_bilap/err_bilaplace;
							std::cout << "frequency x = " << freq_x << " of " << simVars.disc.res[0]/2 << std::endl;
							std::cout << "frequency y = " << freq_y << " of " << simVars.disc.res[1]/2 << std::endl;
							std::cout << "error diff x = " << err_x << std::endl;
							std::cout << "error diff y = " << err_y << std::endl;
							std::cout << "error diff2 x = " << err2_x << std::endl;
							std::cout << "error diff2 y = " << err2_y << std::endl;
							std::cout << "error laplacian = " << err_laplace<< std::endl;
							std::cout << "error bilaplacian = " << err_bilaplace<< std::endl;
							std::cout << "conv x = " << conv_x << std::endl;
							std::cout << "conv y = " << conv_y << std::endl;
							std::cout << "conv2 x = " << conv2_x << std::endl;
							std::cout << "conv2 y = " << conv2_y << std::endl;
							std::cout << "conv lap = " << conv_lap << std::endl;
							std::cout << "conv bilap = " << conv_bilap << std::endl;

							if ( std::min({conv_x, conv_y, conv2_x, conv2_y, conv_lap, conv_bilap}) != 0)
							{
								if (abs(std::min({conv_x, conv_y, conv2_x, conv2_y, conv_lap, conv_bilap})-4.0) > eps_convergence)
								{
									std::cerr << "Cart: Error threshold exceeded, no convergence given!" << std::endl;
									exit(-1);
								}
							}
							prev_error_diff_x = err_x;
							prev_error_diff_y = err_y;
							prev_error_diff2_x = err2_x;
							prev_error_diff2_y = err2_y;
							prev_error_lap=err_laplace;
							prev_error_bilap=err_bilaplace;
							break;

						}
					}
				}
			}

			std::cout << "TEST A: DONE" << std::endl;

			/*
			 * Test * operator and anti-aliasing
			 */

			{
				std::cout << std::endl;
				std::cout << " Testing multiplication and de-aliasing" << std::endl;
				std::cout << " ----------------------------------------" << std::endl;
				DataArray<2> h1(res);
				DataArray<2> h2(res);
				DataArray<2> h12(res);
				DataArray<2> h12_dealiased(res);
				DataArray<2> h12_alias(res);
				DataArray<2> h12_noalias(res);
				DataArray<2> h12_truncated(res);


				Operators2D op(simVars.disc.res, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs);

				//Nyquist freq
				std::size_t nyq=simVars.disc.res[0]/2;
				//Truncated nyquist freq
				std::size_t nyqtrunc=(int) floor(2.0*nyq/3.0);
				//Total num of freqs
				std::size_t n=simVars.disc.res[0];

				//dx, dy
				double dx=simVars.sim.domain_size[0]/simVars.disc.res[0];
				double dy=simVars.sim.domain_size[1]/simVars.disc.res[1];

				for (std::size_t k = 0; k < 7; k++)
				{

					//Frequencies for tests
					// Important: Make sure freq_1 > freq_2 always !!!!!!!!!!!
					double freq_1 = nyq/4.0+k*nyq/8.0; // 2 dx wave
					//double freq_1 = 2*nyq-freq_2-1;2*k;
					double freq_2 = nyq/4.0; // 10 dx wave

					//Product frequencies
					double freq_sum=freq_1+freq_2;
					double freq_sub=freq_1-freq_2;

					//Frequency info
					std::cout << "Freq_1 = " << freq_1 << ", freq_2 = " << freq_2
							<< ", freq_sum = " << freq_sum << ", freq_sub = " << freq_sub
							<< ", Nyq limit = " << nyq << ", Truncated Nyq limit = " << nyqtrunc
							<< std::endl;

					//Check aliasing on original spectrum
					double trunc_sum = 1.0;
					bool alias_present = 0;
					if(freq_sum > (double) nyq) //these modes cannot be represented on the grid, and will be aliased
					{
						trunc_sum=0.0;
						alias_present = 1;
					}

					if(alias_present)
					{
						std::cout << "Frequency " << freq_sum << " is not representable on this grid and will"
						 << " contaminate the mode " << n-freq_sum << " (alias on low frequency)" << std::endl;
					}
					else
					{
						std::cout << "No aliasing present on original multiplication spectrum" << std::endl;
					}

					//Check aliasing on truncated spectrum
					double trunc_freq1=1.0;
					double trunc_freq2=1.0;
					double trunc_freq_sum=1.0;
					bool alias_trunc_present=0;
					bool alias_trunc_mult_present=0;

					if(freq_1 > (double) nyqtrunc)  //these modes cannot be represented on the truncated spectrum
					{
						trunc_freq1=0.0;
						alias_trunc_present=1;
						std::cout << "Frequency " << freq_1 << " is not representable on the truncated spectrum so multiplication will be zero for truncated spectrum" << std::endl;
					}
					if(freq_2 > (double) nyqtrunc) //these modes cannot be represented on the truncated spectrum
					{
						trunc_freq2=0.0;
						alias_trunc_present=1;
						std::cout << "Frequency " << freq_2 << " is not representable on the truncated spectrum so multiplication will be zero for truncated spectrum" << std::endl;
					}
					if(freq_sum > (double) nyqtrunc) //these modes cannot be represented on the truncated spectrum
					{
						trunc_freq_sum=0.0;
						alias_trunc_mult_present=1;
						std::cout << "Frequency " << freq_sum << " is not representable on the truncated spectrum so multiplication will truncate this high mode" << std::endl;
					}
					if(!alias_trunc_present && !alias_trunc_mult_present)
						std::cout << "No aliasing present on truncated multiplication spectrum" << std::endl;

					//cos(a x) cos(b x)  = 1/2 (cos( (a-b) x)+cos( (a+b) x))
					for (std::size_t j = 0; j < simVars.disc.res[1]; j++)
					{
						for (std::size_t i = 0; i < simVars.disc.res[0]; i++)
						{
							double x = (double)i * dx;
							double y = (double)j * dy;
							h1.set(j, i,	cos(2.0*freq_1*M_PIl*x)*cos(2.0*freq_1*M_PIl*y)	);
							h2.set(j, i,	cos(2.0*freq_2*M_PIl*x)*cos(2.0*freq_2*M_PIl*y)	);

							//                            ok                      aliased
							h12.set(j, i,	0.5*(cos(2.0*(freq_sub)*M_PIl*x)+cos(2.0*(freq_sum)*M_PIl*x))
									*0.5*(cos(2.0*(freq_sub)*M_PIl*y)+cos(2.0*(freq_sum)*M_PIl*y)));

							// just the representable wave
							h12_noalias.set(j, i,	0.5*( cos(2.0*(freq_sub)*M_PIl*x)+trunc_sum*cos(2.0*(freq_sum)*M_PIl*x))
									*0.5*( cos(2.0*(freq_sub)*M_PIl*y)+trunc_sum*cos(2.0*(freq_sum)*M_PIl*y)));

							h12_truncated.set(j, i,	trunc_freq1*trunc_freq2*0.5*( cos(2.0*(freq_sub)*M_PIl*x)+trunc_freq_sum*cos(2.0*(freq_sum)*M_PIl*x))
									*0.5*( cos(2.0*(freq_sub)*M_PIl*y)+trunc_freq_sum*cos(2.0*(freq_sum)*M_PIl*y)));
						}
					}

					//Standard multiplication
					double err_mult = (h1*h2-h12).reduce_maxAbs();
					//Multiplication with dealiasing from * operator
					double err_mult_dealias = (h1*h2-h12_noalias).reduce_maxAbs();
					//Multiplication with dealiasing from mult function (truncation)
					double err_mult_dealias2 = (h1.mult(h2)-h12_truncated).reduce_maxAbs();


					// TODO: remove this
//					err_mult_dealias = (h1*h2-h12_truncated).reduce_maxAbs();

					//debug
					//std::cout << "h1" << std::endl;
					//(h1).printSpectrumNonZero();
					//std::cout << "h2" << std::endl;
					//(h2).printSpectrumNonZero();
					//std::cout << "h1*h2" << std::endl;
					//(h1*h2).printSpectrumNonZero();
					//std::cout << "h1.mult(h2)" << std::endl;
					//(h1.mult(h2)).printSpectrumNonZero();
					//std::cout << "h12" << std::endl;
					//(h12).printSpectrumNonZero();
					//std::cout << "h12 noalias" << std::endl;
					//(h12_noalias).printSpectrumIndex();
					//std::cout << "h12 alias" << std::endl;
					//(h12_alias).printSpectrumIndex();


#if 0
					std::cout << "error mult * with possibly aliased exact solution = " << err_mult << std::endl;
					std::cout << "error mult * with respect to dealised exact solution = " << err_mult_dealias << std::endl;
					std::cout << "error mult function with respect to truncated and dealiased exact solution = " << err_mult_dealias2 << std::endl;
#endif

#if SWEET_USE_SPECTRAL_SPACE && SWEET_USE_SPECTRAL_DEALIASING
					/**
					 * Check correct dealiasing.
					 * If error is too high, dealiasing obviously failed
					 */
					if ( err_mult_dealias  > eps)
					{
						std::cout << "error operator*(...) with possibly aliased exact solution = " << err_mult << std::endl;
						std::cout << "error operator*(...) with respect to dealised exact solution = " << err_mult_dealias << std::endl;
						std::cout << "error mult() function with respect to truncated and dealiased exact solution = " << err_mult_dealias2 << std::endl;

						std::cerr << " WARNING: threshold for multiplication * operator too high !" << std::endl;
						if(alias_present)
							std::cerr << "    Multiplication has alias but dealiasing not able to remove it or removed it incorrectly" << std::endl;
						else
							std::cerr << "    Multiplication dealiasing affected spectrum without need" << std::endl;
						std::cout << "    h1*h2 nonzero spectrum entries" << std::endl;
						(h1*h2).printSpectrumNonZero();
						assert(false);
						exit(-1);
					}

					if ( err_mult_dealias2  > eps)
					{
						std::cerr << " WARNING: error for multiplication function 'mult' too high !" << std::endl;
					}

#else
					if ( err_mult_dealias  > eps || err_mult_dealias2 > eps)
						std::cerr << " Turn on de-aliasing on compile time to analyse de-aliasing errors" << std::endl;

					if ( err_mult  > eps)
					{
						std::cerr << " Error threshold for multiplication operator too high !" << std::endl;
						exit(-1);
					}
#endif


					std::cout << " " << std::endl;


				}
			}

			std::cout << "TEST B: DONE" << std::endl;
		}
	}

	std::cout << "SUCCESSFULLY FINISHED -- check warnings" << std::endl;

	return 0;
}
