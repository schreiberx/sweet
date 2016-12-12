/*
 * REXI test optimizer
 */

#include <rexi/REXI.hpp>
#include <libmath/GaussQuadrature.hpp>
#include <sweet/SimulationVariables.hpp>


SimulationVariables simVars;


template <typename T>
class GaussianApproximationOnTheFly
{
	typedef std::complex<T> complex;



public:
	T mu;				///< average
	std::vector<complex> a;	///< weights for approximation
	int L;					///< 2*L+1 = number of weights

	T pi, pi2;
	T error_threshold;

	T currentError()
	{
		auto opt_function = [&](T i_x) -> T
		{
			return std::abs(this->approxGaussian(i_x, 1.0) - std::exp(i_x));
		};

		T error = GaussQuadrature::integrate5_intervals_adaptive_recursive<T>(
				-pi2, pi2,
				opt_function,
				error_threshold
				);

		return error;
	}



	void optimize(
			T &io_opti_value
	)
	{
		// search space
		T left_x_param = -pi2;
		T right_x_param = pi2;

		io_opti_value = left_x_param;
		T left_error = currentError();

		io_opti_value = right_x_param;
		T right_error = currentError();

		// previous error reduction
		T prev_error_change = std::numeric_limits<T>::max();
		T prev_error = std::numeric_limits<T>::max();

		while (true)
		{
			std::cout << "left/right: [" << left_x_param << ", " << left_error << "], [" << right_x_param << ", " << right_error << "]\t";
			// mid point
			T mid_x_param = (T)0.5*(left_x_param + right_x_param);
			io_opti_value = mid_x_param;
			T mid_error = currentError();

			// test for local extrema which is not yet the right solution
			prev_error_change = std::abs(prev_error-mid_error);
			prev_error = mid_error;

			std::cout << "Current error: " << mid_error << std::endl;

			if (prev_error_change < error_threshold)
			{
				std::cout << "Error change less than " << prev_error_change << std::endl;
				break;
			}

			if (left_error < right_error)
			{
				// use left segment
				right_x_param = mid_x_param;
				right_error = mid_error;
			}
			else
			{
				// use right segment
				left_x_param = mid_x_param;
				left_error = mid_error;
			}
		}
	}


	GaussianApproximationOnTheFly(
			int i_L = 11	///< L
	)
	{
		pi = std::acos((T)0)*(T)2.0;
		pi2 = pi*(T)2.0;

		L = i_L;

		mu = 0;

		error_threshold = 1e-6;

		a.resize(L*2+1);
		for (int i = 0; i < a.size(); i++)
		{
			a[i].real((T)rand() / (T)RAND_MAX);
			a[i].imag((T)rand() / (T)RAND_MAX);
		}

		T error = std::numeric_limits<T>::max();

		print();
		do
		{
			for (int i = 0; i < a.size(); i++)
			{
				for (int j = 0; j < 2; j++)
				{
					std::cout << "Optimizing for a [r" << i << ", c" << j << "]" << std::endl;
					T &opti_variable = reinterpret_cast<T(&)[2]>(a[i])[j];

					optimize(opti_variable);
				}
			}

			{
				std::cout << "Optimizing for mu" << std::endl;
				T &opti_variable = mu;

				optimize(opti_variable);
			}

			error = currentError();
		} while (error > error_threshold);

#if 0
		if (L == 11)
		{
			mu = {	-4.315321510875024, 0};

			L = 11;
			a.resize(L*2+1);
			a = {
					{	-1.0845749544592896e-7,		2.77075431662228e-8		},
					{	1.858753344202957e-8,		-9.105375434750162e-7	},
					{	3.6743713227243024e-6,		7.073284346322969e-7	},
					{	-2.7990058083347696e-6,		0.0000112564827639346	},
					{	0.000014918577548849352,	-0.0000316278486761932	},
					{	-0.0010751767283285608,		-0.00047282220513073084	},
					{	0.003816465653840016,		0.017839810396560574	},
					{	0.12124105653274578,		-0.12327042473830248	},
					{	-0.9774980792734348,		-0.1877130220537587		},
					{	1.3432866123333178,			3.2034715228495942		},
					{	4.072408546157305,			-6.123755543580666		},
					{	-9.442699917778205,			0.						},
					{	4.072408620272648,			6.123755841848161		},
					{	1.3432860877712938,			-3.2034712658530275		},
					{	-0.9774985292598916,		0.18771238018072134		},
					{	0.1212417070363373,			0.12326987628935386		},
					{	0.0038169724770333343,		-0.017839242222443888	},
					{	-0.0010756025812659208,		0.0004731874917343858	},
					{	0.000014713754789095218,	0.000031358475831136815	},
					{	-2.659323898804944e-6,		-0.000011341571201752273	},
					{	3.6970377676364553e-6,		-6.517457477594937e-7	},
					{	3.883933649142257e-9,		9.128496023863376e-7	},
					{	-1.0816457995911385e-7,		-2.954309729192276e-8	}
			};
		}
#endif
	}



	/**
	 * directly evaluate basis function which is to be approximated
	 */
	T evalGaussian(
			T x,	///< x-coefficient for Gaussian basis function
			T h	///< h-coefficient for Gaussian basis function
	)
	{
		return std::exp(-(x*x)/((T)4.0*h*h))/std::sqrt((T)4.0*pi);
	}


	/**
	 * evaluate approximation of Gaussian basis function
	 *
	 * with sum of complex rational functions
	 */
	T approxGaussian(
			T x,	///< x-coefficient for Gaussian basis function
			T h	///< h-coefficient for Gaussian basis function
	)
	{
		// scale x, since it depends linearly on h:
		// x^2 ~ h^2
		x /= h;

		T sum = 0;

		for (int l = 0; l < 2*L+1; l++)
		{
			int j = l-L;

			// WORKS with max error 7.15344e-13
			sum += (a[l]/(complex(0, x) + mu + complex(0, j))).real();
		}

		return sum;
	}

	void print()
	{
		std::cout << "mu: " << mu << std::endl;
		for (std::size_t i = 0; i < a.size(); i++)
		{
			std::cout << "a[" << i << "]: " << a[i] << std::endl;
		}
	}
};



int main(int i_argc, char *i_argv[])
{

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr))
		return -1;


	std::cout << std::endl;
	std::cout << "simVars.rexi.rexi_h: " << simVars.rexi.rexi_h << std::endl;
	std::cout << "simVars.rexi.rexi_M: " << simVars.rexi.rexi_M << std::endl;
	std::cout << "simVars.rexi.rexi_L: " << simVars.rexi.rexi_L << std::endl;
	std::cout << "simVars.rexi.rexi_half: " << simVars.rexi.rexi_use_half_poles << std::endl;
	std::cout << std::endl;


	GaussianApproximationOnTheFly<double> ga(11);

	return 0;
}

