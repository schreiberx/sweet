/*
 * SPHTestSolutions_SPH.hpp
 *
 *  Created on: 24 Aug 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk> Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SPH_SPHERETESTSOLUTIONS_SPH_HPP_
#define SRC_INCLUDE_SPH_SPHERETESTSOLUTIONS_SPH_HPP_

#include <cassert>
#include <stdlib.h>

class SphereTestSolutions_SPH
{
	int tn;
	int tm;

public:
	SphereTestSolutions_SPH(int i_n, int i_m)	:
		tn(i_n),
		tm(i_m)
	{
	}


	void test_function__grid_gaussian(double lon, double mu, double &io_data)
	{
		// n=0, m=0  OK
		if (tn == 1 && tm == 0)
		{
			io_data = sqrt(2.0)/2.0;
		}

		// n=1, m=0  OK
		else if (tn == 1 && tm == 0)
		{
			io_data = sqrt(6.0)*mu/2.0;
		}

		// n=2, m=0  OK
		else if (tn == 2 && tm == 0)
		{
			io_data = sqrt(10.0)/4.0*(3.0*mu*mu-1.0);
		}

		// n=2, m=1
		else if (tn == 2 && tm == 1)
		{
			io_data = -mu/2.0*sqrt(-15.0*mu*mu + 15.0);
		}

		// n=2, m=2
		else if (tn == 2 && tm == 2)
		{
			io_data = sqrt(15.0)/4.0*(-mu*mu + 1.0);
		}

		// n=3, m=0  OK
		else if (tn == 3 && tm == 0)
		{
			io_data = sqrt(14.0)/4.0*mu*(5*mu*mu-3.0);
		}

		// n=3, m=1
		else if (tn == 3 && tm == 1)
		{
			io_data = sqrt(-42.0*mu*mu+42.0)*(-5.0*mu*mu/8.0 + 1.0/8.0);
		}

		// n=3, m=2  OK
		else if (tn == 3 && tm == 2)
		{
			io_data = sqrt(105.0)*mu/4.0*(-mu*mu+1.0);
		}

		// n=3, m=3  OK
		else if (tn == 3 && tm == 3)
		{
			io_data = -sqrt(70.0)/8.0*pow(-mu*mu+1.0, 3.0/2.0);
		}

		else
		{
			std::cerr << "NOT SUPPORTED" << std::endl;
			assert(false);
			exit(1);
		}

		io_data *= std::cos(lon*tm);

		io_data *= 1.0/std::sqrt(2.0*M_PI);
	}



	void correct_result_diff_mu__grid_gaussian(double lon, double mu, double &io_data)
	{
		/*
		 * List of results
		 */
		// n=0, m=0  OK
		if (tn == 1 && tm == 0)
		{
			io_data = 0;
		}

		// n=1, m=0  OK
		else if (tn == 1 && tm == 0)
		{
			io_data = sqrt(6.0)/2.0;
		}

		// n=2, m=0  OK
		else if (tn == 2 && tm == 0)
		{
			io_data = sqrt(10.0)/4.0*6.0*mu;
		}

		// n=2, m=1
		else if (tn == 2 && tm == 1)
		{
//							io_data = -mu/2.0*sqrt(-15.0*mu*mu + 15.0);
			io_data = (15.0*mu*mu-7.5)/(sqrt(15.0)*sqrt(1.0-mu*mu));
		}


#if 0						// n=2, m=2
		else if (tn == 2 && tm == 2)
		{
			io_data = sqrt(15.0)/4.0*(-mu*mu + 1.0);
		}

		// n=3, m=0  OK
		else if (tn == 3 && tm == 0)
		{
			io_data = sqrt(14.0)/4.0*mu*(5*mu*mu-3.0);
		}

		// n=3, m=1
		else if (tn == 3 && tm == 1)
		{
			io_data = sqrt(-42.0*mu*mu+42.0)*(-5.0*mu*mu/8.0 + 1.0/8.0);
		}

		// n=3, m=2  OK
		else if (tn == 3 && tm == 2)
		{
			io_data = sqrt(105.0)*mu/4.0*(-mu*mu+1.0);
		}

		// n=3, m=3  OK
		else if (tn == 3 && tm == 3)
		{
			io_data = -sqrt(70.0)/8.0*pow(-mu*mu+1.0, 3.0/2.0);
		}
#endif

		else
		{
			std::cerr << "NOT SUPPORTED" << std::endl;
			assert(false);
			exit(1);
		}

		io_data *= std::cos(lon*tm);

		io_data *= (1.0-mu*mu);

		io_data *= 1.0/std::sqrt(2.0*M_PI);
	}
};



#endif /* SRC_INCLUDE_SPH_SPHERETESTSOLUTIONS_SPH_HPP_ */
