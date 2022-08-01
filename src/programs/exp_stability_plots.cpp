/*
 * test_exp_phi_functions.cpp
 *
 *  Created on: Jul 30, 2022
 *      Author: Martin Schreiber
 */

#include <rexi/EXPFunctions.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>


int main()
{
	/*
	 * We use the split linear ode
	 *
	 * u_t = \lambda_1 u + \lambda_2 u
	 *
	 * where we fix \lambda_1 to be purely oscillatory and the stiff part
	 */

	typedef double T;
	typedef std::complex<T> CT;

	EXPFunctions<T> phi0("phi0");
	EXPFunctions<T> phi1("phi1");
	EXPFunctions<T> phi2("phi2");

	EXPFunctions<T> ups1("ups1");
	EXPFunctions<T> ups2("ups2");
	EXPFunctions<T> ups3("ups3");

	CT lambda1 = CT(0, 1);

	double ext_x = 60;	// time step size
	double ext_y = 1;	// lambda 2
	double dx = ext_x/500;
	double dy = ext_y/500;

	//dx = std::max(dx, dy);
	//dy = std::max(dx, dy);

	std::cout << "dx: " << dx << std::endl;
	std::cout << "dy: " << dy << std::endl;

	for (int timestepping_order = 1; timestepping_order <= 4; timestepping_order++)
	{
		if (timestepping_order == 3)
			continue;

		std::ostringstream ss_stability;
		ss_stability << "output_stability_plot_order_" << timestepping_order << ".csv";
		std::cout << "Writing to " << ss_stability.str() << std::endl;


		std::ofstream ofile_stability(ss_stability.str());
		ofile_stability << std::setprecision(15);

		std::ostringstream ss_errors;
		ss_errors << "output_errors_plot_order_" << timestepping_order << ".csv";
		std::cout << "Writing to " << ss_errors.str() << std::endl;

		std::ofstream ofile_errors(ss_errors.str());
		ofile_errors << std::setprecision(15);

		ofile_stability << "0";
		ofile_errors << "0";

		for (double x = -ext_x; x < ext_x+dx*0.5; x += dx)
		{
			ofile_stability << "\t" << x;
			ofile_errors << "\t" << x;
		}


		ofile_stability << std::endl;
		ofile_errors << std::endl;

		bool first_row = true;
		for (double y = -ext_y; y <= ext_y+dy*0.5; y += dy)
		{
			if (!first_row)
			{
				ofile_stability << std::endl;
				ofile_errors << std::endl;
			}

			if (first_row)
				first_row = false;

			ofile_stability << y;
			ofile_errors << y;

			for (double x = -ext_x; x < ext_x+dx*0.5; x += dx)
			{
				/*
				 * Determine lambda2 coefficient
				 *
				 * We assume a purely oscillatory
				 */
				CT lambda2(0, y);

				/*
				 * The time step size of the nonlinear integration part is given by the x coordinate
				 */
				double dt = x;

				/*
				 * Compute analytical solution
				 */
				CT analytical_u = std::exp(dt*(lambda1 + lambda2));

				/*
				 * Compute ETDnRK solution
				 */
				CT etdnrk_u = 0;

				if (timestepping_order == 1)
				{
					CT Un = 1.0;

					/*
					 * U_{1}=\varphi_{0}(\Delta tL)U_{0}+\Delta t\varphi_{1}(\Delta tL)N(U_{0}).
					 */
					CT phi0_Un = phi0.eval(dt*lambda1)*Un;
					CT FUn = lambda2*Un;
					CT phi1_FUn = phi1.eval(dt*lambda1)*FUn;

					etdnrk_u = phi0_Un + dt*phi1_FUn;
				}
				else if (timestepping_order == 2)
				{
					CT Un = 1.0;

					/*
					 * A_{n}=\phi_0(\Delta tL)U_{n}+\Delta t\varphi_{1}(\Delta tL)F(U_{n},t_{n})
					 *
					 * U_{n+1}=A_{n}+\Delta t\varphi_{2}(\Delta tL)\left(F(A_{n},t_{n}+\Delta t)-F(U_{n},t_{n})\right)
					 */
					CT FUn = lambda2*Un;
					CT An = phi0.eval(dt*lambda1)*Un + dt*phi1.eval(dt*lambda1)*FUn;
					CT FAn = lambda2*An;
					etdnrk_u = An + dt*phi2.eval(dt*lambda1)*(FAn - FUn);
				}
				else if (timestepping_order == 4)
				{
					CT Un = 1.0;

					/*
					 * A_{n} = \phi_0\left(\frac{1}{2}\Delta tL\right)U_{n}+\frac{1}{2}\Delta t\varphi_{1}\left(\frac{1}{2}\Delta tL\right)F(U_{n},t_{n})
					 * B_{n} = \phi_0\left(\frac{1}{2}\Delta tL\right)U_{n}+\frac{1}{2}\Delta t\varphi_{1}\left(\frac{1}{2}\Delta tL\right)F(A_{n},t_{n}+\frac{1}{2}\Delta t)
					 * C_{n} = \phi_0\left(\frac{1}{2}\Delta tL\right)A_{n}+\frac{1}{2}\Delta t\varphi_{1}\left(\frac{1}{2}\Delta tL\right)\left(2F(B_{n},t_{n}+\frac{1}{2}\Delta t)-F(U_{n},t_{n})\right)
					 */
					CT An = phi0.eval(0.5*lambda1)*Un + 0.5*dt*phi1.eval(0.5*lambda1)*lambda2*Un;
					CT Bn = phi0.eval(0.5*lambda1)*Un + 0.5*dt*phi1.eval(0.5*lambda1)*lambda2*An;
					CT Cn = phi0.eval(0.5*lambda1)*An + 0.5*dt*phi1.eval(0.5*lambda1)*(2.0*lambda2*Bn - lambda2*Un);

					/*
					 * R_{0}	=	U_{n}
					 * R_{1}	=	F(U_{n},t_{n})
					 * R_{2}	=	F(A_{n},t_{n}+\frac{1}{2}\Delta t)+F(B_{n},t_{n}+\frac{1}{2}\Delta t)
					 * R_{3}	=	F(C_{n},t_{n}+\Delta t)
					 */
					CT R0 = Un;
					CT R1 = lambda2*Un;
					CT R2 = lambda2*An + lambda2*Bn;
					CT R3 = lambda2*Cn;

					/*
					 * U_{n+1} = \varphi_{0}\left(\Delta tL\right)R_{0}+\Delta t\left(\upsilon_{1}\left(\Delta tL\right)R_{1}+2\upsilon_{2}\left(\Delta tL\right)R_{2}+\upsilon_{3}\left(\Delta tL\right)R_{3}\right)
					 */
					etdnrk_u = phi0.eval(dt*lambda1)*R0 + dt*(ups1.eval(dt*lambda1)*R1 + 2.0*ups2.eval(dt*lambda1)*R2 + ups3.eval(dt*lambda1)*R3);
				}
				else
				{
					std::cerr << "This order is not supported" << std::endl;
					exit(1);
				}


				double stability = std::abs(etdnrk_u) - std::abs(analytical_u);
				double error = std::abs(etdnrk_u - analytical_u);

				/*
				 * OUTPUT
				 */
				ofile_stability << "\t";
				ofile_errors << "\t";

				// Values > 0 refer to an unstable method
				ofile_stability << stability;

#if 1
				ofile_errors << error;
#else
				double z = x + y*30;
				if (z < 1e-1)
					z = 0;
				ofile_errors << z*0.01;
#endif
			}
		}

		ofile_stability.close();
		ofile_errors.close();
	}

	return 0;

}
