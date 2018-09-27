/*
 * swe_spm_simple_vort_div.hpp
 *
 * Example ported from python code of Jeffrey Whitaker
 *
 *  Created on: 19 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *
 */

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataPhysical.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/FatalError.hpp>


#define THIS_DEBUG	0


SphereDataConfig sphereDataConfigInstance;
SphereDataConfig *sphereDataConfig = &sphereDataConfigInstance;

SimulationVariables simVars;


/*
 * This file is based on
 * https://gist.github.com/jswhit/3845307
 *
 * and resembles this code
 *
 * The rewritten python code can be found in
 *     ./swe_sph_simple_vort_div/
 * The output results in physical space should directly match
 */


class Sim
{
public:
	double umax = 80.;
	double phi0 = M_PI/7.;
	double phi1 = 0.5*M_PI - phi0;
	double phi2 = 0.25*M_PI;
	double en = std::exp(-4.0/std::pow((phi1-phi0), 2.0));
	double alpha = 1./3.;
	double beta = 1./15.;
	double hamp = 120.;


	SphereOperators op;


	//double simtime = 66400;
	double simtime = 664;
	double dt = 60;
	int itmax = 6*int(simtime/dt);

	SphereDataConfig *sphereDataConfig;
	//f = 2.*omega*np.sin(lats) # coriolis


#if THIS_DEBUG
	void outputMinMaxSum(
			const SphereDataPhysical i_data,
			const std::string i_prefix
	)
	{
		if (!debug)
			return;

		double min = i_data.physical_reduce_min();
		double max = i_data.physical_reduce_max();
		double sum = i_data.physical_reduce_sum_quad();
		double suminc = i_data.physical_reduce_sum_quad_increasing();

		std::cout << i_prefix << " | outputMinMaxSum:" << std::endl;
		std::cout << "		| min: " << min << std::endl;
		std::cout << "		| max: " << max << std::endl;
		std::cout << "		| sum: " << sum << std::endl;
		std::cout << "		| suminc: " << suminc << std::endl;

	}

	void outputSpecMinMaxSum(
			const SphereData i_data,
			const std::string i_prefix
	)
	{
		if (!debug)
			return;

		std::complex<double> min = i_data.spectral_reduce_min();
		std::complex<double> max = i_data.spectral_reduce_max();
		std::complex<double> sum = i_data.spectral_reduce_sum_quad();
		std::complex<double> suminc = i_data.spectral_reduce_sum_quad_increasing();

		std::cout << i_prefix << " | outputSpecMinMaxSum:" << std::endl;
		std::cout << "		| min: " << min.real() << " " << min.imag() << std::endl;
		std::cout << "		| max: " << max.real() << " " << max.imag() << std::endl;
		std::cout << "		| sum: " << sum.real() << " " << sum.imag() << std::endl;
		std::cout << "		| suminc: " << suminc.real() << " " << suminc.imag() << std::endl;
	}
#else

		#define outputSpecMinMaxSum(a,b)
		#define outputMinMaxSum(a,b)

#endif


	void timestep(
			const SphereData &vrtspec,
			const SphereData &divspec,
			const SphereData &phispec,
			SphereData &o_dvrt_dt,
			SphereData &o_ddiv_dt,
			SphereData &o_dphi_dt
	)
	{
		SphereDataPhysical ug(sphereDataConfig);
		SphereDataPhysical vg(sphereDataConfig);

		SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();
		op.vortdiv_to_uv(vrtspec, divspec, ug, vg);
		SphereDataPhysical phig = phispec.getSphereDataPhysical();

		outputMinMaxSum(vrtg, "vrtg");
		outputMinMaxSum(ug, "ug");
		outputMinMaxSum(vg, "vg");
		outputMinMaxSum(phig, "phig");

		SphereDataPhysical tmpg1 = ug*(vrtg+f);
		SphereDataPhysical tmpg2 = vg*(vrtg+f);
		outputMinMaxSum(tmpg1, "tmpg1");
		outputMinMaxSum(tmpg2, "tmpg2");

		op.uv_to_vortdiv(tmpg1, tmpg2, o_ddiv_dt, o_dvrt_dt);
		outputSpecMinMaxSum(o_ddiv_dt, "o_ddiv_dt");
		outputSpecMinMaxSum(o_dvrt_dt, "o_dvrt_dt");

		o_dvrt_dt *= -1.0;

		SphereDataPhysical tmpg = o_ddiv_dt.getSphereDataPhysical();
		outputMinMaxSum(tmpg, "tmpg");

		tmpg1 = ug*phig;
		tmpg2 = vg*phig;
		outputMinMaxSum(tmpg1, "tmpg1");
		outputMinMaxSum(tmpg2, "tmpg2");

		SphereData tmpspec(sphereDataConfig);
		op.uv_to_vortdiv(tmpg1,tmpg2, tmpspec, o_dphi_dt);

		outputSpecMinMaxSum(tmpspec, "tmpspec");
		outputSpecMinMaxSum(o_dphi_dt, "dphidtspec");

		o_dphi_dt *= -1.0;

		tmpspec = (phig+0.5*(ug*ug+vg*vg));
		tmpspec.request_data_spectral();
		o_ddiv_dt += -op.laplace(tmpspec);

		outputSpecMinMaxSum(o_dvrt_dt, "o_dvrt_dt");
		outputSpecMinMaxSum(o_ddiv_dt, "o_ddiv_dt");
		outputSpecMinMaxSum(o_dphi_dt, "o_dphi_dt");

		outputMinMaxSum(o_dvrt_dt.getSphereDataPhysical(), "o_dvrt_dt");
		outputMinMaxSum(o_ddiv_dt.getSphereDataPhysical(), "o_ddiv_dt");
		outputMinMaxSum(o_dphi_dt.getSphereDataPhysical(), "o_dphi_dt");
	}


	SphereDataPhysical f;
	SphereData phispec;
	SphereData vrtspec;
	SphereData divspec;

	SimulationVariables &simVars;

	Sim(
			SphereDataConfig *sphereDataConfig,
			SimulationVariables &i_simVars
	)	:
		sphereDataConfig(sphereDataConfig),
		f(sphereDataConfig),
		phispec(sphereDataConfig),
		vrtspec(sphereDataConfig),
		divspec(sphereDataConfig),
		simVars(i_simVars)
	{
		setup();
	}

	void setup()
	{
		op.setup(sphereDataConfig, simVars.sim.earth_radius);

		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*simVars.sim.coriolis_omega*mu;
			}
		);
		outputMinMaxSum(f, "f");



		SphereDataPhysical vg(sphereDataConfig);
		vg.physical_set_zero();

	    //u1 = (umax/en)*np.exp(1./((x.lats-phi0)*(x.lats-phi1)))
		SphereDataPhysical u1(sphereDataConfig);
		u1.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				o_data = umax/en*std::exp(1.0/((phi-phi0)*(phi-phi1)));
			}
		);
		outputMinMaxSum(u1, "u1");



		/*
	     ug = np.zeros((nlats),np.float)
	     ug = np.where(np.logical_and(x.lats < phi1, x.lats > phi0), u1, ug)
	     ug.shape = (nlats,1)
	     ug = ug*np.ones((nlats,nlons),dtype=np.float) # broadcast to shape (nlats,nlonss)
	    */

		SphereDataPhysical ug = u1;
		ug.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				if (phi >= phi1 || phi <= phi0)
					o_data = 0;
			}
		);
		outputMinMaxSum(ug, "ug");

		/*
		 * hbump = hamp*np.cos(lats)*np.exp(-(lons/alpha)**2)*np.exp(-(phi2-lats)**2/beta)
		 */

		SphereDataPhysical hbump(sphereDataConfig);
		hbump.physical_update_lambda(
			[&](double lon, double phi, double &o_data)
			{
				o_data = hamp*std::cos(phi)*std::exp(-(lon/alpha)*(lon/alpha))*std::exp(-(phi2-phi)*(phi2-phi)/beta);
			}
		);
		outputMinMaxSum(hbump, "hbump");


		/*
		 * vrtspec, divspec =  x.getvrtdivspec(ug,vg)
		 */
		op.uv_to_vortdiv(ug, vg, vrtspec, divspec);
		outputSpecMinMaxSum(vrtspec, "vrtspec");
		outputSpecMinMaxSum(divspec, "divspec");


		SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();
		outputMinMaxSum(vrtg, "vrtg");
		outputMinMaxSum(ug, "ug");
		outputMinMaxSum(vg, "vg");

		SphereDataPhysical tmpg1 = ug*(vrtg+f);
		SphereDataPhysical tmpg2 = vg*(vrtg+f);

		outputMinMaxSum((vrtg+f), "(vrtg+f)");
		outputMinMaxSum(ug*(vrtg+f), "ug*(vrtg+f)");
		outputMinMaxSum(vg*(vrtg+f), "vg*(vrtg+f)");

		outputMinMaxSum(tmpg1, "tmpg1");
		outputMinMaxSum(tmpg2, "tmpg2");

		SphereData tmpspec1(sphereDataConfig);
		SphereData tmpspec2(sphereDataConfig);
		op.uv_to_vortdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);
		outputSpecMinMaxSum(tmpspec1, "tmpspec1");
		outputSpecMinMaxSum(tmpspec2, "tmpspec2");

		tmpspec2 = (0.5*(ug*ug+vg*vg));
		outputSpecMinMaxSum(tmpspec2, "tmpspec2");

		phispec = op.inv_laplace(tmpspec1) - tmpspec2;
		outputSpecMinMaxSum(phispec, "phispec");


		simVars.sim.h0 = 10000.0;
		SphereDataPhysical phig = simVars.sim.gravitation*(hbump+simVars.sim.h0) + phispec.getSphereDataPhysical();
		outputMinMaxSum(phig, "phig");

		phispec = phig;
		phispec.spectral_truncate();
		outputSpecMinMaxSum(phispec, "phispec");

		outputMinMaxSum(phispec.getSphereDataPhysical(), "phispecphys");
	}



	void run()
	{
		SphereData phispecdt(sphereDataConfig);
		SphereData vrtspecdt(sphereDataConfig);
		SphereData divspecdt(sphereDataConfig);

		std::cout << "ITMAX: " << itmax << std::endl;
		for (int i = 0; i < itmax+1; i++)
		{
			std::cout << "TIMESTEP " << i << "   " << SphereData(phispec).getSphereDataPhysical().physical_reduce_max() << std::endl;

#if 0
			timestep(vrtspec, divspec, phispec, vrtspecdt, divspecdt, phispecdt);
#else
			timestep(vrtspec, divspec, phispec, vrtspecdt, divspecdt, phispecdt);
			timestep(vrtspec+0.5*dt*vrtspecdt, divspec+0.5*dt*divspecdt, phispec+0.5*dt*phispecdt, vrtspecdt, divspecdt, phispecdt);
#endif

			vrtspec += dt*vrtspecdt;
			divspec += dt*divspecdt;
			phispec += dt*phispecdt;
		}

		SphereDataPhysical vrtg = vrtspec.getSphereDataPhysical();
		SphereDataPhysical ug(sphereDataConfig);
		SphereDataPhysical vg(sphereDataConfig);
		op.vortdiv_to_uv(vrtspec, divspec, ug, vg);
		SphereDataPhysical phig = phispec.getSphereDataPhysical();

		outputMinMaxSum(vrtg, "vrtg");
		outputMinMaxSum(ug, "ug");
		outputMinMaxSum(vg, "vg");
		outputMinMaxSum(phig, "phig");

		vrtspec.physical_file_write("vrtspec_final.csv");
		divspec.physical_file_write("divspec_final.csv");
		phispec.physical_file_write("phispec_final.csv");
	}
};




int main(int i_argc, char *i_argv[])
{

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr))
	{
		return -1;
	}

	sphereDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral, simVars.misc.shtns_use_plans);
	std::cout << "N: " << simVars.disc.res_physical[0] << ", " << simVars.disc.res_physical[1] << std::endl;
	std::cout << "Nlm: " << sphereDataConfigInstance.shtns->nlm << std::endl;

	std::cout << std::setprecision(12);

	Sim sim(sphereDataConfig, simVars);
	sim.run();

}
