/*
 * Author: Pedro Peixoto <ppeixoto@usp.br>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP_


#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>



class SWESphereBenchmark_barotropic_vort_modes	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	SWESphereBenchmark_barotropic_vort_modes()
	{
	}


	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
				benchmark_name == "barotropic_vort_modes"			||
				benchmark_name == "barotropic_vort_modes_1"	||
				benchmark_name == "barotropic_vort_modes_2"	||
				benchmark_name == "barotropic_vort_modes_3"	
		;
	}


	void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		simVars = i_simVars;
		ops = i_ops;
	}


	std::string get_help()
	{
		std::ostringstream stream;
		stream << "  BOROTROPIC VORTICITY MODES :" << std::endl;
		stream << "     'barotropic_vort_modes'" << std::endl;
		stream << "     'barotropic_vort_modes': standard test" << std::endl;
		stream << "     'barotropic_vort_modes_[N]_[M]': mode choices" << std::endl;
		return stream.str();
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	static void computeGeostrophicBalance_nonlinear(
			const SphereOperators_SphereData *i_ops,
			const SphereData_Spectral &i_vort,
			const SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		const SphereData_Config *sphereDataConfig = i_ops->sphereDataConfig;

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical ug(sphereDataConfig);
		SphereData_Physical vg(sphereDataConfig);

		i_ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical vrtg = i_vort.toPhys();

		SphereData_Physical tmpg1 = ug*(vrtg+i_ops->fg);
		SphereData_Physical tmpg2 = vg*(vrtg+i_ops->fg);

		SphereData_Spectral tmpspec1(sphereDataConfig);
		SphereData_Spectral tmpspec2(sphereDataConfig);

		i_ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		o_phi = i_ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);
	}




	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{

		/*
		 * barotropic vorticity initialization:
		 * uses normal modes, that are spherical harmonic modes
		 *
		 * WARNING: It is designed for the barotropic vorticity equations, but can be used for SWE
		 * considering a balanced solution for the full non-linear equations
		 *           
		 */
		
		if (simVars->timecontrol.current_simulation_time == 0)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
		simVars->sim.gravitation = 9.80616;
		simVars->sim.sphere_radius = 6.37122e6;
		simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

		ops->setup(ops->sphereDataConfig, &(simVars->sim));
	

		double a = simVars->sim.sphere_radius;
		double omega = simVars->sim.sphere_rotating_coriolis_omega;
		double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;


		double freq_multiplier = 1.0;

		if (benchmark_name == "geostrophic_balance_2")
			freq_multiplier = 2.0;
		else if (benchmark_name == "geostrophic_balance_4")
			freq_multiplier = 4.0;

		

		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
		
		SphereData_Physical ug, vg;

		std::complex<double> val = 1;
		o_vrt.spectral_set(1,1,val);


		//Set phi for steady state for SWE
		std::cout << "[MULE] barotropic_vort_analytical_setup: 1" << std::endl;

		computeGeostrophicBalance_nonlinear(
				ops,
				o_vrt,
				o_div,
				o_phi_pert
		);

	}
};

#endif
