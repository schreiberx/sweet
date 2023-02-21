/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@Gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_HELPERS_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_HELPERS_HPP_


#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>


class SWESphereBenchmarks_helpers
{

	SimulationVariables *simVars;
	SphereOperators_SphereData *ops;

public:
	SWESphereBenchmarks_helpers(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)	:
		simVars(i_simVars),
		ops(i_ops)
	{

	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_nonlinear(
			SphereData_Spectral &i_vort,
			SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical ug(o_phi.sphereDataConfig);
		SphereData_Physical vg(o_phi.sphereDataConfig);

		ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		SphereData_Physical vrtg = i_vort.toPhys();

		SphereData_Physical tmpg1 = ug*(vrtg+ops->fg);
		SphereData_Physical tmpg2 = vg*(vrtg+ops->fg);

		SphereData_Spectral tmpspec1(o_phi.sphereDataConfig);
		SphereData_Spectral tmpspec2(o_phi.sphereDataConfig);

		ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		SphereData_Spectral phispec = ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);

		o_phi = phispec;
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			SphereData_Spectral &i_vort,
			SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		const SphereData_Config *sphereDataConfig = o_phi.sphereDataConfig;
		/*
		 * Setup Coriolis effect
		 */
		SphereData_Physical f(sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*simVars->sim.sphere_rotating_coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		SphereData_Physical u(sphereDataConfig);
		SphereData_Physical v(sphereDataConfig);

		ops->vrtdiv_to_uv(i_vort, i_div, u, v);

		SphereData_Physical vrtg = i_vort.toPhys();

		SphereData_Physical tmpg1 = u*f;
		SphereData_Physical tmpg2 = v*f;

		SphereData_Spectral tmpspec1(sphereDataConfig);
		SphereData_Spectral tmpspec2(sphereDataConfig);

		ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		SphereData_Spectral phispec = ops->inv_laplace(tmpspec1);

		o_phi = simVars->sim.h0*simVars->sim.gravitation + phispec;
	}


};

#endif
