/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_HELPERS_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_HELPERS_HPP_


#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include "../ShackPDESWESphere.hpp"


class PDESWESphereBenchmarks_HelperGeostropicBalance
{

	sweet::ShackDictionary *shackDict;
	sweet::SphereOperators *ops;

	ShackPDESWESphere *shackPDESWESphere;

	sweet::SphereData_Physical fg;

public:
	PDESWESphereBenchmarks_HelperGeostropicBalance()	:
		shackDict(nullptr),
		ops(nullptr),
		shackPDESWESphere(nullptr)
	{
	}

	void shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPDESWESphere = shackDict->getAutoRegistration<ShackPDESWESphere>();
	}

	void setup(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;

		if (shackPDESWESphere->sphere_use_fsphere)
			fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_nonlinear(
			sweet::SphereData_Spectral &i_vort,
			sweet::SphereData_Spectral &i_div,
			sweet::SphereData_Spectral &o_phi
	)
	{
		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::SphereData_Physical ug(o_phi.sphereDataConfig);
		sweet::SphereData_Physical vg(o_phi.sphereDataConfig);

		ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		sweet::SphereData_Physical vrtg = i_vort.toPhys();

		using namespace sweet;
		sweet::SphereData_Physical tmpg1 = ug*(vrtg+fg);
		sweet::SphereData_Physical tmpg2 = vg*(vrtg+fg);

		sweet::SphereData_Spectral tmpspec1(o_phi.sphereDataConfig);
		sweet::SphereData_Spectral tmpspec2(o_phi.sphereDataConfig);

		ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		sweet::SphereData_Spectral phispec = ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);

		o_phi = phispec;
	}



	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	void computeGeostrophicBalance_linear(
			sweet::SphereData_Spectral &i_vort,
			sweet::SphereData_Spectral &i_div,
			sweet::SphereData_Spectral &o_phi
	)
	{
		const sweet::SphereData_Config *sphereDataConfig = o_phi.sphereDataConfig;
		/*
		 * Setup Coriolis effect
		 */
		sweet::SphereData_Physical f(sphereDataConfig);
		f.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = 2.0*shackPDESWESphere->sphere_rotating_coriolis_omega*mu;
			}
		);

		/*
		 * Compute vorticity and divergence from velocities
		 */
		sweet::SphereData_Physical u(sphereDataConfig);
		sweet::SphereData_Physical v(sphereDataConfig);

		ops->vrtdiv_to_uv(i_vort, i_div, u, v);

		sweet::SphereData_Physical vrtg = i_vort.toPhys();

		sweet::SphereData_Physical tmpg1 = u*f;
		sweet::SphereData_Physical tmpg2 = v*f;

		sweet::SphereData_Spectral tmpspec1(sphereDataConfig);
		sweet::SphereData_Spectral tmpspec2(sphereDataConfig);

		ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		sweet::SphereData_Spectral phispec = ops->inv_laplace(tmpspec1);

		o_phi = shackPDESWESphere->h0*shackPDESWESphere->gravitation + phispec;
	}


	void clear()
	{
		fg.clear();
	}
};

#endif
