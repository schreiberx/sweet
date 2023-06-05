/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_IRK_HPP_



#include <complex>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "../timeHelpers/SWESphBandedMatrixPhysicalReal.hpp"
#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_lg_erk.hpp"


class PDESWESphereTS_lg_irk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

public:
	bool setup(
		const sweet::SphereOperators *io_ops,
		int i_timestep_order,
		double i_timestepSize
	);

public:
	bool setup_main(
		const sweet::SphereOperators *io_ops,
		int i_timestep_order,
		double i_timestepSize,
		double i_crank_nicolson_damping_factor
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;


private:
	PDESWESphereTS_lg_erk lg_erk;

	/// alpha/beta (time step related component for implicit solver)
	double alpha;
	double beta;

	/// Crank-Nicolson damping factor
	double crank_nicolson_damping_factor = 0.5;

	/// timestep size
	double timestep_size;

	/// earth radius
	double r;

	/// inverse of earth radius
	double inv_r;

	/// Average geopotential
	double gh;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		lg_erk.shackRegistration(io_shackDict);
		return true;
	}


public:
	PDESWESphereTS_lg_irk();


public:
	void update_coefficients();


public:
	void runTimestep(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vort,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt = 0,
		double i_simulation_timestamp = -1
	) override;

	virtual ~PDESWESphereTS_lg_irk();
};


#endif
