/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_IRK_HPP_


#include <complex>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "helpers/SWESphBandedMatrixPhysicalReal.hpp"
#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_erk.hpp"
#include "PDESWESphereTS_lg_erk.hpp"



/**
 * Implicit solver
 */
class PDESWESphereTS_l_irk	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		);

	bool setup(
			sweet::SphereOperators *io_ops,
			int i_timestep_order,
			double i_timestep_size
	);

	bool setup_main(
			sweet::SphereOperators *io_ops,
			int i_timestep_order,
			double i_timestep_size,
			double i_crank_nicolson_damping_factor,
			bool i_no_coriolis
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method);
	std::string getIDString();

	std::string timestepping_method;

private:
	PDESWESphereTS_lg_erk swe_sphere_ts_lg_erk;
	PDESWESphereTS_l_erk swe_sphere_ts_l_erk;

	SphBandedMatrixPhysicalReal< std::complex<double> > sphSolverDiv;

	double crank_nicolson_damping_factor;

	/// timestep size
	double timestep_size;

	/// individual time step size
	double dt_explicit = -1;
	double dt_implicit = -1;

	/// earth radius
	double sphere_radius;

	bool use_f_sphere;

	bool no_coriolis;

	/// f0
	double f0;

	/// Coriolis effect
	double two_coriolis;

	sweet::SphereData_Physical mug;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		swe_sphere_ts_lg_erk.shackRegistration(io_shackDict);
		swe_sphere_ts_l_erk.shackRegistration(io_shackDict);
		return true;
	}



public:
	PDESWESphereTS_l_irk();

private:
	void update_coefficients(double i_timestep_size);

public:
	void clear();

public:
	void runTimestep(
			sweet::SphereData_Spectral &io_phi,
			sweet::SphereData_Spectral &io_vrt,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void solveImplicit(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double dt
	);


	virtual ~PDESWESphereTS_l_irk();
};


#endif
