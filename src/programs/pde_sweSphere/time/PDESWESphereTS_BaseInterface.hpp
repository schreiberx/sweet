/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_BASE_INTERFACE_NEW_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_BASE_INTERFACE_NEW_HPP_


#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/sphere/Sphere.hpp>

// Shacks
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/time/ShackTimesteppingSemiLagrangianSphereData.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWESphereTimeDiscretization.hpp"
#include "../benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "../ShackPDESWESphere.hpp"



#if SWEET_PARAREAL || SWEET_XBRAID
#include <sweet/parareal/Parareal_GenericData.hpp>
#endif


class PDESWESphereTS_BaseInterface
{
public:
	sweet::ErrorBase error;

	/*
	 * These are just some default shacks we provide to each time stepping method
	 */

	sweet::ShackDictionary *shackDict;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackExpIntegration *shackExpIntegration;
	sweet::ShackTimesteppingSemiLagrangianSphereData *shackTimesteppingSemiLagrangianSphereData;

	ShackPDESWESphereTimeDiscretization *shackPDESWETimeDisc;
	ShackPDESWESphereBenchmarks *shackPDESWEBenchmark;
	ShackPDESWESphere *shackPDESWESphere;

	sweet::SphereOperators *ops;

	std::string timestepping_method;
	int timestepping_order;
	int timestepping_order2;

	sweet::SphereData_Physical fg;

	PDESWESphereTS_BaseInterface()	:
		shackDict(nullptr),
		shackTimestepControl(nullptr),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr),
		shackExpIntegration(nullptr),
		shackPDESWETimeDisc(nullptr),
		shackPDESWEBenchmark(nullptr),
		shackPDESWESphere(nullptr),
		ops(nullptr),
		timestepping_order(-1),
		timestepping_order2(-1)
	{

	}

	virtual ~PDESWESphereTS_BaseInterface()
	{
	}


	bool setupFG()
	{
		assert(shackPDESWESphere != nullptr);
		if (shackPDESWESphere->sphere_use_fsphere)
			fg = ops->getFG_fSphere(shackPDESWESphere->sphere_fsphere_f0);
		else
			fg = ops->getFG_rotatingSphere(shackPDESWESphere->sphere_rotating_coriolis_omega);

		return true;
	}


	virtual bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackTimestepControl = io_shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackSphereDataOps = io_shackDict->getAutoRegistration<sweet::ShackSphereDataOps>();
		shackIOData = io_shackDict->getAutoRegistration<sweet::ShackIOData>();
		shackExpIntegration = io_shackDict->getAutoRegistration<sweet::ShackExpIntegration>();
		shackTimesteppingSemiLagrangianSphereData = io_shackDict->getAutoRegistration<sweet::ShackTimesteppingSemiLagrangianSphereData>();

		shackPDESWETimeDisc = io_shackDict->getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
		shackPDESWEBenchmark = io_shackDict->getAutoRegistration<ShackPDESWESphereBenchmarks>();
		shackPDESWESphere = io_shackDict->getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	virtual bool setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
	)
	{
		timestepping_method = i_timestepping_method;
		ops = io_ops;
		return true;
	}


public:
	virtual void runTimestep(
			sweet::SphereData_Spectral &io_h_pert,
			sweet::SphereData_Spectral &io_u,
			sweet::SphereData_Spectral &io_v,

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;


	virtual bool implementsTimesteppingMethod(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string getIDString() = 0;

	virtual void printHelp()
	{
	}


	virtual void printImplementedTimesteppingMethods(
		std::ostream &o_ostream,
		const std::string &i_prefix
	)
	{
		o_ostream << i_prefix << "TODO" << std::endl;
	}

#if 0

/*
 * Martin@Joao Please let's discuss this.
 */
// needed for parareal (instead of using directly shackDict.disc.timestepping_method)
protected:
	std::string timestepping_method;
	int timestepping_order;
	int timestepping_order2;

public:

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void runTimestep(
			sweet::SphereData_Spectral &io_h,
			sweet::SphereData_Spectral &io_u,
			sweet::SphereData_Spectral &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	) = 0;

	virtual ~PDESWESphereTS_BaseInterface() {}

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE)
	void runTimestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		sweet::SphereData_Spectral h = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]);
		sweet::SphereData_Spectral u = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]);
		sweet::SphereData_Spectral v = *(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]);

		runTimestep(h, u, v,
				i_dt,
				i_sim_timestamp
			);

		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[0]) = h;
		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[1]) = u;
		*(io_data->get_pointer_to_data_SphereData_Spectral()->simfields[2]) = v;

	}



	// for parareal SL
	virtual void set_previous_solution(
				sweet::SphereData_Spectral &i_phi_prev,
				sweet::SphereData_Spectral &i_vrt_prev,
				sweet::SphereData_Spectral &i_div_prev
	)
	{
	};

	// for parareal SL
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		sweet::SphereData_Spectral phi_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[0];
		sweet::SphereData_Spectral vrt_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[1];
		sweet::SphereData_Spectral div_prev = *i_data->get_pointer_to_data_SphereData_Spectral()->simfields[2];

		set_previous_solution(phi_prev, vrt_prev, div_prev);
	};
#endif

#endif
};

#endif
