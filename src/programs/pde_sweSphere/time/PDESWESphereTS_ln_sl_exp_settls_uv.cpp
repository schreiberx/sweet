/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_ln_sl_exp_settls_uv.hpp"


bool PDESWESphereTS_ln_sl_exp_settls_uv::implementsTimesteppingMethod(const std::string &i_timestepping_method)
{
	/*
	 * Should contain _exp and _settls as well as _uv to indicate vorticity-divergence formulation
	 */
	timestepping_method = i_timestepping_method;
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	return (
		(i_timestepping_method.find("_exp") != std::string::npos)		&&
		(i_timestepping_method.find("_settls") != std::string::npos)	&&
		(i_timestepping_method.find("_uv") != std::string::npos)		&&
		!(i_timestepping_method.find("_only") != std::string::npos)		&&
		true
	);

	return false;
}



bool PDESWESphereTS_ln_sl_exp_settls_uv::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	PDESWESphereTS_ln_sl_exp_settls_uv::LinearCoriolisTreatment_enum _linear_coriolis_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_IGNORE;
	PDESWESphereTS_ln_sl_exp_settls_uv::NLRemainderTreatment_enum _nonlinear_remainder_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::NL_REMAINDER_IGNORE;

	bool _original_linear_operator_sl_treatment = true;


	// Search for Coriolis
	if (timestepping_method.find("l_exp") != std::string::npos)
		_linear_coriolis_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_LINEAR;
	else if (timestepping_method.find("lc_na_sl") != std::string::npos)
		_linear_coriolis_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_SEMILAGRANGIAN;
	else if (timestepping_method.find("lc_") != std::string::npos)
		_linear_coriolis_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_NONLINEAR;


	// Search for nonlinear remainder term
	if (timestepping_method.find("_nr_") != std::string::npos)
		_nonlinear_remainder_treatment = PDESWESphereTS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR;

	if (timestepping_method.find("_ver1") != std::string::npos)
		_original_linear_operator_sl_treatment = false;
	else if (timestepping_method.find("_ver0") != std::string::npos)
		_original_linear_operator_sl_treatment = true;
	else
		_original_linear_operator_sl_treatment = true;

#if 1
	string_id_storage = "";

	if (_linear_coriolis_treatment == PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_LINEAR)
		string_id_storage += "l";
	else
		string_id_storage += "lg";

	string_id_storage += "_exp";

	if (_linear_coriolis_treatment == PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_SEMILAGRANGIAN)
		string_id_storage += "_lc";

	string_id_storage += "_na";

	string_id_storage += "_sl";

	if (_linear_coriolis_treatment == PDESWESphereTS_ln_sl_exp_settls_uv::CORIOLIS_NONLINEAR)
		string_id_storage += "_lc";

	if (_nonlinear_remainder_treatment == PDESWESphereTS_ln_sl_exp_settls_uv::NL_REMAINDER_NONLINEAR)
		string_id_storage += "_nr";

	string_id_storage += "_settls";

	if (!_original_linear_operator_sl_treatment)
		string_id_storage += "_ver1";

	std::string string_id_storage_ = string_id_storage + "_uv";

	if (timestepping_method != string_id_storage_)
	{
		if (!_original_linear_operator_sl_treatment)
		{
			std::cerr << "Provided time stepping method: "+timestepping_method << std::endl;
			std::cerr << "Detected time stepping method: "+string_id_storage_ << std::endl;
			SWEETError("Autodetection of parts of time stepping methods failed!");
		}

		std::string string_id_storage2 = string_id_storage+"_ver0"+"_uv";
		if (timestepping_method != string_id_storage2)
		{
			std::cerr << "Provided time stepping method: "+timestepping_method << std::endl;
			std::cerr << "Detected time stepping method (failed): "+string_id_storage_ << std::endl;
			std::cerr << "Detected alternative time stepping method (failed): "+string_id_storage2 << std::endl;
			SWEETError("Autodetection of parts of time stepping methods failed!");
		}
	}
#endif

	return setup_main(
			io_ops,
			timestepping_order,
			_linear_coriolis_treatment,					// Coriolis treatment
			_nonlinear_remainder_treatment,		// Nonlinear divergence treatment
			_original_linear_operator_sl_treatment			// original SL linear operator treatment
		);
}

bool PDESWESphereTS_ln_sl_exp_settls_uv::setup_main(
		sweet::SphereOperators *io_ops,
		int i_timestepping_order,
		LinearCoriolisTreatment_enum i_coriolis_treatment,
		NLRemainderTreatment_enum i_nonlinear_remainder_treatment,
		bool i_original_linear_operator_sl_treatment
)
{
	ops = io_ops;

	coriolis_treatment = i_coriolis_treatment;
	nonlinear_remainder_treatment = i_nonlinear_remainder_treatment;
	timestepping_order = i_timestepping_order;
	original_linear_operator_sl_treatment = i_original_linear_operator_sl_treatment;


	if (timestepping_order != 2)
		SWEETError("Invalid time stepping order, only 2nd order supported");

	// Setup semi-Lagrangian
	semiLagrangian.setup(ops->sphereDataConfig, shackTimesteppingSemiLagrangianSphereData, timestepping_order);

	swe_sphere_ts_ln_erk_split_uv.setup_main(ops, 1, true, true, true, true, false);

	swe_sphere_ts_l_rexi.setup_variant_50(
			io_ops,
			shackExpIntegration,
			"phi0",
			shackTimestepControl->current_timestep_size,
			shackPDESWESphere->sphere_use_fsphere,
			!(coriolis_treatment == CORIOLIS_LINEAR),
			timestepping_order
		);

	return true;
}

std::string PDESWESphereTS_ln_sl_exp_settls_uv::getIDString()
{
	return string_id_storage;
}


void PDESWESphereTS_ln_sl_exp_settls_uv::runTimestep(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt,	
		double i_simulation_timestamp
)
{
	if (timestepping_order == 2)
	{
		run_timestep_2nd_order(io_phi, io_vrt, io_div, i_fixed_dt, i_simulation_timestamp);
		return;
	}

	SWEETError("Only 2nd order TS method supported!");
}



void PDESWESphereTS_ln_sl_exp_settls_uv::run_timestep_2nd_order(
		sweet::SphereData_Spectral &io_U_phi,
		sweet::SphereData_Spectral &io_U_vrt,
		sweet::SphereData_Spectral &io_U_div,

		double i_dt,		
		double i_simulation_timestamp
)
{
	const sweet::SphereData_Spectral &U_phi = io_U_phi;
	const sweet::SphereData_Spectral &U_vrt = io_U_vrt;
	const sweet::SphereData_Spectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
#if !SWEET_PARAREAL
		/*
		 * First time step:
		 * Simply backup existing fields for multi-step parts of this algorithm.
		 */
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
#endif
	}

	/*
	 * Step 1) Compute departure points
	 *
	 * Step 2) Compute U_D at departure points
	 *
	 * Step 3) Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
	 *
	 * Step 4) Compute update U(t+dt) = exp(dt L)(U_D + dt * N*(t+0.5dt))
	 */

	/*
	 *************************************************************************************************
	 * Step 1) Compute departure points
	 *************************************************************************************************
	 */
	sweet::SphereData_Physical U_u_lon_prev, U_v_lat_prev;
	ops->vrtdiv_to_uv(U_vrt_prev, U_div_prev, U_u_lon_prev, U_v_lat_prev);

	sweet::SphereData_Physical U_u_lon, U_v_lat;
	ops->vrtdiv_to_uv(U_vrt, U_div, U_u_lon, U_v_lat);

	double dt_div_radius = i_dt / shackSphereDataOps->sphere_radius;

	// Calculate departure points
	sweet::ScalarDataArray pos_lon_d, pos_lat_d;
	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_lon_prev, dt_div_radius*U_v_lat_prev,
			dt_div_radius*U_u_lon, dt_div_radius*U_v_lat,
			pos_lon_d, pos_lat_d		// OUTPUT
		);


	/*
	 *************************************************************************************************
	 * Step 2) Compute U_D at departure points
	 *************************************************************************************************
	 */

	/*
	 * Compute X_D
	 */
	sweet::SphereData_Spectral U_phi_D, U_vrt_D, U_div_D;
	semiLagrangian.apply_sl_timeintegration_uv(
			ops,
			U_phi, U_vrt, U_div,
			pos_lon_d, pos_lat_d,
			U_phi_D, U_vrt_D, U_div_D
		);


	/*
	 *************************************************************************************************
	 * Step 3) Compute N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
	 *************************************************************************************************
	 */
	{
		/*
		 * N*(t+0.5dt) = 1/2 ([ 2*N(t) - exp(dtL) N(t-dt) ]_D + N(t))
		 */

		/*
		 * N(t-dt)
		 */
		const sweet::SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;
		sweet::SphereData_Spectral N_U_phi_prev_nr(sphereDataConfig, 0);
		sweet::SphereData_Spectral N_U_vrt_prev_nr(sphereDataConfig, 0);
		sweet::SphereData_Spectral N_U_div_prev_nr(sphereDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
			);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi_prev, U_vrt_prev, U_div_prev,
					N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
					i_simulation_timestamp-i_dt
				);
		}

		/*
		 * exp(dtL) N(t-dt)
		 */
		swe_sphere_ts_l_rexi.runTimestep(
			N_U_phi_prev_nr, N_U_vrt_prev_nr, N_U_div_prev_nr,
			i_dt,
			i_simulation_timestamp
		);


		/*
		 * N(t)
		 */
		sweet::SphereData_Spectral N_U_phi_nr(sphereDataConfig, 0);
		sweet::SphereData_Spectral N_U_vrt_nr(sphereDataConfig, 0);
		sweet::SphereData_Spectral N_U_div_nr(sphereDataConfig, 0);

		if (nonlinear_remainder_treatment == NL_REMAINDER_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv.euler_timestep_update_nr(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
				);
		}

		if (coriolis_treatment == CORIOLIS_NONLINEAR)
		{
			swe_sphere_ts_ln_erk_split_uv.euler_timestep_update_lc(
					U_phi, U_vrt, U_div,
					N_U_phi_nr, N_U_vrt_nr, N_U_div_nr,
					i_simulation_timestamp
				);
		}

		/*
		 * N(t+dt)_D = [ 2*N(t) - N(t-dt) ]_D
		 */
		sweet::SphereData_Spectral N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D;
		semiLagrangian.apply_sl_timeintegration_uv(
				ops,
				2.0 * N_U_phi_nr - N_U_phi_prev_nr,
				2.0 * N_U_vrt_nr - N_U_vrt_prev_nr,
				2.0 * N_U_div_nr - N_U_div_prev_nr,

				pos_lon_d, pos_lat_d,
				N_U_phi_next_D, N_U_vrt_next_D, N_U_div_next_D
			);

		/*
		 * Compute
		 * 		N*(t+0.5dt) = 1/2 ([ 2*N(t) - N(t-dt) ]_D + N(t))
		 * and
		 * 		U_D + dt * N*(t+0.5dt)
		 */
		U_phi_D += (i_dt * 0.5) * (N_U_phi_next_D + N_U_phi_nr);
		U_vrt_D += (i_dt * 0.5) * (N_U_vrt_next_D + N_U_vrt_nr);
		U_div_D += (i_dt * 0.5) * (N_U_div_next_D + N_U_div_nr);
	}


	/*
	 *************************************************************************************************
	 * Step 4) Compute update U(t+dt) = exp(dt L)(U_D + dt * N*(t+0.5dt))
	 *************************************************************************************************
	 */
	swe_sphere_ts_l_rexi.runTimestep(
		U_phi_D, U_vrt_D, U_div_D,
		i_dt,
		i_simulation_timestamp
	);

	/*
	 * Backup previous variables for multi-step SL method
	 */
	U_phi_prev = io_U_phi;
	U_vrt_prev = io_U_vrt;
	U_div_prev = io_U_div;

	/*
	 * Return new fields stored in R_*
	 */
	io_U_phi = U_phi_D;
	io_U_vrt = U_vrt_D;
	io_U_div = U_div_D;
}




void PDESWESphereTS_ln_sl_exp_settls_uv::printHelp()
{
	std::cout << "	Exponential SETTLS UV:" << std::endl;
	std::cout << "		+ lg_exp_na_sl_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ lg_exp_na_sl_nr_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ l_exp_na_sl_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ l_exp_na_sl_nr_settls_uv[_ver0]" << std::endl;
	std::cout << "		+ lg_exp_na_sl_settls_uv_ver1" << std::endl;
	std::cout << "		+ lg_exp_na_sl_nr_settls_uv_ver1" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_settls_uv_ver1" << std::endl;
	std::cout << "		+ lg_exp_na_sl_lc_nr_settls_uv_ver1" << std::endl;
	std::cout << "		+ l_exp_na_sl_settls_uv_ver1" << std::endl;
	std::cout << "		+ l_exp_na_sl_nr_settls_uv_ver1" << std::endl;
}



PDESWESphereTS_ln_sl_exp_settls_uv::PDESWESphereTS_ln_sl_exp_settls_uv()
{
}

PDESWESphereTS_ln_sl_exp_settls_uv::~PDESWESphereTS_ln_sl_exp_settls_uv()
{
}

