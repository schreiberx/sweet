/*
 * PDESWESphereTS_l_irk_n_erk.hpp
 *
 *  Created on: 30 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_IMEX_SDC_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_IMEX_SDC_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <vector>

#include <sweet/core/dict/DictArrayND.hpp>
#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_erk_n_erk.hpp"
#include "PDESWESphereTS_l_irk.hpp"

// Class to store references (pointers) to the initial solution for each time steps
class SWE_Variables_Pointers {
public:
	sweet::SphereData_Spectral* phi;
	sweet::SphereData_Spectral* vort;
	sweet::SphereData_Spectral* div;
	
	// Default constructor
	SWE_Variables_Pointers() : phi(nullptr), vort(nullptr), div(nullptr) {}

	// Set pointers to solution data
	void setPointer(
			sweet::SphereData_Spectral& phi,
			sweet::SphereData_Spectral& vort,
			sweet::SphereData_Spectral& div)
	{
		this->phi = &phi;
		this->vort = &vort;
		this->div = &div;
	}
};

/*
 * Class to store solution data at one node
 */
class SWE_VariableVector {
public:
	sweet::SphereData_Spectral phi;
	sweet::SphereData_Spectral vrt;
	sweet::SphereData_Spectral div;

public:
	// Only constructor
	SWE_VariableVector(const sweet::SphereData_Config* sphereDataConfig):
		phi(sphereDataConfig), 
		vrt(sphereDataConfig), 
		div(sphereDataConfig)
	{

	}

	// Empty constructor
	SWE_VariableVector()
	{
	}

	// Fill values for phi, vort and div
	SWE_VariableVector& operator=(const SWE_VariableVector& u) {
		phi = u.phi;
		vrt = u.vrt;
		div = u.div;
		return *this;
	}

	SWE_VariableVector& operator=(const SWE_Variables_Pointers& u) {
		phi = *(u.phi);
		vrt = *(u.vort);
		div = *(u.div);
		return *this;
	}

	bool setup(const sweet::SphereData_Config *sphere_data_config)
	{
		phi.setup(sphere_data_config);
		vrt.setup(sphere_data_config);
		div.setup(sphere_data_config);

		return true;
	}

	void swap(SWE_VariableVector &io_value)
	{
		phi.swap(io_value.phi);
		vrt.swap(io_value.vrt);
		div.swap(io_value.div);
	}
};


// Class to store all the solution data to each nodes and two iterations
class SDC_NodeStorage {
	std::vector<std::vector<SWE_VariableVector>> v;
public:
	SDC_NodeStorage(
			const sweet::SphereData_Config* sphereDataConfig,
			size_t num_nodes
	){
		std::vector<SWE_VariableVector> nodeValsK;
		std::vector<SWE_VariableVector> nodeValsK1;
		for (size_t i = 0; i < num_nodes; i++) {
			nodeValsK.push_back(SWE_VariableVector(sphereDataConfig));
			nodeValsK1.push_back(SWE_VariableVector(sphereDataConfig));
		}
		v.push_back(nodeValsK);
		v.push_back(nodeValsK1);
	}

	void swapIterate() {
		std::swap(v[0], v[1]);
//		std::rotate(v.begin(), v.begin()+1, v.end());
	}

	SWE_VariableVector& getK0(size_t i) {
		return v[0][i];
	}

	SWE_VariableVector& getK1(size_t i) {
		return v[1][i];
	}
};

// Class to store all the solution data to each nodes and two iterations
class SDC_NodeStorage_ {
	std::vector<SWE_VariableVector> data;


public:
	SDC_NodeStorage_()
	{
	}

public:
	SDC_NodeStorage_(
			const sweet::SphereData_Config* sphereDataConfig,
			size_t num_nodes
	){
		data.resize(num_nodes);

		for (size_t i = 0; i < num_nodes; i++)
		{
			data[i].setup(sphereDataConfig);
		}
	}

	SWE_VariableVector& operator[](int i)
	{
		return data[i];
	}

	void swap(SDC_NodeStorage_ &i_value)
	{
		for (size_t i = 0; i < data.size(); i++)
			data[i].swap(i_value.data[i]);
	}
};


class PDESWESphereTS_ln_imex_sdc	: public PDESWESphereTS_BaseInterface
{
public:
	bool setup_auto(
			const std::string &i_timestepping_method,
			sweet::SphereOperators *io_ops
		) override;

	bool setup(
			sweet::SphereOperators *io_ops,
			int i_order,	///< order of RK time stepping method for linear parts
			int i_order2,	///< order of RK time stepping method for non-linear parts
			int i_version_id
	);

public:
	bool implementsTimesteppingMethod(const std::string &i_timestepping_method) override;
	std::string getIDString() override;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	PDESWESphereTS_l_irk timestepping_l_irk;

	/*
	 * Non-linear time steppers
	 */
	PDESWESphereTS_l_erk_n_erk timestepping_l_erk_n_erk;

	int timestepping_order;


public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	) override
	{
		PDESWESphereTS_BaseInterface::shackRegistration(io_shackDict);

		timestepping_l_irk.shackRegistration(io_shackDict);
		timestepping_l_erk_n_erk.shackRegistration(io_shackDict);
		return true;
	}

private:
	/*
	 * SDC specific attributes
 	 */
	int nNodes;
	int nIter;

	std::string initialSweepType;  // Type of initial sweep

	bool diagonal;       // Whether or not using the diagonal implementation
	bool useEndUpdate;  // Whether or not use collocation update for end point

	typedef sweet::DictArrayND<1, double> Vec;
	typedef sweet::DictArrayND<2, double> Mat;

	Vec tau;
	Vec weights;
	Mat qMat;
	Mat qMatDeltaI;
	Mat qDeltaE;
	Mat qMatDelta0;

	/*
	 * Variables used as temporary storage locations during the time step
	 */

	// To store pointers to initial time steps
	SWE_VariableVector ts_u0;

	SDC_NodeStorage_ ts_linear_tendencies_k0;  		// linear term evaluations
	SDC_NodeStorage_ ts_nonlinear_tendencies_k0;	// non-linear term evaluations
	SDC_NodeStorage_ ts_linear_tendencies_k1;  		// linear term evaluations
	SDC_NodeStorage_ ts_nonlinear_tendencies_k1;	// non-linear term evaluations

	SWE_VariableVector ts_tmp_state;  // temporary variable


	// start of current time step
	double t0;

	// timestep size
	double dt;

	/*
	 * SDC specific methods
 	 */

	// Wrapper evaluating linear terms and storing them in separate variables (eval)
	void eval_linear(const SWE_VariableVector& u, SWE_VariableVector& eval, double t=-1);

	// Wrapper evaluating non-linear terms and storing them in separate variables (eval)
	void eval_nonlinear(const SWE_VariableVector& u, SWE_VariableVector& eval, double t=-1);

	/* Wrapper solving the implicit system built from the linear term :
	 u - dt*L(u) = rhs
	 LHS is always updated, independently on the dt value
	 WARNING : rhs variables are overwritten with u 
	*/ 
	void solveImplicit(SWE_VariableVector& rhs, double dt);

	// Perform y <- a*x + y
	void axpy(double a, const SWE_VariableVector& x, SWE_VariableVector& io_y);

	// Initialize nodes values
	void init_sweep();

	// Perform one sweep
	void sweep(size_t k);

	// Compute end-point solution and update step variables
	void computeEndPoint();

public:
	PDESWESphereTS_ln_imex_sdc();


	void runTimestep(
			sweet::SphereData_Spectral &io_phi_pert,
			sweet::SphereData_Spectral &io_vort,
			sweet::SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	) override;


	virtual ~PDESWESphereTS_ln_imex_sdc();
};

#endif // SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_IMEX_SDC_HPP_
