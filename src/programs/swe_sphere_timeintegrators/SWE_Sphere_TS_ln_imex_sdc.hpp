/*
 * SWE_Sphere_TS_l_irk_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

# pragma once

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include <vector>
#include <array>
#include <algorithm>
using std::vector;
using std::array;

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"

// Class to store a reference to the initial solution for each time steps
class SWE_Variables_Ref {
public:
	SphereData_Spectral* phi;
	SphereData_Spectral* vort;
	SphereData_Spectral* div;
	// Default constructor
	SWE_Variables_Ref() : phi(nullptr), vort(nullptr), div(nullptr) {}

	// Set pointers to solution data
	void setRef(SphereData_Spectral& phi, SphereData_Spectral& vort, SphereData_Spectral& div) {
		this->phi = &phi;
		this->vort = &vort;
		this->div = &div;
	}
	// Copy semantic OFF (just in case ...)
	SWE_Variables_Ref(const SWE_Variables_Ref& other) = delete;
	SWE_Variables_Ref& operator=(const SWE_Variables_Ref& other) = delete;
};

// Class to store solution data at one node
class SWE_Variables {
public:
	SphereData_Spectral phi;
	SphereData_Spectral vort;
	SphereData_Spectral div;
public:
	// Only constructor
	SWE_Variables(const SphereData_Config* sphereDataConfig): 
		phi(sphereDataConfig), 
		vort(sphereDataConfig), 
		div(sphereDataConfig) {}

	// Copy semantic OFF
	SWE_Variables(const SWE_Variables& other) = delete;
	SWE_Variables& operator=(const SWE_Variables& other) = delete;

	// Move semantic ON
	SWE_Variables(SWE_Variables&& other) = default;
	SWE_Variables& operator=(SWE_Variables&& other) = default;

	// Copy values for phi, vort and div
	void copyValues(SWE_Variables& u) {
		this->phi = u.phi;
		this->vort = u.vort;
		this->div = u.div;
	}
	void copyValues(SWE_Variables_Ref& u) {
		this->phi = *u.phi;
		this->vort = *u.vort;
		this->div = *u.div;
	}

};

// Class to store all the solution data to each nodes and two iterations
template<size_t nNodes>
class SDC_NodeStorage {
	vector<vector<SWE_Variables>> v;
public:
	SDC_NodeStorage(const SphereData_Config* sphereDataConfig){
		vector<SWE_Variables> nodeValsK;
		vector<SWE_Variables> nodeValsKp1;
		for (size_t i = 0; i < nNodes; i++) {
			nodeValsK.push_back(SWE_Variables(sphereDataConfig));
			nodeValsKp1.push_back(SWE_Variables(sphereDataConfig));
		}
		v.push_back(std::move(nodeValsK));
		v.push_back(std::move(nodeValsKp1));
	}

	void swapIterate() {
		std::rotate(v.begin(), v.begin()+1, v.end());
	}

	SWE_Variables& getK(size_t i) {
		return v[0][i];
	}

	SWE_Variables& getK1(size_t i) {
		return v[1][i];
	}
};


class SWE_Sphere_TS_ln_imex_sdc	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);
	std::string string_id();

	SimulationVariables& simVars;
	SphereOperators_SphereData& op;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_irk timestepping_l_irk;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_l_erk_n_erk timestepping_l_erk_n_erk;

	/*
	 * To be specified ...
 	 */
	int version_id;

	int timestepping_order;
	int timestepping_order2;

private:
	/*
	 * SDC specific attributes
 	 */
	const static size_t nNodes = 3;
	const static size_t nIter = 3;
	typedef array<double, nNodes> Vec;
	typedef array<Vec, nNodes> Mat;
	typedef array<array<double, nNodes>, 2> Storage;

	// Nodes, Quadrature matrix and QDelta approximation
	const Vec tau {0.0, 0.0, 0.0};
	const Mat qMat {{
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}}
	};
	const Mat qDeltaI {{
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}}
	};
	const Mat qDeltaE {{
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}}
	};

	// Variable storage required for SDC sweeps
	SDC_NodeStorage<nNodes> lTerms;  	// linear term evaluations
	SDC_NodeStorage<nNodes> nTerms;	// non-linear term evaluations
	SWE_Variables rhs;
	SphereData_Spectral tmp;  // for axpy

	// To store pointers to initial time steps
	SWE_Variables_Ref u0;
	double t0;
	double dt;

	/*
	 * SDC specific methods
 	 */

	// Wrapper evaluating linear terms and storing them in separate variables
	void evalLinearTerms(const SWE_Variables_Ref& u, SWE_Variables& eval, double t=-1);
	void evalLinearTerms(const SWE_Variables& u, SWE_Variables& eval, double t=-1);

	// Wrapper evaluating non-linear terms and storing them in separate variables
	void evalNonLinearTerms(const SWE_Variables_Ref& u, SWE_Variables& eval, double t=-1);
	void evalNonLinearTerms(const SWE_Variables& u, SWE_Variables& eval, double t=-1);

	/* Wrapper solving the implicit system built from the linear term :
	 u - dt*L(u) = rhs
	 LHS is always updated, independently on the dt value
	 WARNING : rhs variables are overwritten with u 
	*/ 
	void solveImplicit(SWE_Variables& rhs, double dt);

	// Perform y <- a*x + y
	void axpy(double a, SWE_Variables& x, SWE_Variables& y);

	// Initialize nodes values
	void initSweep();

	// Perform one sweep
	void sweep();

public:
	SWE_Sphere_TS_ln_imex_sdc(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method for linear parts
			int i_order2,	///< order of RK time stepping method for non-linear parts
			int i_version_id
	);

	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_ln_imex_sdc();
};
