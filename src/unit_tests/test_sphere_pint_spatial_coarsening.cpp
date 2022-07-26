/*
 * test_sphere_pint_spatial_coarsening.cpp
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 *
 *  Created on: 26 Jul 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#include <iostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/MemBlockAlloc.hpp>
#include <sweet/sphere/SphereData_Config.hpp>
#include <sweet/SWEETError.hpp>
#include "../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <parareal/Parareal_GenericData_SphereData_Spectral.hpp>

SimulationVariables simVars;

SphereData_Config sphereDataConfigInstance_highres;
SphereData_Config *sphereDataConfig_highres = &sphereDataConfigInstance_highres;

SphereData_Config sphereDataConfigInstance_lowres;
SphereData_Config *sphereDataConfig_lowres = &sphereDataConfigInstance_lowres;


void printTest(std::string test_name)
{
	std::cout << std::endl;
	std::cout << " ******************************************* " << std::endl;
	std::cout << " ****** " << test_name << std::endl;
	std::cout << " ******************************************* " << std::endl;
}

void printError(double val, std::string msg = "")
{
	std::cout << "    ----> ERROR "<< msg << " : " << val << std::endl << std::endl;
}


class Data
{

public:
	SphereData_Config* sphereDataConfig;
	SphereOperators_SphereData* op;
	BenchmarksSphereSWE sphereBenchmarks;
	Parareal_GenericData* data = nullptr;

public:
	Data(
		SphereData_Config* i_sphereDataConfig,
		SphereOperators_SphereData* i_op
		)
		:
			sphereDataConfig(i_sphereDataConfig),
			op(i_op)
	{
	}

	~Data()
	{
		if (data)
		{
			delete data;
			data = nullptr;
		}
	}

	void setup()
	{

		this->data = new Parareal_GenericData_SphereData_Spectral<3>;
		this->data->setup_data_config(this->sphereDataConfig);
		this->data->allocate_data();
		SphereData_Spectral* phi_pert = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[0];
		SphereData_Spectral* vrt = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[1];
		SphereData_Spectral* div = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[2];

		this->sphereBenchmarks.setup(simVars, *this->op);
		this->sphereBenchmarks.master->get_initial_state(*phi_pert, *vrt, *div);

	}

	void set_zero()
	{
		for (int i = 0; i < 3; i++)
			this->data->get_pointer_to_data_SphereData_Spectral()->simfields[i]->spectral_set_zero();
	}

	void restrict(Data* i_data)
	{
		this->data->restrict(*i_data->data);
	}

	void pad_zeros(Data* i_data)
	{
		this->data->pad_zeros(*i_data->data);
	}
};



int main(
		int i_argc,
		char *const i_argv[]
)
{

	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		std::cout << std::endl;
		return -1;
	}

	int N_highres = 256;
	int N_lowres = 32;

	int N_physical[2] = {-1, -1};
	int N_spectral_highres[2] = {N_highres, N_highres};
	int N_spectral_lowres[2] = {N_lowres, N_lowres};

	sphereDataConfigInstance_highres.setupAuto(
						N_physical,
						N_spectral_highres,
						simVars.misc.reuse_spectral_transformation_plans
					);

	sphereDataConfigInstance_lowres.setupAuto(
						N_physical,
						N_spectral_lowres,
						simVars.misc.reuse_spectral_transformation_plans
					);

	SphereOperators_SphereData op_highres(sphereDataConfig_highres, &(simVars.sim));
	SphereOperators_SphereData op_lowres(sphereDataConfig_lowres, &(simVars.sim));

	Data* data_highres = new Data(sphereDataConfig_highres, &op_highres);
	data_highres->setup();

	Data* data_lowres = new Data(sphereDataConfig_lowres, &op_lowres);
	data_lowres->setup();

	// Error storage
	Parareal_GenericData* error_highres = new Parareal_GenericData_SphereData_Spectral<3>;
	Parareal_GenericData* error_lowres = new Parareal_GenericData_SphereData_Spectral<3>;
	error_highres->setup_data_config(sphereDataConfig_highres);
	error_lowres->setup_data_config(sphereDataConfig_lowres);
	error_highres->allocate_data();
	error_lowres->allocate_data();




	// Tests 1 and 2: High res -> high res (no restriction nor padding zeros)
	Data* data_highres_to_highres = new Data(sphereDataConfig_highres, &op_highres);
	data_highres_to_highres->setup();
	data_highres_to_highres->set_zero();

	printTest("Test 1: high res -> high res (dummy restriction) ");
	data_highres_to_highres->restrict(data_highres);
	*error_highres = *data_highres->data;
	*error_highres -= *data_highres_to_highres->data;
	printError(error_highres->spectral_reduce_maxAbs());

	data_highres_to_highres->set_zero();
	printTest("Test 2: high res -> high res (dummy prolongation) ");
	data_highres_to_highres->pad_zeros(data_highres);
	*error_highres = *data_highres->data;
	*error_highres -= *data_highres_to_highres->data;
	printError(error_highres->spectral_reduce_maxAbs());

	delete data_highres_to_highres;


	// Test 3:
	Data* data_highres_to_lowres = new Data(sphereDataConfig_lowres, &op_lowres);
	data_highres_to_lowres->setup();
	printTest("Test 3: high res -> low res (restriction) ");
	data_highres_to_lowres->set_zero();
	data_highres_to_lowres->restrict(data_highres);
	*error_lowres = *data_lowres->data;
	*error_lowres -= *data_highres_to_lowres->data;
	printError(error_lowres->spectral_reduce_maxAbs());

	// Test 4: 
	Data* data_highres_to_lowres_to_highres = new Data(sphereDataConfig_highres, &op_highres);
	data_highres_to_lowres_to_highres->setup();
	printTest("Test 3: high res -> low res -> high res (restriction + prolongation) ");
	data_highres_to_lowres_to_highres->set_zero();
	data_highres_to_lowres_to_highres->pad_zeros(data_highres_to_lowres);
	*error_highres = *data_highres->data;
	*error_highres -= *data_highres_to_lowres_to_highres->data;
	printError(error_highres->spectral_reduce_maxAbs(N_lowres - 1),"(Up to mode " + std::to_string(N_lowres - 1) + ")" );
	printError(error_highres->spectral_reduce_maxAbs(N_lowres), "(Up to mode " + std::to_string(N_lowres) + ")" );
	delete data_highres_to_lowres;
	delete data_highres_to_lowres_to_highres;

	// Test 5: 




	delete data_highres;
	delete data_lowres;

	delete error_highres;
	delete error_lowres;
}
