/*
 * test_sphere_pint_spatial_coarsening.cpp
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 *
 *  Created on: 26 Jul 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#include <iostream>
#include <cassert>
#include <sweet/core/SimulationVariables.hpp>
#include <sweet/core/MemBlockAlloc.hpp>
#include <sweet/core/sphere/SphereData_Config.hpp>
#include <sweet/core/SWEETError.hpp>
#include "../programs/swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <sweet/parareal/Parareal_GenericData_SphereData_Spectral.hpp>

SimulationVariables simVars;

sweet::SphereData_Config sphereDataConfigInstance_H;
sweet::SphereData_Config *sphereDataConfig_H = &sphereDataConfigInstance_H;

sweet::SphereData_Config sphereDataConfigInstance_L;
sweet::SphereData_Config *sphereDataConfig_L = &sphereDataConfigInstance_L;


void printTest(std::string test_name)
{
	std::cout << std::endl;
	std::cout << " ******************************************* " << std::endl;
	std::cout << " ****** " << test_name << std::endl;
	std::cout << " ******************************************* " << std::endl;
}

void printError(double val, std::string msg = "")
{
	std::cout << "    ----> ERROR " << msg << " : " << val << std::endl << std::endl;
}


class Data
{

public:
	SphereData_Config* sphereDataConfig;
	sweet::SphereOperators* op;
	BenchmarksSphereSWE sphereBenchmarks;
	Parareal_GenericData* data = nullptr;

public:
	Data(
		SphereData_Config* i_sphereDataConfig,
		sweet::SphereOperators* i_op
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
		sweet::SphereData_Spectral* phi_pert = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[0];
		sweet::SphereData_Spectral* vrt = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[1];
		sweet::SphereData_Spectral* div = this->data->get_pointer_to_data_SphereData_Spectral()->simfields[2];

		// To avoid warning message in benchmark
		simVars.timecontrol.current_simulation_time = 1;

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

	double eps = 1e-15;


	int N_Hs[5] = {16, 32, 64, 128, 256};
	int N_Ls[5] = {8, 16, 32, 64, 128};

	for (int i_H = 0; i_H < 5; i_H++)
	{
		for (int i_L = 0; i_L < 5; i_L++)
		{

			int N_H = N_Hs[i_H];
			int N_L = N_Ls[i_L];

			if (N_L >= N_H)
				continue;

			std::cout << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << " TESTING FOR N_H = " << N_H << "; " << "N_L = " << N_L << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;


			int N_physical[2] = {-1, -1};
			int N_spectral_H[2] = {N_H, N_H};
			int N_spectral_L[2] = {N_L, N_L};

			sphereDataConfigInstance_H.setupAuto(
								N_physical,
								N_spectral_H,
								simVars.misc.reuse_spectral_transformation_plans,
								simVars.misc.verbosity,
								simVars.parallelization.num_threads_space
							);

			sphereDataConfigInstance_L.setupAuto(
								N_physical,
								N_spectral_L,
								simVars.misc.reuse_spectral_transformation_plans,
								simVars.misc.verbosity,
								simVars.parallelization.num_threads_space
							);

			sweet::SphereOperators op_H(sphereDataConfig_H, &(simVars.sim));
			sweet::SphereOperators op_L(sphereDataConfig_L, &(simVars.sim));

			Data* data_H = new Data(sphereDataConfig_H, &op_H);
			data_H->setup();

			Data* data_L = new Data(sphereDataConfig_L, &op_L);
			data_L->setup();

			// Error storage
			Parareal_GenericData* error_H = new Parareal_GenericData_SphereData_Spectral<3>;
			Parareal_GenericData* error_L = new Parareal_GenericData_SphereData_Spectral<3>;
			error_H->setup_data_config(sphereDataConfig_H);
			error_L->setup_data_config(sphereDataConfig_L);
			error_H->allocate_data();
			error_L->allocate_data();


			// Tests 1 and 2:
			Data* data_H_to_H = new Data(sphereDataConfig_H, &op_H);
			data_H_to_H->setup();
			data_H_to_H->set_zero();

			printTest("Test 1: high res -> high res (dummy restriction) ");
			data_H_to_H->restrict(data_H);
			*error_H = *data_H->data;
			*error_H -= *data_H_to_H->data;
			double error = error_H->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			data_H_to_H->set_zero();
			printTest("Test 2: high res -> high res (dummy prolongation) ");
			data_H_to_H->pad_zeros(data_H);
			*error_H = *data_H->data;
			*error_H -= *data_H_to_H->data;
			error = error_H->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			delete data_H_to_H;


			// Test 3:
			Data* data_H_to_L = new Data(sphereDataConfig_L, &op_L);
			data_H_to_L->setup();
			printTest("Test 3: high res -> low res (restriction) ");
			data_H_to_L->set_zero();
			data_H_to_L->restrict(data_H);
			*error_L = *data_L->data;
			*error_L -= *data_H_to_L->data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			// Test 4:
			Data* data_H_to_L_to_H = new Data(sphereDataConfig_H, &op_H);
			data_H_to_L_to_H->setup();
			printTest("Test 4: high res -> low res -> high res (restriction + prolongation) ");
			data_H_to_L_to_H->set_zero();
			data_H_to_L_to_H->pad_zeros(data_H_to_L);
			*error_H = *data_H->data;
			*error_H -= *data_H_to_L_to_H->data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			assert(error < eps);
			delete data_H_to_L;
			delete data_H_to_L_to_H;

			// Test 5:
			Data* data_L_to_H = new Data(sphereDataConfig_H, &op_H);
			data_L_to_H->setup();
			printTest("Test 5: low res -> high res (prolongation) ");
			data_L_to_H->set_zero();
			data_L_to_H->restrict(data_L);
			*error_H = *data_H->data;
			*error_H -= *data_L_to_H->data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			assert(error < eps);

			// Test 6:
			Data* data_L_to_H_to_L = new Data(sphereDataConfig_L, &op_L);
			data_L_to_H_to_L->setup();
			printTest("Test 6: low res -> high res -> low res (prolongation + restriction) ");
			data_L_to_H_to_L->set_zero();
			data_L_to_H_to_L->pad_zeros(data_L_to_H);
			*error_L = *data_L->data;
			*error_L -= *data_L_to_H_to_L->data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);
			delete data_L_to_H;
			delete data_L_to_H_to_L;



			delete data_H;
			delete data_L;

			delete error_H;
			delete error_L;

			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << " TEST SUCCESSFUL" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << std::endl;
		}
	}

	std::cout << std::endl;
	std::cout << " !!!!!!!!!!!!!!!!!! " << std::endl;
	std::cout << "  ALL TESTS PASSED " << std::endl;
	std::cout << " !!!!!!!!!!!!!!!!!! " << std::endl;
	std::cout << std::endl;
}
