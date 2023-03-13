/*
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/PDESWESphere_BenchmarksCombined.cpp
 *
 * MULE_SCONS_OPTIONS: --parareal-sphere=enable
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 *
 *  Created on: 26 Jul 2022
 * Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */


#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/time/TimesteppingSemiLagrangianSphereData.hpp>
#include <sweet/core/time/ShackTimesteppingSemiLagrangianSphereData.hpp>


#include "../programs/pde_advectionSphere/PDEAdvectionSphereBenchmarksCombined.hpp"
#include "../programs/pde_advectionSphere/PDEAdvectionSphereTimeSteppers.hpp"
#include "../programs/pde_advectionSphere/time/ShackPDEAdvectionSphereTimeDiscretization.hpp"

#include "../programs/pde_advectionSphere/ProgramPDEAdvectionSphere.hpp"


#include "../programs/pde_sweSphere/PDESWESphere_BenchmarksCombined.hpp"
#include <sweet/parareal/Parareal_GenericData_SphereData_Spectral.hpp>


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
	sweet::ErrorBase error;

	sweet::SphereData_Config* sphereDataConfig;
	sweet::SphereOperators* ops;
	PDESWESphere_BenchmarksCombined sphereBenchmarks;
	sweet::Parareal_GenericData* data = nullptr;

public:
	Data()
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

	bool setup_1_shackRegister(sweet::ShackDictionary *io_shackDict)
	{
		sphereBenchmarks.setup_1_registerAllBenchmark();
		sphereBenchmarks.setup_2_shackRegistration(io_shackDict);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphereBenchmarks);
		return true;
	}

	bool setup_2_benchmarkDetection()
	{
		sphereBenchmarks.setup_3_benchmarkDetection();
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphereBenchmarks);
		return true;
	}

	bool setup_3_data(
			sweet::SphereData_Config* i_sphereDataConfig,
			sweet::SphereOperators* i_ops
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		ops = i_ops;

		data = new sweet::Parareal_GenericData_SphereData_Spectral<3>;
		data->setup_data_config(sphereDataConfig);
		data->allocate_data();
		sweet::SphereData_Spectral* phi_pert = data->get_pointer_to_data_SphereData_Spectral()->simfields[0];
		sweet::SphereData_Spectral* vrt = data->get_pointer_to_data_SphereData_Spectral()->simfields[1];
		sweet::SphereData_Spectral* div = data->get_pointer_to_data_SphereData_Spectral()->simfields[2];

		sphereBenchmarks.setup_4_benchmarkSetup_1_withoutOps();
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphereBenchmarks);

		sphereBenchmarks.setup_5_benchmarkSetup_2_withOps(ops);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphereBenchmarks);

		sphereBenchmarks.benchmark->getInitialState(*phi_pert, *vrt, *div);
		ERROR_FORWARD_ALWAYS_RETURN_BOOLEAN(sphereBenchmarks);

		return true;
	}

	void set_zero()
	{
		for (int i = 0; i < 3; i++)
			data->get_pointer_to_data_SphereData_Spectral()->simfields[i]->spectral_set_zero();
	}

	void restrict(Data* i_data)
	{
		data->restrict(*i_data->data);
	}

	void restrict(Data& i_data)
	{
		data->restrict(*i_data.data);
	}

	void pad_zeros(Data* i_data)
	{
		data->pad_zeros(*i_data->data);
	}
	void pad_zeros(Data& i_data)
	{
		data->pad_zeros(*i_data.data);
	}
};




int main(int i_argc, char *i_argv[])
{

	double eps = 1e-15;


	int N_Hs[5] = {16, 32, 64, 128, 256};
	int N_Ls[5] = {8, 16, 32, 64, 128};

	for (int i_H = 0; i_H < 5; i_H++)
	{
		for (int i_L = 0; i_L < 5; i_L++)
		{
			sweet::SphereData_Config sphereDataConfig_H;
			sweet::SphereData_Config sphereDataConfig_L;


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



			sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
			shackProgArgDict.setup();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

			sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);


			Data data_H;
			Data data_L;
			Data data_H_to_H;
			Data data_H_to_L;
			Data data_H_to_L_to_H;
			Data data_L_to_H;
			Data data_L_to_H_to_L;

			data_H.setup_1_shackRegister(&shackProgArgDict);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(data_H);
			data_L.setup_1_shackRegister(&shackProgArgDict);
			data_H_to_H.setup_1_shackRegister(&shackProgArgDict);
			data_H_to_L.setup_1_shackRegister(&shackProgArgDict);
			data_H_to_L_to_H.setup_1_shackRegister(&shackProgArgDict);
			data_L_to_H.setup_1_shackRegister(&shackProgArgDict);
			data_L_to_H_to_L.setup_1_shackRegister(&shackProgArgDict);


			shackProgArgDict.processProgramArguments();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);


			data_H.setup_2_benchmarkDetection();
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(data_H);
			data_L.setup_2_benchmarkDetection();
			data_H_to_H.setup_2_benchmarkDetection();
			data_H_to_L.setup_2_benchmarkDetection();
			data_H_to_L_to_H.setup_2_benchmarkDetection();
			data_L_to_H.setup_2_benchmarkDetection();
			data_L_to_H_to_L.setup_2_benchmarkDetection();



			shackProgArgDict.printShackData();


			sweet::ShackSphereDataOps shackSphereDataOps_H;
			shackSphereDataOps_H = *shackSphereDataOps;
			shackSphereDataOps_H.space_res_physical[0] = -1;
			shackSphereDataOps_H.space_res_physical[1] = -1;
			shackSphereDataOps_H.space_res_spectral[0] = N_H;
			shackSphereDataOps_H.space_res_spectral[1] = N_H;

			sweet::ShackSphereDataOps shackSphereDataOps_L;
			shackSphereDataOps_L = *shackSphereDataOps;
			shackSphereDataOps_L.space_res_physical[0] = -1;
			shackSphereDataOps_L.space_res_physical[1] = -1;
			shackSphereDataOps_L.space_res_spectral[0] = N_L;
			shackSphereDataOps_L.space_res_spectral[1] = N_L;

			sphereDataConfig_H.setupAuto(shackSphereDataOps_H);
			sphereDataConfig_L.setupAuto(shackSphereDataOps_L);

			sweet::SphereOperators ops_H(&sphereDataConfig_H, &shackSphereDataOps_H);
			sweet::SphereOperators ops_L(&sphereDataConfig_L, &shackSphereDataOps_L);

			data_H.setup_3_data(&sphereDataConfig_H, &ops_H);
			ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(data_H);
			data_L.setup_3_data(&sphereDataConfig_L, &ops_L);
			data_H_to_H.setup_3_data(&sphereDataConfig_H, &ops_H);
			data_H_to_L.setup_3_data(&sphereDataConfig_L, &ops_L);
			data_H_to_L_to_H.setup_3_data(&sphereDataConfig_H, &ops_H);
			data_L_to_H.setup_3_data(&sphereDataConfig_H, &ops_H);
			data_L_to_H_to_L.setup_3_data(&sphereDataConfig_L, &ops_L);

			// Error storage
			sweet::Parareal_GenericData* error_H = new sweet::Parareal_GenericData_SphereData_Spectral<3>;
			sweet::Parareal_GenericData* error_L = new sweet::Parareal_GenericData_SphereData_Spectral<3>;
			error_H->setup_data_config(&sphereDataConfig_H);
			error_L->setup_data_config(&sphereDataConfig_L);
			error_H->allocate_data();
			error_L->allocate_data();


			// Tests 1 and 2:
			data_H_to_H.set_zero();

			printTest("Test 1: high res -> high res (dummy restriction) ");
			data_H_to_H.restrict(data_H);
			*error_H = *data_H.data;
			*error_H -= *data_H_to_H.data;
			double error = error_H->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			data_H_to_H.set_zero();
			printTest("Test 2: high res -> high res (dummy prolongation) ");
			data_H_to_H.pad_zeros(data_H);
			*error_H = *data_H.data;
			*error_H -= *data_H_to_H.data;
			error = error_H->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			// Test 3:
			printTest("Test 3: high res -> low res (restriction) ");
			data_H_to_L.set_zero();
			data_H_to_L.restrict(data_H);
			*error_L = *data_L.data;
			*error_L -= *data_H_to_L.data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

			// Test 4:
			printTest("Test 4: high res -> low res -> high res (restriction + prolongation) ");
			data_H_to_L_to_H.set_zero();
			data_H_to_L_to_H.pad_zeros(data_H_to_L);
			*error_H = *data_H.data;
			*error_H -= *data_H_to_L_to_H.data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			assert(error < eps);

			// Test 5:
			printTest("Test 5: low res -> high res (prolongation) ");
			data_L_to_H.set_zero();
			data_L_to_H.restrict(data_L);
			*error_H = *data_H.data;
			*error_H -= *data_L_to_H.data;
			error = error_H->spectral_reduce_maxAbs(N_L - 1);
			printError(error, "(Up to mode " + std::to_string(N_L - 1) + ")" );
			printError(error_H->spectral_reduce_maxAbs(N_L), "(Up to mode " + std::to_string(N_L) + ")" );
			assert(error < eps);

			// Test 6:
			printTest("Test 6: low res -> high res -> low res (prolongation + restriction) ");
			data_L_to_H_to_L.set_zero();
			data_L_to_H_to_L.pad_zeros(data_L_to_H);
			*error_L = *data_L.data;
			*error_L -= *data_L_to_H_to_L.data;
			error = error_L->spectral_reduce_maxAbs();
			printError(error);
			assert(error < eps);

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
