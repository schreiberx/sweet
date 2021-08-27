/*
 * Author: Pedro Peixoto <ppeixoto@usp.br>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_BAROTROPIC_VORT_MODES_HPP_


#include "SWESphereBenchmarks_helpers.hpp"
#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereData_Config.hpp>

class SWESphereBenchmark_barotropic_vort_modes	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	SWESphereBenchmark_barotropic_vort_modes()
	{
	}

	// Mode setup
	std::size_t maxmodes;   //number of waves to be added
	static const int maxmodeslimit=100; //max number of waves
	std::size_t nmode[maxmodeslimit], mmode[maxmodeslimit];   // spherical harmonic indexes
	double ampl[maxmodeslimit]; //coefficients of normal modes 

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		bool found_bench = false;
		benchmark_name = i_benchmark_name;
		std::cout << i_benchmark_name << std::endl;
	    if (benchmark_name.find("barotropic_vort_modes") !=std::string::npos )
			found_bench = true;

		return
				benchmark_name == "barotropic_vort_modes"			||
				found_bench
		;
	}


	void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		simVars = i_simVars;
		ops = i_ops;
	}


	std::string get_help()
	{
		std::ostringstream stream;
		stream << "  BOROTROPIC VORTICITY MODES :" << std::endl;
		stream << "     'barotropic_vort_modes'" << std::endl;
		stream << "     'barotropic_vort_modes': standard test" << std::endl;
		stream << "     'barotropic_vort_modes_[N]_[n1]_[m1]_[v1]_...[nN]_[mN]_[vN]': mode choices" << std::endl;
		return stream.str();
	}

	/*
	 * Compute surface height for geostrophic balance with given velocities
	 *
	 * (Inspired by code of Jeffrey Whitaker)
	 */
	static void computeGeostrophicBalance_nonlinear(
			const SphereOperators_SphereData *i_ops,
			const SphereData_Spectral &i_vort,
			const SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi
	)
	{
		const SphereData_Config *sphereDataConfig = i_ops->sphereDataConfig;

		/*
		 * Compute velocities
		 */
		SphereData_Physical ug(sphereDataConfig);
		SphereData_Physical vg(sphereDataConfig);

		i_ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		std::cout<< "[MULE] benchmark_barotropic_vort_modes.umax: "<< ug.physical_reduce_max_abs() << std::endl;
		std::cout<< "[MULE] benchmark_barotropic_vort_modes.vmax: "<< vg.physical_reduce_max_abs() << std::endl;
		
		SphereData_Physical vrtg = i_vort.toPhys();

		SphereData_Physical tmpg1 = ug*(vrtg+i_ops->fg);
		SphereData_Physical tmpg2 = vg*(vrtg+i_ops->fg);

		SphereData_Spectral tmpspec1(sphereDataConfig);
		SphereData_Spectral tmpspec2(sphereDataConfig);

		i_ops->uv_to_vrtdiv(tmpg1, tmpg2, tmpspec1, tmpspec2);

		o_phi = i_ops->inv_laplace(tmpspec1) - 0.5*(ug*ug+vg*vg);
	}


	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{

		/*
		 * barotropic vorticity initialization:
		 * uses normal modes, that are spherical harmonic modes
		 *
		 * WARNING: It is designed for the barotropic vorticity equations, but can be used for SWE
		 * considering a balanced solution for the full non-linear equations
		 *           
		 */
		
		if (simVars->timecontrol.current_simulation_time == 0)
		{
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Overriding simulation parameters for this benchmark !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
			std::cout << "!!! WARNING: Amplitudes of normal modes will be divided by Earth Radius !!!" << std::endl;
			std::cout << "!!! WARNING !!!" << std::endl;
		}

		simVars->sim.sphere_rotating_coriolis_omega = 7.292e-5;
		simVars->sim.gravitation = 9.80616;
		simVars->sim.sphere_radius = 6.37122e6;
		simVars->sim.h0 = 29400.0/simVars->sim.gravitation;

		ops->setup(ops->sphereDataConfig, &(simVars->sim));
	

		double a = simVars->sim.sphere_radius;
		double omega = simVars->sim.sphere_rotating_coriolis_omega;
		//double u0 = 2.0*M_PI*a/(12.0*24.0*60.0*60.0);
		double alpha = 0;


		extract_bench_info(benchmark_name);
	

		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
		
		//SphereData_Physical ug, vg;

		// m=-M,...M are Fourier modes
		// n=|m|, ..., M are Legendre nodes
		// M max modes
		//                 n,m  (n>=m)
		// Set modes for the vorticity stream function
		//SphereData_Spectral psi(ops->sphereDataConfig); // = inv_laplace(i_vrt);
	
		//psi.spectral_set_zero();

		//Add mode values
		for (int n = 0; n < (int)maxmodes; n++){			
			if(nmode[n] < mmode[n]){
				std::cout<< "Modes: n="<<nmode[n]<<" , m="<<mmode[n]<<std::endl;
				SWEETError("SWESphereBenchmark_barotropic_vort_modes: n cannot be smaller than m");	
			}
			// Only real part of amplitude is considered, so ampl[] is a vector of real values
			o_vrt.spectral_set(nmode[n],mmode[n],ampl[n]/a);	
			//std::cout << nmode[n] <<", "<< mmode[n] <<", "<< ampl[n] << std::endl;
		}
		

		//psi.spectral_print();
		//psi.spectral_structure_print();

		//o_vrt = ops->laplace(psi); 

		//Set phi for steady state for SWE
		// TODO: Not sure if the initial condition is really steady state in all cases
		computeGeostrophicBalance_nonlinear(
				ops,
				o_vrt,
				o_div,
				o_phi_pert
		);
	}


	private:
	void extract_bench_info(const std::string &bcase)
	{
		
		std::string basic_name = "barotropic_vort_modes";
		std::string bcase_code = bcase;

		bcase_code.replace(benchmark_name.find(basic_name),basic_name.length()+1,"");

		if(bcase_code==""){
			SWEETError("SWESphereBenchmark_barotropic_vort_modes: please choose the normal mode case appending to benchmark name the code _[N]_[n1]_[m1]_[v1]_...[nN]_[mN]_[vN]");
		};

		std::cout<< "[MULE] benchmark_barotropic_vort_modes.case:"<< benchmark_name << std::endl;

		//Convert parameter to words
		std::string str = bcase_code;
		std::replace( str.begin(), str.end(), '_', ' ');

		std::stringstream iss(str);
		
		iss >> maxmodes;
		if(maxmodes>maxmodeslimit){
			std::cout<< "Modes:"<<maxmodes<<" , maxmodes hardcoded:"<<maxmodeslimit<<std::endl;
			SWEETError("SWESphereBenchmark_barotropic_vort_modes: Adjust maximun number of waves");	
		}
		std::cout<< "[MULE] benchmark_barotropic_vort_modes.code: "<<bcase_code<< std::endl;
		std::cout<< "[MULE] benchmark_barotropic_vort_modes.maxmodes: " << maxmodes << std::endl;
		
		//loop over waves
		for (int n = 0; n < (int)maxmodes; n++){
			//get a single mode
			iss >> nmode[n];
			iss >> mmode[n];
			iss >> ampl[n];
			std::cout<< "[MULE] benchmark_barotropic_vort_modes."<<n<<".nmode: "<< nmode[n] << std::endl;
			std::cout<< "[MULE] benchmark_barotropic_vort_modes."<<n<<".mmode: "<< mmode[n] << std::endl;
			std::cout<< "[MULE] benchmark_barotropic_vort_modes."<<n<<".ampl: "<< ampl[n] << std::endl;
		}
		
		return;

	}
	
#if 0
public:
	static
	void output_spectral_modes_evol(
			SimulationVariables &i_simVars, // Simulation variables
			PlaneData &i_mode_geo,    //Coeficients multiplying  mode
			PlaneData &i_mode_igwest,    //Coeficients multiplying  mode
			PlaneData &i_mode_igeast    //Coeficients multiplying  mode
	)
	{
		const char* filename_geo = "output_nm_geo_evol.txt";
		const char* filename_igwest = "output_nm_igwest_evol.txt";
		const char* filename_igeast = "output_nm_igeast_evol.txt";

		const PlaneDataConfig *planeDataConfig = i_mode_geo.planeDataConfig;

		std::ofstream file0;	
		std::ofstream file1;	
		std::ofstream file2;	
		if (i_simVars.timecontrol.current_timestep_nr == 0){
			file0.open(filename_geo, std::ofstream::out | std::ofstream::trunc);	
			file1.open(filename_igwest, std::ofstream::out | std::ofstream::trunc);	
			file2.open(filename_igeast, std::ofstream::out | std::ofstream::trunc);	
		}
		else
		{
			file0.open(filename_geo, std::ofstream::out | std::ofstream::app);	
			file1.open(filename_igwest, std::ofstream::out | std::ofstream::app);	
			file2.open(filename_igeast, std::ofstream::out | std::ofstream::app);	
		}

		std::stringstream buffer0, buffer1, buffer2;
		buffer0 << std::setprecision(8);
		buffer1 << std::setprecision(8);
		buffer2 << std::setprecision(8);

		//Headers
		if (i_simVars.timecontrol.current_timestep_nr == 0){
			//header
			buffer0 << "n\t time";
			buffer1 << "n\t time";
			buffer2 << "n\t time";
			T k0, k1;
			for (std::size_t ik1 = 0; ik1 < planeDataConfig->spectral_data_size[1]; ik1++)
			{
				for (std::size_t ik0 = 0; ik0 < planeDataConfig->spectral_data_size[0]; ik0++)
				{
					k0 = (T) ik0;
					if (ik1 <= planeDataConfig->spectral_data_size[1]/2)
						k1 = (T)ik1;
					else
						k1 = (T)((int)ik1-(int)planeDataConfig->spectral_data_size[1]);

					//Geo modes
					buffer0 << "\t("<< k0 <<","<<k1<<")";
					//West modes (adjust y_mode index)
					buffer1 << "\t("<< k0 <<","<<k1<<")";
					//East mode (adjust x_mode index)
					buffer2 << "\t("<< (k0 == 0 ? abs(k0) : -k0)  <<","<< (k1 == 0 ? abs(k1) : -k1) <<")";

				}
			}
			file0 << buffer0.str() << std::endl;
			file1 << buffer1.str() << std::endl;
			file2 << buffer2.str() << std::endl;
		}

		buffer0 = dump_normal_modes(i_simVars, i_mode_geo);
		file0 << buffer0.str() << std::endl;
		buffer0.str(std::string());

		buffer1 = dump_normal_modes(i_simVars, i_mode_igwest);
		file1 << buffer1.str() << std::endl;
		buffer1.str(std::string());

		buffer2 = dump_normal_modes(i_simVars, i_mode_igeast);
		file2 << buffer2.str() << std::endl;
		buffer2.str(std::string());

		file0.close();
		file1.close();
		file2.close();

		return ;
	}

public:
	static
	std::stringstream dump_normal_modes(
			SimulationVariables &i_simVars, // Simulation variables
			PlaneData &i_mode    //Coeficients multiplying  mode
	)
	{
		const PlaneDataConfig *planeDataConfig = i_mode.planeDataConfig;

		std::stringstream buffer;
		buffer << std::setprecision(8);

		buffer << i_simVars.timecontrol.current_timestep_nr;
		buffer << "\t" << i_simVars.timecontrol.current_simulation_time;
		const double zero=0.0;
		const double scale_factor =((double)(planeDataConfig->spectral_array_data_number_of_elements)); 
		//planeDataConfig->physical_data_size[0]*planeDataConfig->physical_data_size[1]));
		for (std::size_t ik1 = 0; ik1 < planeDataConfig->spectral_data_size[1]; ik1++)
		{
			for (std::size_t ik0 = 0; ik0 < planeDataConfig->spectral_data_size[0]; ik0++)
			{
				const std::complex<double> &value = i_mode.p_spectral_get(ik1, ik0);
				double norm = value.real()*value.real()+value.imag()*value.imag();
				norm=std::sqrt(norm/scale_factor);
				if (norm > 1.0e-15)
				{
					buffer << "\t" << norm;
				}
				else
				{
					buffer << "\t" << zero;
				}
			}
		}	
		return buffer;
	}

#endif
};

#endif
