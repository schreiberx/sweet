/*
 * SWE_Plane_Normal_Modes.hpp
 *
 *  Created on: 17 Nov 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *      based on previous implementation by Martin Schreiber in swe_plane.cpp
 *
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_
#define SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <rexi/REXIFunctions.hpp>

#include <functional>
#if SWEET_EIGEN
#include <Eigen/Eigenvalues>
#endif
/**
 * SWE Plane normal mode
 */
class SWE_Plane_Normal_Modes
{
public:

#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif
	typedef std::complex<T> complex;

	static
	void add_normal_mode(
			std::size_t ik0,				//wavenumber in x
			std::size_t ik1,				// wavenumeber in y
			double geo_mode,    //Coeficient multiplying geostrophic mode
			double igwest_mode, //Coeficient multiplying west gravity mode
			double igeast_mode, //Coeficient multiplying east gravity mode
			PlaneData &io_h, // h: surface height (perturbation)
			PlaneData &io_u, // u: velocity in x-direction
			PlaneData &io_v, // v: velocity in y-direction
			SimulationVariables &i_simVars // Simulation variables
	)
	{

		REXIFunctions<T> rexiFunctions;

		const PlaneDataConfig *planeDataConfig = io_h.planeDataConfig;

		if (i_simVars.disc.space_grid_use_c_staggering)
			FatalError("Staggering not supported");
		
		//std::cout<<"Adding mode to fields"<<std::endl;

		io_h.request_data_spectral();
		io_u.request_data_spectral();
		io_v.request_data_spectral();

		//Check if k0 is in correct sprectral area
		//std::cout<< io_h.planeDataConfig->spectral_data_size[1] << std::endl;
		if( ik1<0 || ik1 >= planeDataConfig->spectral_data_size[1]) 
			FatalError("Normal_mode: mode not within reach");

		if( ik0<0 || ik0 >= planeDataConfig->spectral_data_size[0]) 
			FatalError("Normal_mode: mode not within reach");
	
		//Check for mirror effects
		T k0 = (T)ik0;
		T k1;
		if (ik1 < planeDataConfig->spectral_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = (T)((int)ik1-(int)planeDataConfig->spectral_data_size[1]);
		
		complex v[3][3];
		complex lambda[3];

		SWE_Plane_Normal_Modes::sw_eigen_decomp(
				k0,				//wavenumber in x
				k1,				// wavenumeber in y
				i_simVars,  // Input Simulation variables
				false, // Direct EV matrix (not inverse)
				v , // EV matrix
				lambda // output eigen values */
		);

		complex U[3];
		// Set normal mode acording to desired wave type
		// These are weights for the modes
		U[0] = geo_mode;
		U[1] = igwest_mode;
		U[2] = igeast_mode;


		//Define normal mode as combination of eigen vectors
		complex UEV[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				UEV[k] += v[k][j] * U[j];

		//std::cout<<"spectral before"<< std::endl;
		//io_v.print_spectralIndex();

		complex h_add, u_add, v_add;
		h_add = io_h.p_spectral_get(ik1, ik0)+UEV[0];
		u_add = io_u.p_spectral_get(ik1, ik0)+UEV[1];
		v_add = io_v.p_spectral_get(ik1, ik0)+UEV[2];

		/* Add normal mode to data */
		io_h.p_spectral_set(ik1, ik0, h_add);
		io_u.p_spectral_set(ik1, ik0, u_add);
		io_v.p_spectral_set(ik1, ik0, v_add);

		io_h.spectral_zeroAliasingModes();
		io_u.spectral_zeroAliasingModes();
		io_v.spectral_zeroAliasingModes();

		//Request physical data, to ensure that irroring is well well balanced (it may fill in the mirror mode)
		io_h.request_data_physical();
		io_u.request_data_physical();
		io_v.request_data_physical();

/*Debug output*/
#if 0

		std::cout<<"EV matrix"<<std::endl;
		for (int j = 0; j < 3; j++)	{
			for (int i = 0; i < 3; i++)
				std::cout<< v[j][i]<<" "; 
			std::cout<<std::endl;
		}
	
		std::cout<<"Eigen values"<<std::endl;
		for (int j = 0; j < 3; j++)
			std::cout<< lambda[j]<<" "; 
		std::cout<<std::endl;
		
		std::cout<<"Adding normal mode"<<std::endl;
		for (int j = 0; j < 3; j++)
			std::cout<< UEV[j]<<" "; 
		std::cout<<std::endl;
	
		std::cout<<ik0<<" "<<ik1<< " "<< io_v.p_spectral_get(ik1, ik0) << UEV[2] << std::endl;

#endif
		return;
	}

	static
	void project_to_normal_mode(
			std::size_t ik0,				//wavenumber in x
			std::size_t ik1,				// wavenumeber in y
			PlaneData &i_h, // h: surface height (perturbation)
			PlaneData &i_u, // u: velocity in x-direction
			PlaneData &i_v, // v: velocity in y-direction
			SimulationVariables &i_simVars, // Simulation variables
			complex &o_geo_mode,    //Output: Coeficient multiplying geostrophic mode
			complex &o_igwest_mode, //Output: Coeficient multiplying west gravity mode
			complex &o_igeast_mode //Output: Coeficient multiplying east gravity mode
	)
	{

		REXIFunctions<T> rexiFunctions;

		const PlaneDataConfig *planeDataConfig = i_h.planeDataConfig;

		if (i_simVars.disc.space_grid_use_c_staggering)
			FatalError("Staggering not supported");
		
		//std::cout<<"Adding mode to fields"<<std::endl;

		i_h.request_data_spectral();
		i_u.request_data_spectral();
		i_v.request_data_spectral();

		//Check if k0 is in correct sprectral area
		//std::cout<< io_h.planeDataConfig->spectral_data_size[1] << std::endl;
		if( ik1<0 || ik1 >= planeDataConfig->spectral_data_size[1]) 
			FatalError("Normal_mode: mode not within reach");

		if( ik0<0 || ik0 >= planeDataConfig->spectral_data_size[0]) 
			FatalError("Normal_mode: mode not within reach");
	
		//Check for mirror effects
		T k0 = (T)ik0;
		T k1;
		if (ik1 < planeDataConfig->spectral_data_size[1]/2)
			k1 = (T)ik1;
		else
			k1 = (T)((int)ik1-(int)planeDataConfig->spectral_data_size[1]);
		
		complex v[3][3];
		complex lambda[3];

		SWE_Plane_Normal_Modes::sw_eigen_decomp(
				k0,				//wavenumber in x
				k1,				// wavenumeber in y
				i_simVars,  // Input Simulation variables
				true, // Inverse ev matrix
				v , // inverse EV matrix
				lambda // output eigen values */
		);

		complex U[3];
		// Set (h,u,v) spectral coeficients
		// These are weights for the modes
		U[0] = i_h.p_spectral_get(ik1, ik0);
		U[1] = i_u.p_spectral_get(ik1, ik0);
		U[2] = i_v.p_spectral_get(ik1, ik0);


		//Apply inverse EV matrix to obtain data in EV space
		complex UEV[3] = {0.0, 0.0, 0.0};
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				UEV[k] += v[k][j] * U[j];

		//Return the modes
		o_geo_mode=UEV[0];
		o_igwest_mode=UEV[1];
		o_igeast_mode=UEV[2];
		//std::cout<< " Geost: "<< o_geo_mode<<std::endl;
		//std::cout<< " IGWest: "<< o_igwest_mode<<std::endl;
		//std::cout<< " IGEast: "<< o_igeast_mode<<std::endl;

		return;		
	}

	/* Get linear shallow water operator eigen decomposition */
	static
	void sw_eigen_decomp(
			T k0,				//wavenumber in x
			T k1,				// wavenumeber in y
			SimulationVariables &i_simVars, // Input Simulation variables
			bool i_inverse = false, // Input true, returns inverse matriz, false: returns direct
			complex o_v[3][3] = {0}, // output eigen vector (direct or inverse)
			complex o_evalues[3] =  0 // output eigen values (optional)
	)
	{
		REXIFunctions<T> rexiFunctions;
		bool i_evalues = false;
		if (o_evalues){
			i_evalues = true;
			//std::cout<<i_evalues<<" "<< o_evalues[0]   << std::endl;
		}
		else{
			i_evalues = false;
		}
		//std::cout<<o_evalues[1]<<std::endl;
		//std::cout<<i_inverse<<std::endl;

		if (i_simVars.disc.space_grid_use_c_staggering)
			FatalError("Staggering not supported");
		
		complex I(0.0, 1.0);
		//std::cout<<"Calculating EV for mode (" << k0 << ", " << k1 << ")" << std::endl;
		//std::cout<<"hi"<< std::endl;
		T s0 = i_simVars.sim.plane_domain_size[0];
		T s1 = i_simVars.sim.plane_domain_size[1];

		T f = i_simVars.sim.plane_rotating_f0;
		T h = i_simVars.sim.h0;
		T g = i_simVars.sim.gravitation;

		T sqrt_h = rexiFunctions.l_sqrt(h);
		T sqrt_g = rexiFunctions.l_sqrt(g);

		complex b = -k0*I;	// d/dx exp(I*k0*x) = I*k0 exp(I*k0*x)
		complex c = -k1*I;

		b = b*rexiFunctions.pi2/s0;
		c = c*rexiFunctions.pi2/s1;

		/*
		 * Matrix with Eigenvectors (column-wise)
		 */
		complex v[3][3];
		complex v_inv[3][3];

		/*
		 * Eigenvalues
		 */
		complex lambda[3];

		
		if (i_simVars.sim.plane_rotating_f0 == 0)
		{
			/*
			 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,0%7D,%7Bg*c,0,0%7D%7D
			 */
			if (k0 == 0 && k1 == 0)
			{
				v[0][0] = 1;
				v[1][0] = 0;
				v[2][0] = 0;

				v[0][1] = 0;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = 0;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = 0;
					lambda[2] = 0;
				}
			}
			else if (k0 == 0)
			{
				v[0][0] = 0;
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -sqrt_h/sqrt_g;
				v[1][1] = 0;
				v[2][1] = 1;

				v[0][2] = sqrt_h/sqrt_g;
				v[1][2] = 0;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -c*sqrt_g*sqrt_h;
					lambda[2] = c*sqrt_g*sqrt_h;;
				}
			}
			else if (k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,0%7D,%7Bg*c*0,0,0%7D%7D
				 */

				v[0][0] = 0;
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -sqrt_h/sqrt_g;
				v[1][1] = 1;
				v[2][1] = 0;

				v[0][2] = sqrt_h/sqrt_g;
				v[1][2] = 1;
				v[2][2] = 0;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -b*sqrt_g*sqrt_h;
					lambda[2] = b*sqrt_g*sqrt_h;
				}
			}
			else
			{
				v[0][0] = 0;
				v[1][0] = -c/b;
				v[2][0] = 1.0;

				v[0][1] = -(sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
				v[1][1] = b/c;
				v[2][1] = 1.0;

				v[0][2] = (sqrt_h*rexiFunctions.l_sqrtcplx(b*b + c*c))/(c*sqrt_g);
				v[1][2] = b/c;
				v[2][2] = 1.0;

				if (i_evalues){
					lambda[0] = 0.0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b + c*c)*sqrt_h*sqrt_g;
				}
			}
		}
		else
		{
			if (k0 == 0 && k1 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,0,0%7D,%7B0,0,f%7D,%7B0,-f,0%7D%7D
				 */
				v[0][0] = 0;
				v[1][0] = -I;
				v[2][0] = 1;

				v[0][1] = 0;
				v[1][1] = I;
				v[2][1] = 1;

				v[0][2] = 1;
				v[1][2] = 0;
				v[2][2] = 0;

				if (i_evalues){
					lambda[0] = I*f;
					lambda[1] = -I*f;
					lambda[2] = 0;
				}
			}
			else if (k0 == 0)
			{
				/*
				 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b*0,h*c%7D,%7Bg*b*0,0,f%7D,%7Bg*c,-f,0%7D%7D
				 */
				v[0][0] = f/(c*g);
				v[1][0] = 1;
				v[2][0] = 0;

				v[0][1] = -(c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][1] =  -f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][1] = 1;

				v[0][2] = (c*h)/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[1][2] = f/rexiFunctions.l_sqrtcplx(-f*f + c*c*g*h);
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(c*c*g*h-f*f);
				}
			}
			else if (k1 == 0)
			{
					/*
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c*0%7D,%7Bg*b,0,f%7D,%7Bg*c*0,-f,0%7D%7D
					 */
				v[0][0] = -f/(b*g);
				v[1][0] = 0;
				v[2][0] = 1;

				v[0][1] = -(b*h)/f;
				v[1][1] = rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
				v[2][1] = 1;

				v[0][2] = -(b*h)/f;
				v[1][2] = -rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h)/f;
				v[2][2] = 1;

				if (i_evalues){
					lambda[0] = 0;
					lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
					lambda[2] = rexiFunctions.l_sqrtcplx(b*b*g*h-f*f);
				}
			}
			else
				{
					/*
					 * Compute EV's of
					 * Linear operator
					 *
					 * [ 0  hb  hc ]
					 * [ gb  0   f ]
					 * [ gc -f   0 ]
					 *
					 * http://www.wolframalpha.com/input/?i=eigenvector%7B%7B0,h*b,h*c%7D,%7Bg*b,0,f%7D,%7Bg*c,-f,0%7D%7D
					 */

					v[0][0] = -f/(b*g);
					v[1][0] = -c/b;
					v[2][0] = 1.0;

					v[0][1] = -(c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][1] = -(f*f - b*b*g*h)/(b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][1] = 1.0;

					v[0][2] = -(-c*f*h + b*h*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h))/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[1][2] =  -(-f*f + b*b*g*h)/(-b*c*g*h + f*rexiFunctions.l_sqrtcplx(-f*f + b*b*g*h + c*c*g*h));
					v[2][2] = 1.0;

					if (i_evalues){
						lambda[0] = 0.0;
						lambda[1] = -rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
						lambda[2] =  rexiFunctions.l_sqrtcplx(b*b*g*h + c*c*g*h - f*f);
					}
				}
		}

			/*
			 * Invert Eigenvalue matrix
			 */

		if (i_inverse){
			v_inv[0][0] =  (v[1][1]*v[2][2] - v[1][2]*v[2][1]);
			v_inv[0][1] = -(v[0][1]*v[2][2] - v[0][2]*v[2][1]);
			v_inv[0][2] =  (v[0][1]*v[1][2] - v[0][2]*v[1][1]);

			v_inv[1][0] = -(v[1][0]*v[2][2] - v[1][2]*v[2][0]);
			v_inv[1][1] =  (v[0][0]*v[2][2] - v[0][2]*v[2][0]);
			v_inv[1][2] = -(v[0][0]*v[1][2] - v[0][2]*v[1][0]);

			v_inv[2][0] =  (v[1][0]*v[2][1] - v[1][1]*v[2][0]);
			v_inv[2][1] = -(v[0][0]*v[2][1] - v[0][1]*v[2][0]);
			v_inv[2][2] =  (v[0][0]*v[1][1] - v[0][1]*v[1][0]);

			complex s = v[0][0]*v_inv[0][0] + v[0][1]*v_inv[1][0] + v[0][2]*v_inv[2][0];

			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					v_inv[j][i] /= s;
			}
			//Return inverse matrix
			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					o_v[j][i] = v_inv[j][i] ;
			}
		}	
		else{
			//Return direct matrix
			for (int j = 0; j < 3; j++)	{
				for (int i = 0; i < 3; i++)
					o_v[j][i] = v[j][i] ;
			}
		}
		if (i_evalues){
			for (int j = 0; j < 3; j++)	{
				o_evalues[j] = lambda[j] ;
			}
		}
		return;
}


	template <typename TCallbackClass>
	static
	void normal_mode_analysis(
			PlaneData &io_prog_h_pert, // h: surface height (perturbation)
			PlaneData &io_prog_u, // u: velocity in x-direction
			PlaneData &io_prog_v, // v: velocity in y-direction
			int number_of_prognostic_variables,
			SimulationVariables &i_simVars, // Simulation variables
			TCallbackClass *i_class,
			void(TCallbackClass::* const i_run_timestep_method)(void)
	)
	{

		const PlaneDataConfig *planeDataConfig = io_prog_h_pert.planeDataConfig;

		// dummy time step to get time step size
		if (i_simVars.timecontrol.current_timestep_size <= 0)
			FatalError("Normal mode analysis requires setting fixed time step size");

		/*
		 *
		 * Mode-wise normal mode analysis
		 *
		 *
		 */

		if (i_simVars.misc.normal_mode_analysis_generation == 4)
		{
#if SWEET_EIGEN
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
			FatalError("SWE_Plane_Normal_Modes: This test was build for linear or linearized models, so please compile without dealising --plane-spectral-dealiasing=disable.");
#endif


			/*
			 * Setup all output files
			 */
			const char* filename; //general filename
			char buffer_real[1024];

			if (i_simVars.iodata.output_file_name == "")
				filename = "output_%s_t%020.8f.csv";
			else
				filename = i_simVars.iodata.output_file_name.c_str();

			sprintf(buffer_real, filename, "normal_modes_plane", i_simVars.timecontrol.current_timestep_size*i_simVars.iodata.output_time_scale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to files of the form '" << buffer_real << "'" << std::endl;

			//Positive inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_plane_igpos", i_simVars.timecontrol.current_timestep_size*i_simVars.iodata.output_time_scale);
			std::ofstream file_igpos(buffer_real, std::ios_base::trunc);

			//Negative inertia-gravity modes
			sprintf(buffer_real, filename, "normal_modes_plane_igneg", i_simVars.timecontrol.current_timestep_size*i_simVars.iodata.output_time_scale);
			std::ofstream file_igneg(buffer_real, std::ios_base::trunc);

			//Geostrophic modes
			sprintf(buffer_real, filename, "normal_modes_plane_geo", i_simVars.timecontrol.current_timestep_size*i_simVars.iodata.output_time_scale);
			std::ofstream file_geo(buffer_real, std::ios_base::trunc);

			//std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);
			file_igpos << std::setprecision(20);
			file_igneg << std::setprecision(20);
			file_geo << std::setprecision(20);

			file << "# dt " << i_simVars.timecontrol.current_timestep_size << std::endl;
			file << "# g " << i_simVars.sim.gravitation << std::endl;
			file << "# h " << i_simVars.sim.h0 << std::endl;
			file << "# r " << i_simVars.sim.sphere_radius << std::endl;
			file << "# f " << i_simVars.sim.plane_rotating_f0 << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
			file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

			file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
			file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
			file << "# normalmodegeneration " << i_simVars.misc.normal_mode_analysis_generation << std::endl;
			file << "# antialiasing ";
#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
			file << 1;
#else
			file << 0;
#endif
			file << std::endl;

			PlaneData* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			int number_of_prognostic_variables = 3;
			//The basic state is with zero in all variables
			// The only non zero variable in the basic state is the total height
			//    for which the constant is added within run_timestep()
			io_prog_h_pert.physical_set_zero();
			io_prog_u.physical_set_zero();
			io_prog_v.physical_set_zero();

			//int num_timesteps = 1;

			// Timestep and perturbation
			double dt = i_simVars.timecontrol.current_timestep_size;
			double eps = dt;

			//Matrix representing discrete linear operator in spectral space
			Eigen::MatrixXcf A(3,3) ;
			//Eigen solver
			Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
			//Final eigenvalues
			std::complex<double> eval[3];

			//For each spectral mode
			//for (int r = 0; r < 2; r++) //only required to get the symmetric half of the spectrum
			//{
			int r = 0;

			for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
			{
				std::cout << "." << std::flush;
				for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
				{
					//This is the mode to be analysed
					//std::cout << "Mode (i,j)= (" << i << " , " << j <<")"<< std::endl;


					for (int outer_prog_id = 0; outer_prog_id < number_of_prognostic_variables; outer_prog_id++)
					{

						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog[inner_prog_id]->spectral_set_zero();

						// activate mode via real coefficient
						prog[outer_prog_id]->p_spectral_set(j, i, 1.0);
						//Activate the symetric couterpart of the mode (only needed if j>0 )
						if (j > 0)
							prog[outer_prog_id]->p_spectral_set(planeDataConfig->spectral_data_size[1]-j, i, 1.0);

						/*
						 * RUN timestep
						 */
						prog[outer_prog_id]->request_data_physical();
						(i_class->*i_run_timestep_method)();

						/*
						 * compute
						 * 1/dt * (U(t+1) - U(t))
						 */
						prog[outer_prog_id]->request_data_spectral();

						std::complex<double> val = prog[outer_prog_id]->p_spectral_get(j, i);
						val = val - 1.0; //subtract U(0) from mode
						prog[outer_prog_id]->p_spectral_set(j, i, val);

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							(*prog[inner_prog_id]) /= eps;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							A(inner_prog_id,outer_prog_id)=prog[inner_prog_id]->p_spectral_get(j, i);;
						}

					}

					//std::cout << "Lik matrix" << std::endl;
					//std::cout << A << std::endl;

					//std::cout<<"Normal modes" << std::endl;
					ces.compute(A);
					for(int i=0; i<3; i++)
					{
						eval[i]=ces.eigenvalues()[i];
						//std::cout << "Eigenvalue "<< i << " : " << eval[i].real() <<" "<<eval[i].imag() << std::endl;

					}
					/* We will try to separate the modes in 3 types:
					 * -positive inertia-gravity (imag>f) - we will adopt coriolis f to test as if > zero, since the exact freq is sqrt(f^2+cK*K)
					 * -negative inertia-gravity (imag<-f)
					 * -negative inertia-gravity (imag aprox 0) - we will fit all other modes here
					 */
					int count_igpos=0;
					int count_igneg=0;
					int count_geo=0;
					for(int i=0; i<3; i++)
					{
						if(eval[i].imag() > 0.5 * i_simVars.sim.plane_rotating_f0)
						{
							//std::cout<< "IG pos mode: " << eval[i].imag() << std::endl;
							//file_igpos << eval[i].imag();
							file_igpos << eval[i].real()<< "\t" << eval[i].imag();
							file_igpos << "\t";
							count_igpos++;
						}
						if(eval[i].imag() < - 0.5 * i_simVars.sim.plane_rotating_f0)
						{
							//std::cout<< "IG neg mode: " << eval[i].imag() << std::endl;
							//file_igneg << eval[i].imag();
							file_igneg << eval[i].real()<< "\t" << eval[i].imag();
							file_igneg << "\t";
							count_igneg++;
						}
						if(eval[i].imag() >= - 0.5 * i_simVars.sim.plane_rotating_f0 && eval[i].imag() <=  0.5 * i_simVars.sim.plane_rotating_f0 )
						{
							//std::cout<< "IG geo mode: " << eval[i].imag() << std::endl;
							//file_geo << eval[i].imag();
							file_geo << eval[i].real()<< "\t" << eval[i].imag();
							file_geo << "\t";
							count_geo++;
						}
					}
					//Check if we got the correct modes
					if ( count_igpos * count_igneg * count_geo > 0 )
					{
						count_igpos=0;
						count_igneg=0;
						count_geo=0;
					}
					else
					{
						FatalError("SWE_Plane_Normal_Modes: Could not separate modes!!");
					}

					//std::cout<<"-------------------------" << std::endl;
				}
				file_igpos << std::endl;
				file_igneg << std::endl;
				file_geo << std::endl;
			}

			//}
			//std::cout<<"-------------------------" << std::endl;
			//FatalError("still needs work...");
#else
			FatalError("SWE_Plane_Normal_Modes: Cannot test this without Eigen library. Please compile with --eigen=enable");
#endif
		}
		/*
		 * Do a normal mode analysis using perturbation, see
		 * Hillary Weller, John Thuburn, Collin J. Cotter,
		 * "Computational Modes and Grid Imprinting on Five Quasi-Uniform Spherical C Grids"
		 */
		else
		{

			//run_timestep();
			const char* filename;
			char buffer_real[1024];

			if (i_simVars.iodata.output_file_name == "")
				filename = "output_%s_normalmodes.csv";
			else
				filename = i_simVars.iodata.output_file_name.c_str();


			sprintf(buffer_real, filename, "normal_modes_physical", i_simVars.timecontrol.current_timestep_size*i_simVars.iodata.output_time_scale);
			std::ofstream file(buffer_real, std::ios_base::trunc);
			std::cout << "Writing normal mode analysis to file '" << buffer_real << "'" << std::endl;

			std::cout << "WARNING: OUTPUT IS TRANSPOSED!" << std::endl;

			// use very high precision
			file << std::setprecision(20);

			PlaneData* prog[3] = {&io_prog_h_pert, &io_prog_u, &io_prog_v};

			/*
			 * Maximum number of prognostic variables
			 *
			 * Advection e.g. has only one
			 */
			if (number_of_prognostic_variables <= 0)
				FatalError("simVars.pde.number_of_prognostic_variables must be set!");

			if (number_of_prognostic_variables == 3)
			{
				io_prog_h_pert.physical_set_zero();
				io_prog_u.physical_set_zero();
				io_prog_v.physical_set_zero();
			}
			else if (number_of_prognostic_variables == 1)
			{
				io_prog_h_pert.physical_set_zero();
			}
			else
			{
				FatalError("Not yet supported");
			}

#if 0
			if (i_simVars.disc.timestepping_method == SimulationVariables::Discretization::LEAPFROG_EXPLICIT)
			{
				FatalError("Not yet tested and supported");
				std::cout << "WARNING: Leapfrog time stepping doesn't make real sense since 1st step is based on RK-like method" << std::endl;
				std::cout << "We'll do two Leapfrog time steps here to take the LF errors into account!" << std::endl;
				std::cout << "Therefore, we also halve the time step size here" << std::endl;

				i_simVars.timecontrol.current_timestep_size = 0.5*i_simVars.sim.CFL;
				i_simVars.sim.CFL = -i_simVars.timecontrol.current_timestep_size;
			}
#endif

			int num_timesteps = 1;
			if (i_simVars.misc.normal_mode_analysis_generation >= 10)
			{
				if (i_simVars.timecontrol.max_timesteps_nr > 0)
					num_timesteps = i_simVars.timecontrol.max_timesteps_nr;
			}

			if (i_simVars.timecontrol.max_simulation_time > 0)
				file << "# t " << i_simVars.timecontrol.max_simulation_time << std::endl;
			else
				file << "# t " << (num_timesteps*(-i_simVars.timecontrol.current_timestep_size)) << std::endl;

			file << "# g " << i_simVars.sim.gravitation << std::endl;
			file << "# h " << i_simVars.sim.h0 << std::endl;
//			file << "# r " << i_simVars.sim.sphere_radius << std::endl;
			file << "# f " << i_simVars.sim.plane_rotating_f0 << std::endl;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			int specmodes = planeDataConfig->get_spectral_iteration_range_area(0)+planeDataConfig->get_spectral_iteration_range_area(1);
			file << "# specnummodes " << specmodes << std::endl;
			file << "# specrealresx " << planeDataConfig->spectral_real_modes[0] << std::endl;
			file << "# specrealresy " << planeDataConfig->spectral_real_modes[1] << std::endl;
#endif

			file << "# physresx " << planeDataConfig->physical_res[0] << std::endl;
			file << "# physresy " << planeDataConfig->physical_res[1] << std::endl;
			file << "# normalmodegeneration " << i_simVars.misc.normal_mode_analysis_generation << std::endl;
			file << "# antialiasing ";

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
			file << 1;
#else
			file << 0;
#endif

			file << std::endl;


			// iterate over all prognostic variables
			for (int outer_prog_id = 0; outer_prog_id < number_of_prognostic_variables; outer_prog_id++)
			{
				if (i_simVars.misc.normal_mode_analysis_generation == 1 || i_simVars.misc.normal_mode_analysis_generation == 11)
				{
					// iterate over physical space
					for (std::size_t outer_i = 0; outer_i < planeDataConfig->physical_array_data_number_of_elements; outer_i++)
					{
						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

						std::cout << "." << std::flush;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog[inner_prog_id]->physical_set_zero();

						// activate mode
						prog[outer_prog_id]->request_data_physical();
						prog[outer_prog_id]->physical_space_data[outer_i] = 1;

						/*
						 * RUN timestep
						 */

						(i_class->*i_run_timestep_method)();

						if (i_simVars.misc.normal_mode_analysis_generation == 1)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog[outer_prog_id]->request_data_physical();
							prog[outer_prog_id]->physical_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								(*prog[inner_prog_id]) /= i_simVars.timecontrol.current_timestep_size;
						}

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							prog[inner_prog_id]->request_data_physical();
							for (std::size_t k = 0; k < planeDataConfig->physical_array_data_number_of_elements; k++)
							{
								file << prog[inner_prog_id]->physical_space_data[k];
								if (inner_prog_id != number_of_prognostic_variables-1 || k != planeDataConfig->physical_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
#if 1
				else if (i_simVars.misc.normal_mode_analysis_generation == 3 || i_simVars.misc.normal_mode_analysis_generation == 13)
				{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
					FatalError("Only available with if plane spectral space is activated during compile time!");
#else

					// iterate over spectral space
					for (int r = 0; r < 2; r++)
					{

						for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
						{
							for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
							{
								// reset time control
								i_simVars.timecontrol.current_timestep_nr = 0;
								i_simVars.timecontrol.current_simulation_time = 0;

								std::cout << "." << std::flush;

								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
									prog[inner_prog_id]->spectral_set_zero();

								// activate mode via real coefficient
								prog[outer_prog_id]->p_spectral_set(j, i, 1.0);

								/*
								 * RUN timestep
								 */
								(i_class->*i_run_timestep_method)();


								if (i_simVars.misc.normal_mode_analysis_generation == 3)
								{
									/*
									 * compute
									 * 1/dt * (U(t+1) - U(t))
									 */
									prog[outer_prog_id]->request_data_spectral();

									std::complex<double> val = prog[outer_prog_id]->p_spectral_get(j, i);
									val = val - 1.0;
									prog[outer_prog_id]->p_spectral_set(j, i, val);

									for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
										(*prog[inner_prog_id]) /= i_simVars.timecontrol.current_timestep_size;
								}


								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								{
									prog[inner_prog_id]->request_data_spectral();

									/*
									 * REAL
									 */

									for (int r = 0; r < 2; r++)
									{
										for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
										{
											for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
											{
												file << prog[inner_prog_id]->p_spectral_get(j, i).real();
												file << "\t";
											}
										}
									}
								}


								for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								{
									/*
									 * IMAG
									 */
									int c = 0;
									for (int r = 0; r < 2; r++)
									{
										for (std::size_t j = planeDataConfig->spectral_data_iteration_ranges[r][1][0]; j < planeDataConfig->spectral_data_iteration_ranges[r][1][1]; j++)
										{
											for (std::size_t i = planeDataConfig->spectral_data_iteration_ranges[r][0][0]; i < planeDataConfig->spectral_data_iteration_ranges[r][0][1]; i++)
											{
												file << prog[inner_prog_id]->p_spectral_get(j, i).imag();

												if (inner_prog_id != number_of_prognostic_variables-1 || c != specmodes-1)
													file << "\t";
												else
													file << std::endl;

												c++;
											}
										}
									}
								}
							}
						}
					}
#endif
				}
#else
				else if (i_simVars.misc.normal_mode_analysis_generation == 3 || i_simVars.misc.normal_mode_analysis_generation == 13)
				{
					PlaneDataComplex t1(planeDataConfig);
					PlaneDataComplex t2(planeDataConfig);
					PlaneDataComplex t3(planeDataConfig);
					PlaneDataComplex* prog_cplx[3] = {&t1, &t2, &t3};

					// iterate over spectral space
					for (std::size_t outer_i = 0; outer_i < planeDataConfig->spectral_complex_array_data_number_of_elements; outer_i++)
					{
						// reset time control
						i_simVars.timecontrol.current_timestep_nr = 0;
						i_simVars.timecontrol.current_simulation_time = 0;

						std::cout << "." << std::flush;

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
							prog_cplx[inner_prog_id]->spectral_set_zero();

						// activate mode via real coefficient
						prog_cplx[outer_prog_id]->request_data_spectral();
						prog_cplx[outer_prog_id]->spectral_space_data[outer_i].real(1);

						// convert PlaneDataComplex to PlaneData
						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							*prog[inner_prog_id] = Convert_PlaneDataComplex_To_PlaneData::physical_convert(*prog_cplx[inner_prog_id]);
							prog[inner_prog_id]->spectral_zeroAliasingModes();
						}


						/*
						 * RUN timestep
						 */
						(i_class->*i_run_timestep_method)();

						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							prog[inner_prog_id]->spectral_zeroAliasingModes();
#warning "update this physical_convert maybe to spectral_convert"

							*prog_cplx[inner_prog_id] = Convert_PlaneData_To_PlaneDataComplex::physical_convert(*prog[inner_prog_id]);

							prog_cplx[inner_prog_id]->request_data_spectral();
						}

						if (i_simVars.misc.normal_mode_analysis_generation == 3)
						{
							/*
							 * compute
							 * 1/dt * (U(t+1) - U(t))
							 */
							prog_cplx[outer_prog_id]->request_data_spectral();
							prog_cplx[outer_prog_id]->spectral_space_data[outer_i] -= 1.0;

							for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
								prog_cplx[inner_prog_id]->operator*=(1.0/i_simVars.timecontrol.current_timestep_size);
						}


						// convert PlaneDataComplex to PlaneData
						for (int inner_prog_id = 0; inner_prog_id < number_of_prognostic_variables; inner_prog_id++)
						{
							prog_cplx[inner_prog_id]->request_data_spectral();

							/*
							 * REAL
							 */
							for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
							{
								file << prog_cplx[inner_prog_id]->spectral_space_data[k].real();
								file << "\t";
							}

							/*
							 * IMAG
							 */
							for (std::size_t k = 0; k < planeDataConfig->spectral_complex_array_data_number_of_elements; k++)
							{
								file << prog_cplx[inner_prog_id]->spectral_space_data[k].imag();

								if (inner_prog_id != number_of_prognostic_variables-1 || k != planeDataConfig->spectral_complex_array_data_number_of_elements-1)
									file << "\t";
								else
									file << std::endl;
							}
						}
					}
				}
#endif
			}
		}
	}


	~SWE_Plane_Normal_Modes()
	{

	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_NORMAL_MODES_HPP_ */
