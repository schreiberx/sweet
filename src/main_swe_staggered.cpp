
#include "DataArray.hpp"
#include "libgl/VisualizationEngine.hpp"
#include <unistd.h>

#define ENABLE_GUI	1


#if ENABLE_GUI
	#include "libgl/draw/GlDrawQuad.hpp"
	#include "libgl/draw/GlDrawCube.hpp"
	#include "libgl/core/CGlTexture.hpp"
	#include "libgl/shaders/shader_blinn/CShaderBlinn.hpp"
#endif


//// problem size
std::size_t N = 64;

// resolution
std::size_t res[2] = {N,N};

double h0 = 1000.0;

// gravitation
double g = 9.81;

// cfl condition
double CFL = 0.01;

// viscosity
double viscocity = 0.0;

// viscosity
double hyper_viscocity = 0.0;

// domain length
double domain_length = 1000;

// setup scenario
int setup_scenario = 0;


class SimulationSWEStaggered
#if ENABLE_GUI
		:
		public VisualizationEngine::ProgramCallbacks
#endif
{
public:
	double domain_size = 1000;

	// prognostics
	DataArray<2> P, u, v;

	// diagnostics
	DataArray<2> H, U, V;
	DataArray<2> eta;

	DataArray<2> P_t, u_t, v_t;

	DataArray<2> diff1_x, diff1_y;
	DataArray<2> diff2_x, diff2_y;

	DataArray<2> op_d_f_x, op_d_f_y;
	DataArray<2> op_d_b_x, op_d_b_y;
	DataArray<2> op_avg_f_x, op_avg_f_y;
	DataArray<2> op_avg_b_x, op_avg_b_y;

	// cell size
	double hx, hy;

	// number of simulated time steps
	int timestep_nr = 0;

	// mass
	double mass = 0;
	// energy
	double energy = 0;
	// potential enstropy
	double potential_entrophy = 0;

	/**
	 * See "The Dynamics of Finite-Difference Models of the Shallow-Water Equations", Robert Sadourny
	 *
	 * Prognostic:
	 *     V_t + \eta N x (P V) + grad(P + 0.5 V.V) = 0
	 *     P_t + div(P V) = 0
	 *
	 * Potential vorticity:
	 *     \eta = rot (V) / P
	 *
	 *   ______u______
	 *   |           |
	 *   |           |
	 *   v     P     v
	 *   |           |
	 *   |_____u_____|
	 */
public:
	SimulationSWEStaggered(
	)	:
		P(res),	// density/pressure
		u(res),	// velocity (x-direction)
		v(res),	// velocity (y-direction)

		H(res),	//
		U(res),	// mass flux (x-direction)
		V(res),	// mass flux (y-direction)
		eta(res),
		P_t(res),
		u_t(res),
		v_t(res),
		diff1_x(res),
		diff1_y(res),
		diff2_x(res),
		diff2_y(res),

		op_d_f_x(res),
		op_d_f_y(res),
		op_d_b_x(res),
		op_d_b_y(res),
		op_avg_f_x(res),
		op_avg_f_y(res),
		op_avg_b_x(res),
		op_avg_b_y(res)
	{
		CFL = ::CFL;

		hx = domain_size/(double)res[0];
		hy = domain_size/(double)res[1];

		{
			P.data_setall(h0);

			double center_x = 0.7;
			double center_y = 0.6;

			if (setup_scenario == 0)
			{
				/*
				 * radial dam break
				 */
				double radius = 0.2;
				for (std::size_t j = 0; j < res[1]; j++)
				{
					for (std::size_t i = 0; i < res[0]; i++)
					{
						double x = ((double)i+0.5)/(double)res[0];
						double y = ((double)j+0.5)/(double)res[1];

						double dx = x-center_x;
						double dy = y-center_y;

						if (radius*radius >= dx*dx+dy*dy)
							P.getDataRef(j,i) += 1.0;
					}
				}
			}

			if (setup_scenario == 1)
			{
				/*
				 * fun with Gaussian
				 */
				for (std::size_t j = 0; j < res[1]; j++)
				{
					for (std::size_t i = 0; i < res[0]; i++)
					{
						double x = ((double)i+0.5)/(double)res[0];
						double y = ((double)j+0.5)/(double)res[1];

						double dx = x-center_x;
						double dy = y-center_y;

						P.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
					}
				}
			}
		}


		u.data_setall(0);
		v.data_setall(0);


		double op_avg_f_x_kernel[3][3] = {
				{0,0,0},
				{0,1,1},
				{0,0,0},
		};
		op_avg_f_x.setup_kernel(op_avg_f_x_kernel, 0.5);

		double op_avg_f_y_kernel[3][3] = {
				{0,1,0},
				{0,1,0},
				{0,0,0},
		};
		op_avg_f_y.setup_kernel(op_avg_f_y_kernel, 0.5);

		double op_avg_b_x_kernel[3][3] = {
				{0,0,0},
				{1,1,0},
				{0,0,0},
		};
		op_avg_b_x.setup_kernel(op_avg_b_x_kernel, 0.5);

		double op_avg_b_y_kernel[3][3] = {
				{0,0,0},
				{0,1,0},
				{0,1,0},
		};
		op_avg_b_y.setup_kernel(op_avg_b_y_kernel, 0.5);



		double op_d_f_x_kernel[3][3] = {
				{0,0,0},
				{0,-1,1},
				{0,0,0}
		};
		op_d_f_x.setup_kernel(op_d_f_x_kernel, 1.0/hx);

		double op_d_f_y_kernel[3][3] = {
				{0,1,0},
				{0,-1,0},
				{0,0,0}
		};
		op_d_f_y.setup_kernel(op_d_f_y_kernel, 1.0/hy);


		double op_d_b_x_kernel[3][3] = {
				{0,0,0},
				{-1,1,0},
				{0,0,0}
		};
		op_d_b_x.setup_kernel(op_d_b_x_kernel, 1.0/hx);

		double op_d_b_y_kernel[3][3] = {
				{0,0,0},
				{0,1,0},
				{0,-1,0},
		};
		op_d_b_y.setup_kernel(op_d_b_y_kernel, 1.0/hy);


		double diff1_x_kernel[3][3] = {
				{0,0,0},
				{-1.0,0,1.0},
				{0,0,0}
		};
		diff1_x.setup_kernel(diff1_x_kernel, 1.0/(2.0*hx));

		double diff1_y_kernel[3][3] = {
				{0,-1.0,0},	// lower y coordinate
				{0,0,0},
				{0,1.0,0}	// higher y coordinate
		};
		diff1_y.setup_kernel(diff1_y_kernel, 1.0/(2.0*hy));

		double diff2_x_kernel[3][3] = {
				{0,0,0},
				{1.0,-2.0,1.0},
				{0,0,0}
			};
		diff2_x.setup_kernel(diff2_x_kernel, 1.0/(hx*hx));

		double diff2_y_kernel[3][3] = {
				{0,1.0,0},
				{0,-2.0,0},
				{0,1.0,0}
		};
		diff2_y.setup_kernel(diff2_y_kernel, 1.0/(hy*hy));
	}


	void run_timestep()
	{
		U = op_avg_f_x(P)*u;
		V = op_avg_f_y(P)*v;

		H = P + 0.5*(op_avg_b_x(u*u) + op_avg_b_y(v*v));

		eta = (op_d_f_x(v) - op_d_f_y(u)) / op_avg_f_x(op_avg_f_y(P));

		// mass
		mass = P.reduce_sum();
		// energy
		energy = 0.5*(P*P + P*op_avg_b_x(u*u) + P*op_avg_b_y(v*v)).reduce_sum();
		// potential enstropy
		potential_entrophy = 0.5*(eta*op_avg_f_x(op_avg_f_y(P))).reduce_sum();

		u_t = op_avg_b_y(eta*op_avg_f_x(V)) - op_d_f_x(H);
		v_t = op_avg_b_x(eta*op_avg_f_y(U)) - op_d_f_y(H);
		P_t = -op_d_b_x(U) - op_d_b_y(V);

#if 1
		if (viscocity > 0)
		{
			// TODO: is this correct?
			v_t -= (diff2_y(u) + diff2_y(v))*viscocity;
			u_t -= (diff2_x(u) + diff2_x(v))*viscocity;
		}

		if (hyper_viscocity > 0)
		{
			// TODO: is this correct?
			u_t -= (diff2_x(diff2_x(u)) + diff2_x(diff2_x(v)))*viscocity;
			v_t -= (diff2_y(diff2_y(u)) + diff2_y(diff2_y(v)))*viscocity;
		}
#endif

		double limit_speed = std::max(hx/u.reduce_maxAbs(), hy/v.reduce_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
 //           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
		double limit_gh = std::min(hx, hy)/std::sqrt(g*P.reduce_maxAbs());

//        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
        double dt = CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

		P += dt*P_t;
		u += dt*u_t;
		v += dt*v_t;

		timestep_nr++;
	}


#if ENABLE_GUI

	CGlDrawQuad *cGlDrawQuad;
	CGlTexture *cGlTexture;
	unsigned char *texture_data;

	VisualizationEngine *visualizationEngine;

	void vis_setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;

		cGlTexture = new CGlTexture(GL_TEXTURE_2D, GL_RGBA, GL_BLUE, GL_UNSIGNED_BYTE);
		cGlTexture->bind();
		cGlTexture->resize(res[1], res[0]);
		cGlTexture->unbind();

		texture_data = new unsigned char[P.array_data_cartesian_length];

		cGlDrawQuad = new CGlDrawQuad;

	}

	void vis_render()
	{
		// execute simulation time step
		run_timestep();

		P.requestDataInCartesianSpace();

		double scale_d = 1.0/(P-h0).reduce_maxAbs();
//		double scale_d = 0.5;

#pragma omp parallel for simd
		for (std::size_t i = 0; i < P.array_data_cartesian_length; i++)
		{
			double value;
			// average height
			value = P.array_data_cartesian_space[i]-h0;

			// scale
			value *= scale_d;

			// [-1;1] -> [0;255]
			value = (value+1.0)*0.5*255.0;

			texture_data[i] = value;
		}


		visualizationEngine->engineState->commonShaderPrograms.shaderTexturize.use();
		visualizationEngine->engineState->commonShaderPrograms.shaderTexturize.pvm_matrix_uniform.set(visualizationEngine->engineState->matrices.pvm);

			cGlTexture->bind();
			cGlTexture->setData(texture_data);

				cGlDrawQuad->render();

			cGlTexture->unbind();
		visualizationEngine->engineState->commonShaderPrograms.shaderBlinn.disable();
	}

	const char* vis_getStatusString()
	{
		static char title_string[1024];
		sprintf(title_string, "Timestep: %i, Mass: %f, Energy: %f, Potential Entrophy: %f", timestep_nr, mass, energy, potential_entrophy);
//		std::cout << "Timestep: " << timestep_nr << std::endl;
		return title_string;
	}

	void vis_viewportChanged(int i_width, int i_height)
	{

	}

	void vis_keypress(char i_key)
	{

	}

	void vis_shutdown()
	{
		delete [] texture_data;
		delete cGlTexture;
		delete cGlDrawQuad;

	}
#endif
};




int main(int i_argc, char *i_argv[])
{
	if (i_argc > 1)
	{
		res[0] = atoi(i_argv[1]);
		res[1] = res[0];
	}

	if (i_argc > 2)
		CFL = atof(i_argv[2]);

	if (i_argc > 3)
		viscocity = atof(i_argv[3]);

	if (i_argc > 4)
		hyper_viscocity = atof(i_argv[4]);

	if (i_argc > 5)
		setup_scenario = atoi(i_argv[5]);

	if (i_argc > 6)
		domain_length = atof(i_argv[6]);


	SimulationSWEStaggered *simulationSWE = new SimulationSWEStaggered;

	VisualizationEngine(simulationSWE, "SWE");

	delete simulationSWE;

	return 1;
}
