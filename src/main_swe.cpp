
#include "DataArray.hpp"
#include "libgl/VisualizationEngine.hpp"
#include "libgl/draw/GlDrawQuad.hpp"
#include "libgl/draw/GlDrawCube.hpp"
#include "libgl/core/CGlTexture.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinn.hpp"
#include <unistd.h>


#define DIM 2


// problem size
std::size_t N = 64;

double h0 = 1000.0;

// gravitation
double g = 9.81;

// cfl condition
double CFL = 0.01;

// viscosity
double viscocity = 0.0;

// viscosity
double hyper_viscocity = 0.0;

// setup scenario
int setup_scenario = 0;

class SimulationSWE
{
public:
	double domain_size = 1000;

	DataArray<2> h, u, v;
	DataArray<2> h_t, u_t, v_t;

	DataArray<2> diff1_x, diff1_y;
	DataArray<2> diff2_x, diff2_y;


	// cell size
	double hx, hy;


	// number of simulated time steps
	int timestep_nr = 0;


public:
	SimulationSWE(
			std::size_t i_res[2],
			double domain_size
	)	:
		h(i_res),
		u(i_res),
		v(i_res),
		h_t(i_res),
		u_t(i_res),
		v_t(i_res),
		diff1_x(i_res),
		diff1_y(i_res),
		diff2_x(i_res),
		diff2_y(i_res)
	{
		CFL = ::CFL;

		hx = domain_size/(double)i_res[0];
		hy = domain_size/(double)i_res[1];

		{
			h.data_setall(h0);

			double center_x = 0.7;
			double center_y = 0.6;

			if (setup_scenario == 0)
			{
				/*
				 * radial dam break
				 */
				double radius = 0.2;
				for (std::size_t j = 0; j < i_res[1]; j++)
				{
					for (std::size_t i = 0; i < i_res[0]; i++)
					{
						double x = ((double)i+0.5)/(double)i_res[0];
						double y = ((double)j+0.5)/(double)i_res[1];

						double dx = x-center_x;
						double dy = y-center_y;

						if (radius*radius >= dx*dx+dy*dy)
							h.getDataRef(j,i) += 1.0;
					}
				}
			}

			if (setup_scenario == 1)
			{
				/*
				 * fun with Gaussian
				 */
				for (std::size_t j = 0; j < i_res[1]; j++)
				{
					for (std::size_t i = 0; i < i_res[0]; i++)
					{
						double x = ((double)i+0.5)/(double)i_res[0];
						double y = ((double)j+0.5)/(double)i_res[1];

						double dx = x-center_x;
						double dy = y-center_y;

						h.getDataRef(j,i) += std::exp(-50.0*(dx*dx + dy*dy));
					}
				}
			}
		}

		if (0)
		{
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
		else
		{
			double diff1_x_kernel[5][5] = {
					{0,0,0,0,0},{0,0,0,0,0},
					{1.0, -8.0, 0, 8.0, -1.0},
					{0,0,0,0,0},{0,0,0,0,0}
			};
			diff1_x.setup_kernel(diff1_x_kernel, 1.0/(12*hx));

			double diff1_y_kernel[5][5] = {
					{0,0, -1.0, 0,0},
					{0,0,  8.0, 0,0},
					{0,0,  0.0, 0,0},
					{0,0, -8.0, 0,0},
					{0,0,  1.0, 0,0}
			};
			diff1_y.setup_kernel(diff1_y_kernel, 1.0/(12*hy));

			double diff2_x_kernel[5][5] = {
					{0,0,0,0,0},{0,0,0,0,0},
					{-1.0, 16.0, -30.0, 16.0, -1.0},
					{0,0,0,0,0},{0,0,0,0,0}
			};
			diff2_x.setup_kernel(diff2_x_kernel, 1.0/(12*hx));

			double diff2_y_kernel[5][5] = {
					{0,0,  -1.0, 0,0},
					{0,0,  16.0, 0,0},
					{0,0, -30.0, 0,0},
					{0,0,  16.0, 0,0},
					{0,0,  -1.0, 0,0}
			};
			diff2_y.setup_kernel(diff2_y_kernel, 1.0/12*hx);
		}

		u.data_setall(0);
		v.data_setall(0);
	}

	void run_timestep()
	{
		std::cout << "Timestep: " << timestep_nr << std::endl;

		/*
		 * non-conservative formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y
		 *	v_t = -g * h_y - u * v_x - v * v_y
		 */
		h_t = -diff1_x(u*h) - diff1_y(v*h);
		u_t = -g*diff1_x(h) - u*diff1_x(u) - v*diff1_y(u);
		v_t = -g*diff1_y(h) - u*diff1_x(v) - v*diff1_y(v);

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

		double limit_speed = std::max(hx/u.get_maxAbs(), hy/v.get_maxAbs());

        // limit by re
        double limit_visc = limit_speed;
//        if (viscocity > 0)
 //           limit_visc = (viscocity*0.5)*((hx*hy)*0.5);

        // limit by gravitational acceleration
		double limit_gh = std::min(hx, hy)/std::sqrt(g*h.get_maxAbs());

        std::cout << limit_speed << ", " << limit_visc << ", " << limit_gh << std::endl;
        double dt = CFL*std::min(std::min(limit_speed, limit_visc), limit_gh);

//       dt = 0.05;

		h += dt*h_t;
		u += dt*u_t;
		v += dt*v_t;

		timestep_nr++;
	}
};




VisualizationEngine *visualizationEngine;


class VisualizationSimulationSWE
{
public:
	static
	int &getGlobalArgCSingleton()
	{
		static int global_argc = 1;
		return global_argc;
	}

	static
	char **&getGlobalArgVSingleton()
	{
		static char **global_argv = nullptr;
		return global_argv;
	}


	SimulationSWE *simulationSWE;

	CGlDrawQuad *cGlDrawQuad;
	CGlTexture *cGlTexture;

	unsigned char *texture_data;

public:
	VisualizationSimulationSWE()
	{
		std::size_t i_res[2] = {N,N};
		double domain_size = 1024;

		simulationSWE = new SimulationSWE(i_res, domain_size);

		cGlTexture = new CGlTexture(GL_TEXTURE_2D, GL_RGBA, GL_RED, GL_UNSIGNED_BYTE);
		cGlTexture->bind();
		cGlTexture->resize(i_res[1], i_res[0]);
		cGlTexture->unbind();

		texture_data = new unsigned char[simulationSWE->h.array_data_cartesian_length];

		cGlDrawQuad = new CGlDrawQuad;
	}

	~VisualizationSimulationSWE()
	{
		delete [] texture_data;
		delete cGlTexture;
		delete cGlDrawQuad;
		delete simulationSWE;
	}

	void frame_callback()
	{
		simulationSWE->run_timestep();

		simulationSWE->h.requestDataInCartesianSpace();

		double scale_d = 1.0/(simulationSWE->h-h0).get_maxAbs();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < simulationSWE->h.array_data_cartesian_length; i++)
		{
			double value;
			// average height
			value = simulationSWE->h.array_data_cartesian_space[i]-h0;

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
};



VisualizationSimulationSWE *visSimulationSWE;

class ProgramCallbacks	:
		public VisualizationEngine::ProgramCallbacks
{
	void setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;
		visSimulationSWE = new VisualizationSimulationSWE;
	}

	void render()
	{
		visSimulationSWE->frame_callback();
	}

	const char* getStatusString()
	{
		return "";
	}

	void viewportChanged(int i_width, int i_height)
	{

	}

	void keypress(char i_key)
	{

	}

	void shutdown()
	{
	}

};


int main(int i_argc, char *i_argv[])
{
	if (i_argc > 1)
		N = atoi(i_argv[1]);

	if (i_argc > 2)
		CFL = atof(i_argv[2]);

	if (i_argc > 3)
		viscocity = atof(i_argv[3]);

	if (i_argc > 4)
		hyper_viscocity = atof(i_argv[4]);

	if (i_argc > 5)
		setup_scenario = atoi(i_argv[5]);

	ProgramCallbacks programCallbacks;

	VisualizationSimulationSWE::getGlobalArgCSingleton() = i_argc;
	VisualizationSimulationSWE::getGlobalArgVSingleton() = i_argv;

	VisualizationEngine(&programCallbacks, "SWE");

	delete visSimulationSWE;

	return 1;
}
