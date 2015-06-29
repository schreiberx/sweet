
#include "DataArray.hpp"
#include "libgl/VisualizationEngine.hpp"
#include "libgl/draw/GlDrawQuad.hpp"
#include "libgl/draw/GlDrawCube.hpp"
#include "libgl/core/CGlTexture.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinn.hpp"
       #include <unistd.h>


#define DIM 2


void runTests(
	std::size_t size_y,
	std::size_t size_x
)
{
	std::size_t size[2] = {size_x,size_y};

	DataArray<2> h(size);
	DataArray<2> hu(size);
	DataArray<2> hv(size);

	h.data_setall(10);
	hu.data_setall(0);
	hv.data_setall(0);

	int c = 1;
	for (std::size_t j = 0; j < size[1]; j++)
		for (std::size_t i = 0; i < size[0]; i++)
		{
			h.getDataRef(j, i) = c;
			c++;
		}

	// shift left test
	if (1)
	{
		DataArray<2> h_t(size);
		DataArray<2> u_t(size);
		DataArray<2> v_t(size);

		DataArray<2> op_shift_left(size);
		double shift_left_kernel[3][3] = {{0,0,0},{0,0,1},{0,0,0}};
		op_shift_left.setup_kernel(shift_left_kernel);

		std::cout << "H 0" << std::endl;
		std::cout << h << std::endl;

		h_t = op_shift_left(h);
		std::cout << "H 1" << std::endl;
		std::cout << h_t << std::endl;

		std::cout << h_t << std::endl;
		u_t = op_shift_left(h_t);
		std::cout << "H 2" << std::endl;
		std::cout << u_t << std::endl;

		h = op_shift_left(h);
		h = op_shift_left(h);
		std::cout << "H 2" << std::endl;
		std::cout << h << std::endl;

		return;
	}

	// shift up
	if (1)
	{
		DataArray<2> op_shift_left(size);
		double shift_left_kernel[3][3] = {{0,1,0},{0,0,0},{0,0,0}};
		op_shift_left.setup_kernel(shift_left_kernel);
		DataArray<2> shifted_result = op_shift_left(h);

		std::cout << "H" << std::endl;
		std::cout << h << std::endl;
		std::cout << "SHIFT UP" << std::endl;
		std::cout << shifted_result << std::endl;
	}

#if 0
	DataArray<2> test(size);
	double test_kernel[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
	//double test_kernel[5][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15},{16,17,18,19,20},{21,22,23,24,25}};
	test.setup_kernel(test_kernel);
	std::cout << test << std::endl;
	test = h;
	test.test_fftw();
	std::cout << test << std::endl;
	return;
#endif
}


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

	// gravitation
	double g = 9.81;

	// cfl condition
	double CFL = 0.01;

	// viscosity
	double viscocity = 0.5;

	// viscosity
	double hyper_viscocity = 50.0;

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
		hx = domain_size/(double)i_res[0];
		hy = domain_size/(double)i_res[1];

		{
			h.data_setall(10);

			double center_x = 0.6;
			double center_y = 0.7;

			if (1)
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

			if (0)
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

		double diff1_x_kernel[3][3] = {{0,0,0},{-1.0/(2.0*hx),0,1.0/(2.0*hx)},{0,0,0}};
		diff1_x.setup_kernel(diff1_x_kernel);
		double diff1_y_kernel[3][3] = {{0,-1.0/(2.0*hy),0},{0,0,0},{0,1.0/(2.0*hy),0}};
		diff1_y.setup_kernel(diff1_y_kernel);
		double diff2_x_kernel[3][3] = {{0,0,0},{1.0/(hx*hx),-2.0/(hx*hx),1.0/(hx*hx)},{0,0,0}};
		diff2_x.setup_kernel(diff2_x_kernel);
		double diff2_y_kernel[3][3] = {{0,1.0/(hy*hy),0},{0,-2.0/(hy*hy),0},{0,1.0/(hy*hy),0}};
		diff2_y.setup_kernel(diff2_y_kernel);

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
			v_t -= (diff2_y(u) + diff2_y(v))*viscocity;
			u_t -= (diff2_x(u) + diff2_x(v))*viscocity;
		}

		if (hyper_viscocity > 0)
		{
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

//        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> " << dt << std::endl;
 //       dt = 0.05;

		h += dt*h_t;
		u += dt*u_t;
		v += dt*v_t;

		timestep_nr++;
	}
};




VisualizationEngine *visualizationEngine;

std::size_t N = 64;

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

		double scale_d = 1.0/(simulationSWE->h-10.0).get_maxAbs();

#pragma omp parallel for simd
		for (std::size_t i = 0; i < simulationSWE->h.array_data_cartesian_length; i++)
		{
			double value;
			// average height
			value = simulationSWE->h.array_data_cartesian_space[i]-10.0;

			// scale
			value *= scale_d;

			// [-1;1] -> [0;255]
			value = (value+1.)*0.5*255.0;

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
#if 0
	runTests(8, 8);

//	runSWE(8, 12);
//	int N = 1024;
//	runSWE(N,N);

#else

	N = 128;

	ProgramCallbacks programCallbacks;

	VisualizationSimulationSWE::getGlobalArgCSingleton() = i_argc;
	VisualizationSimulationSWE::getGlobalArgVSingleton() = i_argv;

	VisualizationEngine(&programCallbacks, "SWE");

	delete visSimulationSWE;

#endif
	return 1;
}
