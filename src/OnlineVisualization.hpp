/*
 * Visualization.hpp
 *
 *  Created on: 29 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef ONLINE_VISUALIZATION_HPP_
#define ONLINE_VISUALIZATION_HPP_


#include "libgl/VisualizationEngine.hpp"
#include "libgl/draw/GlDrawQuad.hpp"
#include "libgl/draw/GlDrawCube.hpp"
#include "libgl/core/CGlTexture.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinn.hpp"


class OnlineVisualization	:
		public VisualizationEngine::ProgramCallbacks
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
	VisualizationEngine **getVisualizationEngineSingletonPtr()
	{
		VisualizationEngine *visualizationEngine = nullptr;
		return &visualizationEngine;
	}

public:
	OnlineVisualization()
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

	~OnlineVisualization()
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


	OnlineVisualization *onlineVisualization;

	void setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;
		onlineVisualization = new OnlineVisualization;
	}

	void render()
	{
		onlineVisualization->frame_callback();
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





#endif /* VISUALIZATION_HPP_ */
