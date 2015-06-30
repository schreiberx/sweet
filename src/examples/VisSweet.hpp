/*
 * VisSweet.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_VISSWEET_HPP_
#define SRC_EXAMPLES_VISSWEET_HPP_

#include <sweet/DataArray.hpp>

#include "libgl/VisualizationEngine.hpp"
#include "libgl/draw/GlDrawQuad.hpp"
#include "libgl/draw/GlDrawCube.hpp"
#include "libgl/core/CGlTexture.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinn.hpp"



/**
 * A visualization class specifically designed for SWEET applications
 */

template <typename SimT>
class VisSweet	:
		public VisualizationEngine::ProgramCallbacks
{
	SimT *simulation;

	CGlDrawQuad *cGlDrawQuad;
	CGlTexture *cGlTexture = nullptr;
	unsigned char *texture_data = nullptr;

	int sim_runs_per_frame = 1;

	VisualizationEngine *visualizationEngine;

	double vis_min = 0;
	double vis_max = 0;

	void vis_setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;

		cGlDrawQuad = new CGlDrawQuad;

	}

	void vis_render()
	{
		DataArray<2> &visData = simulation->vis_get_vis_data_array();

		if (cGlTexture == nullptr)
		{
			cGlTexture = new CGlTexture(GL_TEXTURE_2D, GL_RGBA, GL_BLUE, GL_UNSIGNED_BYTE);
			cGlTexture->bind();
			cGlTexture->resize(visData.resolution[1], visData.resolution[0]);
			cGlTexture->unbind();

			texture_data = new unsigned char[visData.array_data_cartesian_length];
		}

		visData.requestDataInCartesianSpace();

		vis_min = visData.reduce_min();
		vis_max = visData.reduce_max();

		vis_max = std::max(vis_max, vis_min+0.000001);	/// avoid numerical issues if min == max

		double real_delta = vis_max-vis_min;
		vis_min -= real_delta*0.1;
		vis_max += real_delta*0.1;

		double inv_delta = 1.0/(vis_max-vis_min);

#pragma omp parallel for simd
		for (std::size_t i = 0; i < visData.array_data_cartesian_length; i++)
		{
			double value = (visData.array_data_cartesian_space[i]-vis_min)*inv_delta;
			value *= 255.0;

			texture_data[i] = value;
		}


		visualizationEngine->engineState->commonShaderPrograms.shaderTexturize.use();
		visualizationEngine->engineState->commonShaderPrograms.shaderTexturize.pvm_matrix_uniform.set(visualizationEngine->engineState->matrices.pvm);

			cGlTexture->bind();
			cGlTexture->setData(texture_data);

				cGlDrawQuad->render();

			cGlTexture->unbind();
		visualizationEngine->engineState->commonShaderPrograms.shaderBlinn.disable();

		// execute simulation time step
		simulation->vis_post_frame_processing(sim_runs_per_frame);

	}

	const char* vis_getStatusString()
	{
		static char title_string[1024];
		sprintf(title_string, "%s, vis min/max: %f/%f", simulation->vis_get_status_string(), vis_min, vis_max);
		return title_string;
	}

	void vis_viewportChanged(int i_width, int i_height)
	{

	}

	void vis_keypress(char i_key)
	{
		switch(i_key)
		{
		case '1':
			sim_runs_per_frame = 1;
			break;

		case '2':
			sim_runs_per_frame = 10;
			break;

		case 'r':
			simulation->reset();
			break;

		case ' ':
			simulation->vis_pause();
			break;

		case 'j':
			simulation->run_timestep();
			break;
		}
	}

	void vis_shutdown()
	{
		delete [] texture_data;
		delete cGlTexture;
		delete cGlDrawQuad;

	}

public:
	VisSweet(SimT *i_simulation)
	{
		simulation = i_simulation;
		VisualizationEngine(this, "SWE");
	}

};



#endif /* SRC_EXAMPLES_VISSWEET_HPP_ */
