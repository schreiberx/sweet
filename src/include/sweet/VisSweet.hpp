/*
 * VisSweet.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_VISSWEET_HPP_
#define SRC_EXAMPLES_VISSWEET_HPP_

#include <sweet/VisSweetHUD.hpp>
#include "../libgl/draw/GlDrawCube.hpp"
#include "../libgl/shaders/shader_blinn/CShaderBlinn.hpp"
#include "../libgl/VisualizationEngine.hpp"
#include "../libgl/core/GlTexture.hpp"
#include "../libgl/draw/GlDrawQuad.hpp"

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
	#include "../libgl/draw/GlDrawSphereSph.hpp"
#endif
#include "../libgl/hud/GlFreeType.hpp"
#include "../libgl/hud/GlRenderOStream.hpp"
#include <sweet/plane/PlaneData.hpp>




/**
 * A visualization class specifically designed for SWEET applications
 */
template <typename SimT>
class VisSweet	:
		public VisualizationEngine::ProgramCallbacks
{
	/**
	 * Simulation class from the SWEET-using application
	 *
	 * Certain interfaces have to be implemented, see other programs for these interfaces
	 */
	SimT *simulation;

	GlDrawQuad *glDrawQuad;

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
	GlDrawSphereSph *glDrawSphereSph;
#endif

	GlTexture *glTexture = nullptr;
	unsigned char *texture_data = nullptr;

	bool hud_visible = true;
	VisSweetHUD *visSweetHUD;
	GlFreeType *glFreeType;
	GlRenderOStream *glRenderOStream;

	int sim_runs_per_frame = 1;

	VisualizationEngine *visualizationEngine;

	double vis_min = 0;
	double vis_max = 0;


	void vis_setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;

		glDrawQuad = new GlDrawQuad;
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		glDrawSphereSph = nullptr;
#endif

		glFreeType = new GlFreeType;
		glFreeType->loadFont(14, false);
		if (!glFreeType->valid)
		{
//			std::cerr << cGlFreeType.error << std::endl;
			exit(-1);
		}
		//CError_AppendReturn(cGlFreeType);

		glRenderOStream = new GlRenderOStream(*glFreeType);

		visSweetHUD = new VisSweetHUD;
		visSweetHUD->setup();
		visSweetHUD->assembleGui(*glFreeType, *glRenderOStream);
		visSweetHUD->setHudVisibility(hud_visible);

	}



	bool should_quit()
	{
		return simulation->should_quit();
	}



	void vis_render()
	{
		const PlaneData *ro_visPlaneData;
		double aspect_ratio = 0;
		int render_primitive = 0;
		void *bogus_data;

		simulation->vis_get_vis_data_array(
				&ro_visPlaneData,
				&aspect_ratio,
				&render_primitive,
				&bogus_data
		);

		PlaneData &visData = (PlaneData&)*ro_visPlaneData;


		if (glTexture == nullptr)
		{
			glTexture = new GlTexture(GL_TEXTURE_2D, GL_RED, GL_RED, GL_UNSIGNED_BYTE);
			glTexture->bind();
			glTexture->resize(visData.planeDataConfig->physical_data_size[0], visData.planeDataConfig->physical_data_size[1]);
			glTexture->unbind();

			texture_data = new unsigned char[visData.planeDataConfig->physical_array_data_number_of_elements];
		}

		visData.request_data_physical();

		vis_min = visData.reduce_min();
		vis_max = visData.reduce_max();

		vis_max = std::max(vis_max, vis_min+1e-20);	//< avoid numerical issues if min == max

		double real_delta = vis_max-vis_min;
//		vis_min -= real_delta*0.03;
//		vis_max += real_delta*0.03;

		double inv_delta = 1.0/real_delta;

//#pragma omp parallel for OPENMP_PAR_SIMD
		for (std::size_t i = 0; i < visData.planeDataConfig->physical_array_data_number_of_elements; i++)
		{
			double value = (visData.physical_space_data[i]-vis_min)*inv_delta;
			value *= 255.0;

			texture_data[i] = (unsigned char)std::min(255.0, std::max(0.0, value));
		}

		glTexture->bind();
		glTexture->setData(texture_data);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

			visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.use();
//			visualizationEngine->engineState->commonShaderPrograms.shaderBlinn.use();

			if (render_primitive == 0)
			{
				double scale_x = 1.0;
				double scale_y = aspect_ratio;

				if (aspect_ratio > 1.0)
				{
					scale_x /= aspect_ratio;
					scale_y = 1.0;
				}

				visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.pvm_matrix_uniform.set(
						visualizationEngine->engineState->matrices.pvm*GLSL::scale((float)scale_x, (float)scale_y, (float)1.0)
				);

//				visualizationEngine->engineState->commonShaderPrograms.shaderBlinn.pvm_matrix_uniform.set(
//						visualizationEngine->engineState->matrices.pvm*GLSL::scale((float)scale_x, (float)scale_y, (float)1.0)
//				);

				glDrawQuad->renderWithoutProgram();
			}
			else
			{
#if SWEET_USE_SPHERE_SPECTRAL_SPACE


				VisualizationEngine::EngineState::Matrices &m = visualizationEngine->engineState->matrices;
				visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.pvm_matrix_uniform.set(
						m.pvm
				);

				if (glDrawSphereSph == nullptr)
				{
					SphereDataConfig *sphereDataConfig = (SphereDataConfig*)bogus_data;

					glDrawSphereSph = new GlDrawSphereSph;
					glDrawSphereSph->initSphere(sphereDataConfig);
				}

				glDrawSphereSph->renderWithoutProgram();

#endif
			}

			visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.disable();

		glTexture->unbind();

		if (hud_visible)
		{
			glFreeType->viewportChanged(visualizationEngine->renderWindow->window_width, visualizationEngine->renderWindow->window_height);

			std::string status_string = simulation->vis_get_status_string();
			std::replace(status_string.begin(), status_string.end(), ',', '\n');

			glFreeType->setPosition(10, visualizationEngine->renderWindow->window_height-16);
			glFreeType->renderString(status_string.c_str());
			visSweetHUD->render();
		}


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
		if (i_key >= '1' && i_key <= '9')
		{
			sim_runs_per_frame = std::pow(2, i_key-'1');
			return;
		}

		switch(i_key)
		{
		case 'r':
			simulation->reset();
			break;

		case ' ':
			simulation->vis_pause();
			break;

		case SDLK_BACKSPACE:
			hud_visible = !hud_visible;
			visSweetHUD->setHudVisibility(hud_visible);
			break;

		case 'j':
			simulation->run_timestep();
			break;

		default:
			simulation->vis_keypress(i_key);
			break;
		}
	}


	void vis_shutdown()
	{
		delete visSweetHUD;
		delete glRenderOStream;
		delete glFreeType;

		delete [] texture_data;
		delete glTexture;

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		delete glDrawSphereSph;
#endif
		delete glDrawQuad;
	}



public:
	VisSweet(SimT *i_simulation)	:
		glDrawQuad(nullptr),
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		glDrawSphereSph(nullptr),
#endif

		visSweetHUD(nullptr),
		glFreeType(0),
		glRenderOStream(nullptr)
	{
		simulation = i_simulation;

		VisualizationEngine(this, "SWEET");

	}

};



#endif /* SRC_EXAMPLES_VISSWEET_HPP_ */
