/*
 * VisSweet.hpp
 *
 *  Created on: 30 Jun 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_EXAMPLES_VISSWEET_HPP_
#define SRC_EXAMPLES_VISSWEET_HPP_

#include <limits>
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

	double viz_min = 0;
	double viz_max = 0;

	int viewport_width = -1, viewport_height = -1;

	double font_size = -1;

	void vis_setup(VisualizationEngine *i_visualizationEngine)
	{
		visualizationEngine = i_visualizationEngine;

		glDrawQuad = new GlDrawQuad;
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		glDrawSphereSph = nullptr;
#endif

		glFreeType = new GlFreeType;

		update_font();

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

		viz_min = std::numeric_limits<double>::infinity();
		viz_max = std::numeric_limits<double>::infinity();
		bool reset = false;

		simulation->vis_get_vis_data_array(
				&ro_visPlaneData,
				&aspect_ratio,
				&render_primitive,
				&bogus_data,
				&viz_min,
				&viz_max,
				&reset
		);

		PlaneData &visData = (PlaneData&)*ro_visPlaneData;

		visData.request_data_physical();

		if (std::isinf(viz_min))
		{
			viz_min = visData.reduce_min();
			viz_max = visData.reduce_max();

			viz_max = std::max(viz_max, viz_min+1e-20);	//< avoid numerical issues if min == max
		}


		if (glTexture == nullptr || reset)
		{
			delete glTexture;
			glTexture = new GlTexture(GL_TEXTURE_2D, GL_RED, GL_RED, GL_UNSIGNED_BYTE);
			glTexture->bind();
			glTexture->resize(visData.planeDataConfig->physical_data_size[0], visData.planeDataConfig->physical_data_size[1]);
			glTexture->unbind();

			texture_data = new unsigned char[visData.planeDataConfig->physical_array_data_number_of_elements];
		}


		double real_delta = viz_max-viz_min;

		double inv_delta = 1.0/real_delta;

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (	std::size_t i = 0;
				i < visData.planeDataConfig->physical_array_data_number_of_elements;
				i++
		)
		{
			double value = (visData.physical_space_data[i]-viz_min)*inv_delta;
			value *= 255.0;

			texture_data[i] = (unsigned char)std::min(255.0, std::max(0.0, value));//			texture_data[i] = 128;
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
					SphereData_Config *sphereDataConfig = (SphereData_Config*)bogus_data;

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

			glFreeType->setPosition(10, visualizationEngine->renderWindow->window_height-font_size-10);
			glFreeType->renderString(status_string.c_str());

			visSweetHUD->render();
		}


		// execute simulation time step
		simulation->vis_post_frame_processing(sim_runs_per_frame);
	}



	const char* vis_getStatusString()
	{
		static char title_string[4096];
		sprintf(title_string, "%s, viz min/max: %.12e/%.12e", simulation->vis_get_status_string(), viz_min, viz_max);
		return title_string;
	}

	void update_font()
	{
		// reload font
		font_size = std::max(viewport_width, viewport_height);
		font_size /= 40.0;

		font_size = std::min(50.0, font_size);
		font_size = std::max(10.0, font_size);

		glFreeType->loadFont(font_size, false);
		if (!glFreeType->valid)
		{
			SWEETError("Loading font failed");
			exit(-1);
		}
	}


	void vis_viewportChanged(int i_width, int i_height)
	{
		viewport_width = i_width;
		viewport_height = i_height;

		update_font();
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
