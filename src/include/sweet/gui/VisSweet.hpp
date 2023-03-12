/*
 * VisSweet.hpp
 *
 *  Created on: 30 Jun 2015
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_VISSWEET_HPP_
#define SRC_EXAMPLES_VISSWEET_HPP_

#include <limits>
#include <string>
#include <cmath>
#include <sweet/gui/VisSweetHUD.hpp>
#include "../libgl/draw/GlDrawCube.hpp"
#include "../libgl/shaders/shader_blinn/CShaderBlinn.hpp"
#include "../libgl/VisualizationEngine.hpp"
#include "../libgl/core/GlTexture.hpp"
#include "../libgl/draw/GlDrawQuad.hpp"

#ifndef SWEET_USE_SPHERE_SPECTRAL_SPACE
#define SWEET_USE_SPHERE_SPECTRAL_SPACE 1
#endif

#if SWEET_USE_SPHERE_SPECTRAL_SPACE
	#include "../libgl/draw/GlDrawSphereSph.hpp"
#endif
#include "../libgl/hud/GlFreeType.hpp"
#include "../libgl/hud/GlRenderOStream.hpp"

#include <sweet/core/plane/PlaneData_Physical.hpp>


class SimulationGUICallbacks
{
public:
	/**
	 * postprocessing of frame: do time stepping
	 */
	virtual
	void vis_post_frame_processing(int i_num_iterations) = 0;

	virtual
	void vis_getDataArray(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *vis_reset
	) = 0;

	virtual
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline) = 0;


	virtual
	void vis_pause() = 0;


	virtual
	void vis_keypress(int i_key) = 0;

	virtual
	bool should_quit() = 0;

	virtual
	bool runTimestep() = 0;

	virtual
	bool reset() = 0;

	~SimulationGUICallbacks() {}
};



/**
 * A visualization class specifically designed for SWEET applications
 */
class VisSweet	:
		public VisualizationEngine::ProgramCallbacks
{
	/**
	 * Simulation class from the SWEET-using application
	 *
	 * Certain interfaces have to be implemented, see other programs for these interfaces
	 */
	SimulationGUICallbacks *simCallbacks;

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



	bool vis_shouldQuit()
	{
		return simCallbacks->should_quit();
	}



	void vis_render()
	{
		const sweet::PlaneData_Physical *ro_visPlaneData;
		double aspect_ratio = 0;
		int render_primitive = 0;
		void *bogus_data;

		vis_min = std::numeric_limits<double>::infinity();
		vis_max = std::numeric_limits<double>::infinity();
		bool reset = false;


		simCallbacks->vis_getDataArray(
				&ro_visPlaneData,
				&aspect_ratio,
				&render_primitive,
				&bogus_data,
				&vis_min,
				&vis_max,
				&reset
		);

		sweet::PlaneData_Physical &visData = (sweet::PlaneData_Physical&)*ro_visPlaneData;

		if (std::isinf(vis_min))
		{
			vis_min = visData.physical_reduce_min();
			vis_max = visData.physical_reduce_max();

			vis_max = std::max(vis_max, vis_min+1e-20);	//< avoid numerical issues if min == max
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


		double real_delta = vis_max-vis_min;

		double inv_delta = 1.0/real_delta;

		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (	std::size_t i = 0;
				i < visData.planeDataConfig->physical_array_data_number_of_elements;
				i++
		)
		{
			double value = (visData.physical_space_data[i]-vis_min)*inv_delta;
			value *= 255.0;

			texture_data[i] = (unsigned char)std::min(255.0, std::max(0.0, value));//			texture_data[i] = 128;
		}

		glTexture->bind();
		glTexture->setData(texture_data);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

			visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.use();

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

				glDrawQuad->renderWithoutProgram();
			}
			else
			{
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
				VisualizationEngine::EngineState::Matrices &m = visualizationEngine->engineState->matrices;
				visualizationEngine->engineState->commonShaderPrograms.shaderTexturizeRainbowmap.pvm_matrix_uniform.set(
						m.pvm
				);

				/*
				 * Initialize on the first time we do this
				 */
				if (glDrawSphereSph == nullptr)
				{
					sweet::SphereData_Config *sphereDataConfig = (sweet::SphereData_Config*)bogus_data;

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

			bool o_replace_comma_with_linebreaks = true;
			std::string status_string = simCallbacks->vis_getStatusString(o_replace_comma_with_linebreaks);
			if (o_replace_comma_with_linebreaks)
				std::replace(status_string.begin(), status_string.end(), ',', '\n');

			glFreeType->setPosition(10, visualizationEngine->renderWindow->window_height-font_size-10);
			glFreeType->renderString(status_string.c_str());

			visSweetHUD->render();
		}


		// execute simulation time step
		simCallbacks->vis_post_frame_processing(sim_runs_per_frame);
	}


	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		std::ostringstream ss;
		ss <<  simCallbacks->vis_getStatusString(o_replace_commas_with_newline) << ", vis min/max: " << vis_min << ", " << vis_max;
		return ss.str();
	}

	void update_font()
	{
		// Approx. desired font size on standard display and standard viewport height
		font_size = 14.0;

		// rescale with viewport size
		font_size *= std::max(viewport_width, viewport_height)/(2.0*800.0);

		// Limit font size
		font_size = std::min(50.0, font_size);
		font_size = std::max(14.0, font_size);

		// Incorporate HiDPI information
		float ddpi;
		SDL_GetDisplayDPI(0, &ddpi, nullptr, nullptr);
		font_size *= ddpi/100.0;

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
			simCallbacks->reset();
			break;

		case ' ':
			simCallbacks->vis_pause();
			break;

		case SDLK_BACKSPACE:
			hud_visible = !hud_visible;
			visSweetHUD->setHudVisibility(hud_visible);
			break;

		case 'j':
			simCallbacks->runTimestep();
			break;

		default:
			simCallbacks->vis_keypress(i_key);
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
	VisSweet(SimulationGUICallbacks &i_simCallbacks)	:
		glDrawQuad(nullptr),
#if SWEET_USE_SPHERE_SPECTRAL_SPACE
		glDrawSphereSph(nullptr),
#endif

		visSweetHUD(nullptr),
		glFreeType(0),
		glRenderOStream(nullptr)
	{
		simCallbacks = &i_simCallbacks;

		VisualizationEngine(this, "SWEET");

	}

};



#endif
