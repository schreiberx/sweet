/*
 * Copyright 2010 Martin Schreiber
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * CConfig.hpp
 *
 *  Created on: Mar 22, 2010
 * Author: Martin SCHREIBER <schreiberx@gmail.com> (schreiberx@gmail.com)
 */

//#include "simulations/hyperbolic_common/CTsunamiConfig.hpp"

#ifndef VISSWEET_HUD_HPP_
#define VISSWEET_HUD_HPP_

#include <sweet/libgl/hud/GlFreeType.hpp>
#include <sweet/libgl/hud/GlHudConfig.hpp>
#include <sweet/libgl/hud/GlRenderOStream.hpp>
#include <sweet/libgl/hud/GlWindow.hpp>

/*
 * callback handler helper precompiler definitions
 *
 * usage:
 * {	// open scope for each callback handler
 * CCONFIG_CALLBACK_START(CMainClass)
 * }
 */


/*
 * create prefix code for callbacks
 * \param TMainClass	class to which the user pointer is casted to
 *
 * the variable c is created to make the the class of type TMainClass accessible
 */
#define CGUICONFIG_CALLBACK_START(TMainClass)				\
	class CConfigCallback								\
	{													\
		public:											\
		static void callbackHandler(void *p_user_ptr)	\
		{												\
			TMainClass &c = *(TMainClass*)p_user_ptr;

/*
 * create postfix code for callbacks
 */
#define CGUICONFIG_CALLBACK_END()							\
		}												\
	};


/*
 * precompiler directive to handler
 */
#define CGUICONFIG_CALLBACK_HANDLER	(CConfigCallback::callbackHandler)


#define CGUICONFIG_CALLBACK_INSTALL(config_variable)		\
	cGuiConfig.setCallbackHandler(&config_variable, CGUICONFIG_CALLBACK_HANDLER, this);



/**
 * main gui configuration section
 *
 * this class is intended to hold all configurable parameters for the visualization.
 */
class VisSweetHUD
{
	typedef float T;

public:
	GlHudConfig cGlHudConfigMainLeft;
	GlHudConfig cGlHudConfigMainRight;
	GlHudConfig cGlHudConfigVisualization;
	GlHudConfig cGlHudConfigInfo;

private:
	GlWindow windowMain;
	GlWindow windowLights;
	GlWindow windowInfo;


	bool main_drag_n_drop_active;
	bool lights_drag_n_drop_active;
	bool info_drag_n_drop_active;

public:
	/////////////////////////////////////////////////////////////////////

	bool reset;
	bool quit;
	/**
	 * run or pause the simulatino
	 */
	bool run_simulation_timesteps;

	bool render_skybox;

	bool refraction_active;

	bool take_screenshot;
	bool take_screenshot_series;

	bool output_fps_and_parameters;
	bool output_world_position_at_mouse_cursor;

	bool hud_visible;



	/************************************
	 * VISUALIZATION
	 ************************************/
	bool visualization_enabled;
	int dofs_visualization_method;
	int boundary_visualization_method;
	int visualization_render_cluster_borders;
	int visualization_render_wireframe;
	int visualization_simulation_steps_per_frame;


	/************************************
	 * LIGHT
	 ************************************/
	bool light0_rotation;
	GLSL::vec3 light0_world_position;

	bool light0_enabled;

	GLSL::vec3 light0_ambient_color3;
	GLSL::vec3 light0_diffuse_color3;
	GLSL::vec3 light0_specular_color3;

	// water material
	GLSL::vec3 water_ambient_color3;
	GLSL::vec3 water_diffuse_color3;
	GLSL::vec3 water_specular_color3;
	T water_specular_exponent;
	T water_reflectance_at_normal_incidence;
	T water_alpha;

	T refraction_index;

	/////////////////////////////////////////////////////////////////////
private:
	int mouse_x;
	int mouse_y;

	GlFreeType *cGlFreeType = nullptr;

	GlRenderOStream *glRenderOStream = nullptr;


public:
	VisSweetHUD()	:
		main_drag_n_drop_active(false),
		lights_drag_n_drop_active(false),
		info_drag_n_drop_active(false),
		run_simulation_timesteps(true)
	{
		setHudVisibility(false);
	}

	~VisSweetHUD()
	{
		free();
	}

	void free()
	{
		if (cGlFreeType != nullptr)
		{
			delete cGlFreeType;
			cGlFreeType = nullptr;
		}

		if (glRenderOStream != nullptr)
		{
			delete glRenderOStream;
			glRenderOStream = nullptr;
		}
	}


	void setHudVisibility(bool p_render_hud)
	{
		hud_visible = p_render_hud;
	}


	void setup()
	{
		mouse_x = -1;
		mouse_y = -1;

		GlHudConfig::COption o;

		// main window
		windowMain.setPosition(GLSL::ivec2(3, 3));
		windowMain.setSize(GLSL::ivec2(530, 594));
		windowMain.setBackgroundColor(GLSL::vec4(74.0/256.0, 96.0/256.0, 150.0/256.0, 0.7));

		cGlHudConfigMainLeft.insert(o.setupText("|--- MAIN CONTROL ---|"));
		cGlHudConfigMainLeft.insert(o.setupBoolean("Run simulation", &run_simulation_timesteps));
		cGlHudConfigMainLeft.insert(o.setupBoolean("Reset visualization", &reset));		reset = false;
		cGlHudConfigMainLeft.insert(o.setupBoolean("Quit program", &quit));				quit = false;

		cGlHudConfigMainLeft.insert(o.setupLinebreak());
		cGlHudConfigMainLeft.insert(o.setupText("|--- BOGUS STUFF ---|"));
		cGlHudConfigMainLeft.insert(o.setupBoolean("Output FPS and parameters", &output_fps_and_parameters));	output_fps_and_parameters = true;
		cGlHudConfigMainLeft.insert(o.setupBoolean("Output world position at mouse cursor",						&output_world_position_at_mouse_cursor));				output_world_position_at_mouse_cursor = false;
		cGlHudConfigMainLeft.insert(o.setupBoolean("Screenshot (screenshot.bmp)", &take_screenshot));			take_screenshot = false;
		cGlHudConfigMainLeft.insert(o.setupBoolean("Screenshot series", &take_screenshot_series));				take_screenshot_series = false;
		cGlHudConfigMainLeft.insert(o.setupText("(screenshots/screenshot_#####.bmp)"));

		cGlHudConfigMainLeft.insert(o.setupLinebreak());
		cGlHudConfigMainLeft.insert(o.setupText("Press [SPACE] to show/hide HUDs!"));

		/************************************
		 * VISUALIZATION
		 ************************************/
		cGlHudConfigMainLeft.insert(o.setupText("|--- Visualization ---|"));
		cGlHudConfigMainLeft.insert(o.setupBoolean("Enabled", &visualization_enabled));
		visualization_enabled = true;
		cGlHudConfigMainLeft.insert(o.setupInt("Surface visualization Method", &dofs_visualization_method, 1, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()));
		dofs_visualization_method = 0;
		cGlHudConfigMainLeft.insert(o.setupInt("Bathymetry Visualization Method", &boundary_visualization_method, 1, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()));
		boundary_visualization_method = 0;
		cGlHudConfigMainLeft.insert(o.setupInt("Render wireframe", &visualization_render_wireframe, 1, 0, 2));
#if SIMULATION_TSUNAMI_1D
		visualization_render_wireframe = 1;
#else
		visualization_render_wireframe = 0;
#endif
		cGlHudConfigMainLeft.insert(o.setupInt("Cluster Borders", &visualization_render_cluster_borders, 1, 0, 2));
		visualization_render_cluster_borders = 0;
		cGlHudConfigMainLeft.insert(o.setupInt("Simulationsteps per frame", &visualization_simulation_steps_per_frame, 1, 1, std::numeric_limits<int>::max()));
		visualization_simulation_steps_per_frame = 1;

#if 0
		/************************
		 * SIMULATION
		 ************************/
		cGlHudConfigMainRight.insert(o.setupLinebreak());
		cGlHudConfigMainRight.insert(o.setupText("|--- SIMULATION CONTROL ---|"));
		cGlHudConfigMainRight.insert(o.setupBoolean("Run simulation", &cSimulationParameters->run_simulation));
		cGlHudConfigMainRight.insert(o.setupBoolean("Reset simulation fluid", &cSimulationParameters->simulation_reset));
		cGlHudConfigMainRight.insert(o.setupBoolean("Insert Random Raindrops", &cSimulationParameters->simulation_random_raindrops_activated));

		/**
		 * Tsunami simulation parameters
		 */
		cGlHudConfigMainRight.insert(o.setupText("|--- Setup / Dataset parameters ---|"));
		cGlHudConfigMainRight.insert(o.setupInt("World scene", &cSimulationParameters->simulation_world_scene_id, 1, -999, 999));
		cGlHudConfigMainRight.insert(o.setupInt("Terrain scene", &cSimulationParameters->simulation_terrain_scene_id, 1, 0, 7));
		cGlHudConfigMainRight.insert(o.setupInt("Water height scene", &cSimulationParameters->simulation_water_surface_scene_id, 1, 0, 7));

		cGlHudConfigMainRight.insert(o.setupText("|--- EdgeComm parameters ---|"));
		cGlHudConfigMainRight.insert(o.setupFloat("Gravitation length (TODO): ", &cSimulationParameters->simulation_parameter_gravitation, 0.1, -100, 0));
		cGlHudConfigMainRight.insert(o.setupFloatHI("Timestep size: ", &cSimulationParameters->simulation_parameter_timestep_size, 0.01));
		//cGlHudConfigMainRight.insert(o.setupFloatHI("CFL: ", &cSimulationParameters->simulation_parameter_cfl, 0.001, 0, 1));



		cGlHudConfigMainRight.insert(o.setupLinebreak());
		cGlHudConfigMainRight.insert(o.setupText("|--- DATASET PARAMETERS ---|"));
		cGlHudConfigMainRight.insert(o.setupText("Generic Parameters:"));
		cGlHudConfigMainRight.insert(o.setupFloatHI(" + Bathymetry default distance: ", &cSimulationParameters->simulation_dataset_terrain_default_distance, 0.1, -std::numeric_max<float>(), std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloat(" + Water Surface default height: ", &cSimulationParameters->simulation_dataset_water_surface_default_displacement, 0.01, -std::numeric_max<float>(), std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupText("Surface setup for column (1):"));
		cGlHudConfigMainRight.insert(o.setupFloat(" + X-Position: ", &cSimulationParameters->simulation_dataset_cylinder_posx, 0.02, -std::numeric_max<float>(), std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloat(" + Y-Position: ", &cSimulationParameters->simulation_dataset_cylinder_posy, 0.02, -std::numeric_max<float>(), std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloatHI(" + Radius: ", &cSimulationParameters->simulation_dataset_cylinder_radius, 0.02, 0.000001, std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloat(" + Height: ", &cSimulationParameters->simulation_dataset_cylinder_height, 0.02, -std::numeric_max<float>(), std::numeric_max<float>()));

		/************************
		 * Coarsen/Refine
		 ************************/
		cGlHudConfigMainRight.insert(o.setupLinebreak());
		cGlHudConfigMainRight.insert(o.setupText("|--- ADAPTIVE PARAMETERS ---|"));
		cGlHudConfigMainRight.insert(o.setupFloatHI("Coarsen height threshold", &cSimulationParameters->coarsen_parameter_0, 0.02, 0.0000001, std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloatHI("Refine height threshold", &cSimulationParameters->refine_parameter_0, 0.02, 0.0000001, std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloatHI("Coarsen slope threshold (TODO)", &cSimulationParameters->coarsen_parameter_1, 0.02, 0.0000001, std::numeric_max<float>()));
		cGlHudConfigMainRight.insert(o.setupFloatHI("Refine slope threshold (TODO)", &cSimulationParameters->refine_parameter_1, 0.02, 0.0000001, std::numeric_max<float>()));



		/************************
		 * PARALLELIZATION
		 ************************/
		cGlHudConfigMainRight.insert(o.setupLinebreak());
		cGlHudConfigMainRight.insert(o.setupText("|--- PARALLELIZATION ---|"));
		cGlHudConfigMainRight.insert(o.setupUInt("Split Elements", &cSimulationParameters->splitting_cluster_split_workload_size, 1, 2, std::numeric_limits<int>::max()));
		cGlHudConfigMainRight.insert(o.setupUInt("Join Elements", &cSimulationParameters->splitting_cluster_join_workload_size, 1, 0, std::numeric_limits<int>::max()));
		cGlHudConfigMainRight.insert(o.setupInt("Number of Threads to use", &simulation_number_of_threads, 1, 1, std::numeric_limits<int>::max()));	// number is set by gui
#endif


#if 0
		/***********************
		 * LIGHT
		 ***********************/
		cGlHudConfigVisualization.insert(o.setupText("|--- LIGHT 0 ---|"));

		cGlHudConfigVisualization.insert(o.setupBoolean("enabled", &light0_enabled));					light0_enabled = true;
		cGlHudConfigVisualization.insert(o.setupFloat("position X", &light0_world_position[0], 0.1));	light0_world_position = GLSL::vec3(0.5, 4.5, 0.5);
		cGlHudConfigVisualization.insert(o.setupFloat("position Y", &light0_world_position[1], 0.1));
		cGlHudConfigVisualization.insert(o.setupFloat("position Z", &light0_world_position[2], 0.1));

		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT R", &light0_ambient_color3[0], 0.01, 0.0, 1.0));	light0_ambient_color3 = GLSL::vec3(1,1,1);
		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT G", &light0_ambient_color3[1], 0.01, 0.0, 1.0));
		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT B", &light0_ambient_color3[2], 0.01, 0.0, 1.0));
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE R", &light0_diffuse_color3[0], 0.01, 0.0, 1.0));	light0_diffuse_color3 = GLSL::vec3(1,1,1);
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE G", &light0_diffuse_color3[1], 0.01, 0.0, 1.0));
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE B", &light0_diffuse_color3[2], 0.01, 0.0, 1.0));
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR R", &light0_specular_color3[0], 0.01, 0.0, 1.0));	light0_specular_color3 = GLSL::vec3(1,1,1);
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR G", &light0_specular_color3[1], 0.01, 0.0, 1.0));
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR B", &light0_specular_color3[2], 0.01, 0.0, 1.0));

		cGlHudConfigVisualization.insert(o.setupBoolean("Light0 rotation", &light0_rotation));				light0_rotation = false;


		// water material
		cGlHudConfigVisualization.insert(o.setupLinebreak());
		cGlHudConfigVisualization.insert(o.setupText("|--- WATER parameters ---|"));
		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT R", &water_ambient_color3[0], 0.01, 0.0, 111.0));		water_ambient_color3 = GLSL::vec3(0.78, 0.88, 1.0);
		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT G", &water_ambient_color3[1], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("AMBIENT B", &water_ambient_color3[2], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE R", &water_diffuse_color3[0], 0.01, 0.0, 111.0));		water_diffuse_color3 = GLSL::vec3(0.0f);
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE G", &water_diffuse_color3[1], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("DIFFUSE B", &water_diffuse_color3[2], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR exponent", &water_specular_exponent, 0.1, 0.0));	water_specular_exponent = 20.0f;
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR R", &water_specular_color3[0], 0.01, 0.0, 111.0));	water_specular_color3 = GLSL::vec3(6.5,6.3,6);
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR G", &water_specular_color3[1], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("SPECULAR B", &water_specular_color3[2], 0.01, 0.0, 111.0));
		cGlHudConfigVisualization.insert(o.setupFloat("Water normal reflection", &water_reflectance_at_normal_incidence, 0.001, 0.0, 1.0));	water_reflectance_at_normal_incidence = 0.07;
		cGlHudConfigVisualization.insert(o.setupFloat("alpha", &water_alpha, 0.01, 0.0, 111.0));	water_alpha = 0.9;

		windowLights.setPosition(GLSL::ivec2(536, 3));
		windowLights.setBackgroundColor(GLSL::vec4(128.0/256.0, 128.0/256.0, 128.0/256.0, 0.7));
		windowLights.setSize(GLSL::ivec2(270, 594));
#endif



		/***********************
		 * CONFIG INFO
		 ***********************/
#if 1
		cGlHudConfigInfo.insert(o.setupText("|--- INFORMATION ---|"));
		cGlHudConfigInfo.insert(o.setupText(""));
		cGlHudConfigInfo.insert(o.setupText(" [backspace]: disable this HUD"));
		cGlHudConfigInfo.insert(o.setupText(" [space]: Pause/continue simulation"));
		cGlHudConfigInfo.insert(o.setupText(" j: single simulation time step"));
		cGlHudConfigInfo.insert(o.setupText(" r: reset simulation"));
		cGlHudConfigInfo.insert(o.setupText(" 1-9: N^2 simulation steps per frame"));
		cGlHudConfigInfo.insert(o.setupText(" q: quit"));
		cGlHudConfigInfo.insert(o.setupText(" e: reset simulation & view"));
		cGlHudConfigInfo.insert(o.setupText(" v/V: increase/decrease visualization id"));
		cGlHudConfigInfo.insert(o.setupText(" a, s, d, w: movement"));

		windowInfo.setPosition(GLSL::ivec2(10, 10));
		windowInfo.setBackgroundColor(GLSL::vec4(256.0/256.0, 256.0/256.0, 128.0/256.0, 0.7));
		windowInfo.setSize(GLSL::ivec2(360, 320));
#endif
	}


	/**
	 * finally assemble GUI
	 */
	void assembleGui(
			GlFreeType &i_cGlFreeType,
			GlRenderOStream &i_cGlRenderOStream
	)
	{
#if 0

		cGlFreeType = &i_cGlFreeType;

#else

		free();

		cGlFreeType = new GlFreeType();
		cGlFreeType->loadFont(18, false);
		if (!cGlFreeType->valid)
		{
			std::cerr << "Loading font failed" << std::endl;
			exit(-1);
		}

		glRenderOStream = new GlRenderOStream(*cGlFreeType);
#endif

		int sx, sy;
#if 0
		sx = 10;
		sy = windowMain.size[1]-10;
		cGlHudConfigMainLeft.setup(*cGlFreeType, i_cGlRenderOStream, sx, sy);
		cGlHudConfigMainRight.setup(*cGlFreeType, i_cGlRenderOStream, sx + cGlHudConfigMainLeft.area_width-20, sy);

		sx = 10;
		sy = windowLights.size[1]-10;
		cGlHudConfigVisualization.setup(*cGlFreeType, i_cGlRenderOStream, sx, sy);
#endif

		sx = 10;
		sy = windowInfo.size[1]-10;
		cGlHudConfigInfo.setup(*cGlFreeType, *glRenderOStream, sx, sy);
	}


	/**
	 * setup a callback handler which is called, when the value changed
	 */
	void setCallbackHandler(
			void *value_ptr,							///< pointer to value
			void (*callback_handler)(void *user_ptr),	///< callback handler
			void *user_ptr
		)
	{
//		cGlHudConfigMainLeft.set_callback(value_ptr, callback_handler, user_ptr);
//		cGlHudConfigMainRight.set_callback(value_ptr, callback_handler, user_ptr);

//		cGlHudConfigVisualization.set_callback(value_ptr, callback_handler, user_ptr);

		cGlHudConfigInfo.set_callback(value_ptr, callback_handler, user_ptr);
	}


	/**
	 * render windowMain with configuration information
	 */
	void render()
	{
		if (hud_visible)
		{
#if 0
			cGlFreeType->viewportChanged(windowLights.size.data);
			windowLights.startRendering();
			cGlHudConfigVisualization.renderConfigContent();
			windowLights.finishRendering();

			cGlFreeType->viewportChanged(windowMain.size.data);
			windowMain.startRendering();
			cGlHudConfigMainLeft.renderConfigContent();
			cGlHudConfigMainRight.renderConfigContent();
			windowMain.finishRendering();
#endif

			cGlFreeType->viewportChanged(windowInfo.size.data);
			windowInfo.startRendering();
			cGlHudConfigInfo.renderConfigContent();
			windowInfo.finishRendering();
		}
	}

	bool mouse_wheel(int x, int y)
	{
		if (hud_visible)
		{
			if (windowInfo.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigInfo.mouse_wheel(x, y);
				return true;
			}
#if 0
			else if (windowMain.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigMainLeft.mouse_wheel(x, y);
				cGlHudConfigMainRight.mouse_wheel(x, y);
				return true;
			}
			else if (windowLights.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigVisualization.mouse_wheel(x, y);
				return true;
			}
#endif
		}

		return false;
	}

	bool mouse_button_down(int button)
	{
		if (hud_visible)
		{
			if (windowInfo.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigInfo.mouse_button_down(button);

				if (cGlHudConfigInfo.active_option == nullptr)
				{
					if (button == RenderWindow::MOUSE_BUTTON_LEFT)
						info_drag_n_drop_active = true;
				}
				return true;
			}
			else if (windowMain.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigMainLeft.mouse_button_down(button);
				cGlHudConfigMainRight.mouse_button_down(button);

				if (cGlHudConfigMainLeft.active_option == nullptr && cGlHudConfigMainRight.active_option == nullptr)
				{
					if (button == RenderWindow::MOUSE_BUTTON_LEFT)
						main_drag_n_drop_active = true;
				}
				return true;
			}
			else if (windowLights.inWindow(mouse_x, mouse_y))
			{
				cGlHudConfigVisualization.mouse_button_down(button);

				if (cGlHudConfigVisualization.active_option == nullptr && cGlHudConfigVisualization.active_option == nullptr)
				{
					if (button == RenderWindow::MOUSE_BUTTON_LEFT)
						lights_drag_n_drop_active = true;
				}
				return true;
			}
		}

		return false;
	}

	void mouse_button_up(int button)
	{
		if (hud_visible)
		{
			cGlHudConfigMainLeft.mouse_button_up(button);
			cGlHudConfigMainRight.mouse_button_up(button);

			cGlHudConfigVisualization.mouse_button_up(button);
			cGlHudConfigInfo.mouse_button_up(button);
		}

		if (button == RenderWindow::MOUSE_BUTTON_LEFT)
		{
			// always disable drag_n_drop, when left mouse button is released
			main_drag_n_drop_active = false;
			lights_drag_n_drop_active = false;
			info_drag_n_drop_active = false;
		}
	}

	void mouse_motion(int x, int y)
	{
		int dx = x - mouse_x;
		int dy = y - mouse_y;

		if (info_drag_n_drop_active)
			windowInfo.setPosition(windowInfo.pos + GLSL::ivec2(dx, dy));

		if (main_drag_n_drop_active)
			windowMain.setPosition(windowMain.pos + GLSL::ivec2(dx, dy));

		if (lights_drag_n_drop_active)
			windowLights.setPosition(windowLights.pos + GLSL::ivec2(dx, dy));

		mouse_x = x;
		mouse_y = y;

		if (hud_visible)
		{
			cGlHudConfigInfo.mouse_motion(mouse_x-windowInfo.pos[0], mouse_y-windowInfo.pos[1]);

			cGlHudConfigMainLeft.mouse_motion(mouse_x-windowMain.pos[0], mouse_y-windowMain.pos[1]);
			cGlHudConfigMainRight.mouse_motion(mouse_x-windowMain.pos[0], mouse_y-windowMain.pos[1]);

			cGlHudConfigVisualization.mouse_motion(mouse_x-windowLights.pos[0], mouse_y-windowLights.pos[1]);
		}
	}
};


#endif
