/*
 * Engine.hpp
 *
 *  Created on: Dec 23, 2014
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_LIBGL_VISUALIZATION_ENGINE_HPP_
#define SRC_INCLUDE_LIBGL_VISUALIZATION_ENGINE_HPP_


#include <libgl/tools/RenderWindow.hpp>

#include <libgl/engine/Lights.hpp>
#include <libgl/engine/InputStateMouse.hpp>
#include <libgl/engine/camera/Camera1stPerson.hpp>
#include <libgl/engine/GlCommonShaderPrograms.hpp>
#include <libgl/engine/EyeBall.hpp>
#include <libgl/engine/Time.hpp>
#include <libgl/hud/GlFreeType.hpp>



class VisualizationEngine
{

public:
	/*
	 * the members of this class have to be overwritten and are called back
	 */
	class ProgramCallbacks
	{
	public:
		virtual void vis_setup(VisualizationEngine *i_visualizationEngine) = 0;

		virtual void vis_render() = 0;

		virtual bool should_quit() = 0;

		virtual const char* vis_getStatusString() = 0;

		virtual void vis_viewportChanged(int i_width, int i_height) = 0;

		virtual void vis_keypress(char i_key) = 0;

		virtual void vis_shutdown() = 0;

		virtual ~ProgramCallbacks()
		{
		}
	};

	ProgramCallbacks *programCallbacks;



public:
	class EngineState
	{
	public:
		/*
		 * quit application?
		 */
		bool quit;

		/*
		 * FoV
		 */
		float perspective_zoom;
		float zoom;
		float far_plane;
		float near_plane;

		float camera_sphere_view_rotation_alpha;
		float camera_sphere_view_rotation_beta;

		float common_scale;


		// parameter for model matrix
		EyeBall<float> modelEyeBall;

		// camera control
		Camera1stPerson<float> camera1stPerson;

		// player's velocity
		std::array<float,3> player_velocity;

		GlFreeType glFreeType;

		InputStateMouse inputStateMouse;

		Time time;

		// scenery lights
		Lights lights;

		// commonly used shader programs
		GlCommonShaderPrograms commonShaderPrograms;


		class Matrices
		{
		public:
			CMatrix4<float> model;
			CMatrix4<float> view;
			CMatrix4<float> projection;

			CMatrix4<float> view_model;
			CMatrix4<float> pvm;
		} matrices;


		EngineState()	:
			quit(false)
		{
		}


		void resetView()
		{
			modelEyeBall.reset();
			camera1stPerson.reset();

			common_scale = 2.0;
			camera1stPerson.setPosition(
					CVector<3,float>(
							0.0*common_scale,
							0.1*common_scale, //domain_visualization_size*0.5,
							1.1*common_scale //domain_visualization_size*1.1
						)
				);

			perspective_zoom = 1.0;

			far_plane = 10.0;
			near_plane = 0.01;

			camera_sphere_view_rotation_alpha = 0;
			camera_sphere_view_rotation_beta = 0;

			camera1stPerson.rotate(0, 0, 0);
		}
	};


	EngineState *engineState;


public:
	/**
	 * Event callbacks for RenderWindow
	 */
	class WindowEventCallbacks	:
		public RenderWindow::WindowEventCallbacks
	{
	public:
		ProgramCallbacks *programCallbacks;
		EngineState *engineState;
		RenderWindow *renderWindow;

		// last mouse position
		float old_mouse_x, old_mouse_y;


		WindowEventCallbacks(
				ProgramCallbacks *i_engineCallback,
				EngineState *i_engineState,
				RenderWindow *i_renderWindow
		)	:
			programCallbacks(i_engineCallback),
			engineState(i_engineState),
			renderWindow(i_renderWindow)
		{
			old_mouse_x = -1;
			old_mouse_y = -1;
		}




		void callback_key_down(int key, int mod, int scancode, int unicode)
		{
			if ((mod & KMOD_RSHIFT) || (mod & KMOD_LSHIFT))
				key += 'A'-'a';

			switch (key)
			{
			case 'q':
			case 'Q':
				engineState->quit = true;
				break;


			case 'a':
				engineState->player_velocity[0] = std::max<float>(engineState->player_velocity[0]-1, -1);	break;

			case 'd':
				engineState->player_velocity[0] = std::min<float>(engineState->player_velocity[0]+1, +1);	break;

			case 'w':
				engineState->player_velocity[2] = std::max<float>(engineState->player_velocity[2]-1, -1);	break;

			case 's':
				engineState->player_velocity[2] = std::min<float>(engineState->player_velocity[2]+1, +1);	break;

			case 'e':
				engineState->resetView();
				break;

			default:
				programCallbacks->vis_keypress(key);
				break;
			}
		}



		void callback_key_up(
				int key,
				int mod,
				int scancode,
				int unicode
		)
		{
			switch (key)
			{
				case 'a':	engineState->player_velocity[0] += 1;	break;
				case 'd':	engineState->player_velocity[0] -= 1;	break;
				case 'w':	engineState->player_velocity[2] += 1;	break;
				case 's':	engineState->player_velocity[2] -= 1;	break;
				default:
					break;
			}
		}


		void callback_quit()
		{
			engineState->quit = true;
		}


		void callback_viewport_changed(
				int i_width,
				int i_height
		)
		{
			engineState->glFreeType.viewportChanged(renderWindow->window_width, renderWindow->window_height);

			programCallbacks->vis_viewportChanged(i_width, i_height);
		}



		void callback_mouse_motion(
				int i_x,
				int i_y
		)
		{
			engineState->inputStateMouse.update((float)i_x*2.0/(float)renderWindow->window_width-1.0, (float)i_y*2.0/(float)renderWindow->window_height-1.0);

			if (old_mouse_x != -1)
			{
				if (engineState->inputStateMouse.mouse_buttons[RenderWindow::MOUSE_BUTTON_LEFT])
				{
					float b = engineState->zoom*10.0*engineState->common_scale;

					engineState->camera1stPerson.rotate(
							-engineState->inputStateMouse.relative_mouse_y*b,
							-engineState->inputStateMouse.relative_mouse_x*b,
							0
						);
				}
				else if (engineState->inputStateMouse.mouse_buttons[RenderWindow::MOUSE_BUTTON_RIGHT])
				{
					double speed = 0.1;
					engineState->modelEyeBall.rotate(speed*(-old_mouse_x + i_x), engineState->modelEyeBall.up);
//					engineState->modelEyeBall.rotate(speed*(-old_mouse_x + i_x), GLSL::vec3(0, 1, 0));
					engineState->modelEyeBall.rotate(speed*(-old_mouse_y + i_y), engineState->modelEyeBall.right);
				}
				else if (engineState->inputStateMouse.mouse_buttons[RenderWindow::MOUSE_BUTTON_MIDDLE])
				{
					engineState->perspective_zoom += (float)(engineState->inputStateMouse.relative_mouse_y);
				}
			}

			old_mouse_x = i_x;
			old_mouse_y = i_y;
		}

		void callback_mouse_button_down(int i_button_nr)
		{
			if (i_button_nr <= 3)
				engineState->inputStateMouse.mouse_buttons[i_button_nr] = true;
		}

		void callback_mouse_button_up(int i_button_nr)
		{
			if (i_button_nr <= 3)
				engineState->inputStateMouse.mouse_buttons[i_button_nr] = false;

			engineState->inputStateMouse.clearRelativeMovement();
		}

		void callback_mouse_wheel(int x, int y)
		{
			engineState->perspective_zoom += (float)y*(-0.1);
		}
	};

	WindowEventCallbacks *windowEventCallbacks;

	RenderWindow *renderWindow;


public:
	template <typename SimT>
	VisualizationEngine(
			SimT *i_simulationClassWithProgramCallbacks,
			const char *i_window_title
		)	:
			programCallbacks(&(ProgramCallbacks&)*i_simulationClassWithProgramCallbacks)
	{
		renderWindow = new RenderWindow(
				"SWEET",
				-1,
				-1
				);


		if (!renderWindow->initialized)
		{
			std::cerr << "Unable to initialize Render Window!" << std::endl;
			return;
		}

		/*
		 * event loop is so far not executed!
		 */


		/*
		 * create engine state (position, light, shaders, etc.)
		 */
		engineState = new EngineState;

		engineState->resetView();

		engineState->player_velocity[0] = 0;
		engineState->player_velocity[1] = 0;
		engineState->player_velocity[2] = 0;

		engineState->common_scale = 1;

		engineState->glFreeType.loadFont(16, true);
		engineState->glFreeType.viewportChanged(renderWindow->window_width, renderWindow->window_height);


		/*
		 * Allocate callbacks
		 */
		windowEventCallbacks = new WindowEventCallbacks(programCallbacks, engineState, renderWindow);


		/*
		 * Attach callbacks to rendering window
		 * NOW WE ARE ABLE TO PROCESS EVENTS!
		 */
		renderWindow->setWindowEventCallback(windowEventCallbacks);

		programCallbacks->vis_setup(this);

		while (!engineState->quit)
		{
			/*
			 * EVENTS
			 */
			renderWindow->eventLoop();

			glClearColor(0,0,0,0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			glEnable(GL_CULL_FACE);
			glEnable(GL_DEPTH_TEST);


			engineState->time.update();

			engineState->inputStateMouse.clearRelativeMovement();


			/*
			 * MOVEMENT
			 */
			float a = engineState->common_scale*(float)engineState->time.frame_elapsed_seconds*(float)1.5;
			CVector<3,float> movement(
					engineState->player_velocity[0]*a,
					engineState->player_velocity[1]*a,
					engineState->player_velocity[2]*a
			);

			engineState->camera1stPerson.moveRelative(movement);

			engineState->zoom = std::exp(engineState->perspective_zoom)*0.2;


			/*
			 * VISUALIZATION
			 */
			engineState->camera1stPerson.computeMatrices();

			float fzoom = engineState->zoom*engineState->near_plane;

			engineState->camera1stPerson.frustum(-fzoom*renderWindow->aspect_ratio, fzoom*renderWindow->aspect_ratio, -fzoom, fzoom, engineState->near_plane, engineState->far_plane);


			/*
			 * MATRICES
			 */
			engineState->modelEyeBall.reconstructRotationMatrix();
			engineState->matrices.model = engineState->modelEyeBall.rotationMatrix;

			engineState->matrices.view = engineState->camera1stPerson.view_matrix;
			engineState->matrices.projection = engineState->camera1stPerson.projection_matrix;

			engineState->matrices.view_model = engineState->matrices.view * engineState->matrices.model;
			engineState->matrices.pvm = engineState->matrices.projection * engineState->matrices.view_model;


			/*
			 * LIGHTING
			 */
			GLSL::vec4 v = engineState->camera1stPerson.view_matrix
					*GLSL::vec3(
							engineState->lights.lights[0].world_pos3[0],
							engineState->lights.lights[0].world_pos3[1],
							engineState->lights.lights[0].world_pos3[2]
					);

			engineState->lights.lights[0].enabled = true;
			engineState->lights.lights[0].light_view_pos3 = {{v[0], v[1], v[2]}};


			engineState->commonShaderPrograms.shaderBlinn.use();
			engineState->commonShaderPrograms.shaderBlinn.texture0_enabled.set1b(false);
			engineState->commonShaderPrograms.shaderBlinn.setupUniformsLights(engineState->lights);
			CGlErrorCheck();

			engineState->commonShaderPrograms.shaderBlinn.view_model_matrix_uniform.set(engineState->matrices.view_model);
			engineState->commonShaderPrograms.shaderBlinn.view_model_normal_matrix3_uniform.set(engineState->matrices.view_model.getInverseTranspose3x3());
			engineState->commonShaderPrograms.shaderBlinn.pvm_matrix_uniform.set(engineState->matrices.pvm);
			CGlErrorCheck();

			engineState->commonShaderPrograms.shaderBlinn.material_ambient_color3_uniform.set(GLSL::vec3(0.1,0.1,0.2));
			engineState->commonShaderPrograms.shaderBlinn.material_diffuse_color3_uniform.set(GLSL::vec3(0.1,0.1,1.0));
			engineState->commonShaderPrograms.shaderBlinn.material_specular_color3_uniform.set(GLSL::vec3(1.0,1.0,1.0));
			engineState->commonShaderPrograms.shaderBlinn.material_specular_exponent_uniform.set1f(40.0f);
			engineState->commonShaderPrograms.shaderBlinn.disable();

			programCallbacks->vis_render();

			engineState->quit |= programCallbacks->should_quit();

			renderWindow->setWindowTitle(programCallbacks->vis_getStatusString());

			renderWindow->swapBuffer();
		}

		programCallbacks->vis_shutdown();

		delete engineState;
		delete windowEventCallbacks;
		delete renderWindow;
	}


};

#endif /* SRC_INCLUDE_LIBGL_ENGINE_ENGINE_HPP_ */
