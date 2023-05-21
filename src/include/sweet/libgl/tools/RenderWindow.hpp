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
 * info: OpenGL 3.2 initialization code taken from khronos website
 */

#ifndef RENDER_WINDOW_HPP_
#define RENDER_WINDOW_HPP_

#include <sweet/libgl/incgl3.h>

extern "C" {
	#include <SDL2/SDL.h>
	#include <SDL2/SDL_thread.h>
}


/*
 * std has to be included before SDL!!!
 */
#include <cmath>

#include <iostream>
#include <sstream>
#include <signal.h>

#include <sweet/libgl/tools/Bitmap.hpp>



/**
 * \brief initialize gui and rendering context
 */
class RenderWindow
{
	/**
	 * \brief callback handlers for gui events
	 */
public:
	class WindowEventCallbacks
	{
	public:
		virtual void callback_key_down(int key, int mod, int scancode, int unicode) = 0;	///< key down event
		virtual void callback_key_up(int key, int mod, int scancode, int unicode) = 0;		///< key up event
		virtual void callback_quit() = 0;							///< quit event
		virtual void callback_mouse_motion(int x, int y) = 0;		///< mouse position in absolute coordinates
		virtual void callback_mouse_button_down(int button) = 0;	///< button press event
		virtual void callback_mouse_button_up(int button) = 0;		///< button release event
		virtual void callback_viewport_changed(int width, int height) = 0;	///< viewport changed (window resized)
		virtual void callback_mouse_wheel(int x, int y) = 0;		///< viewport changed (window resized)

		virtual ~WindowEventCallbacks()
		{
		}
	};




	bool initialized;

	WindowEventCallbacks *windowEventCallbacksImplementation;

	SDL_Thread *save_bitmap_thread;
	std::string save_bitmap_filename;	///< bitmap filename when storing bitmaps using threads
	Bitmap24 save_bitmap;


	SDL_Window *window;
	SDL_GLContext glContext;


public:
	bool fullscreen_active;		///< true, if fullscreen is active

	int window_width;			///< width of viewport/window
	int window_height;			///< height of viewport/window
	float aspect_ratio;			///< current aspect ratio (window/height)

	enum
	{
		MOUSE_BUTTON_LEFT = SDL_BUTTON_LEFT,
		MOUSE_BUTTON_RIGHT = SDL_BUTTON_RIGHT,
		MOUSE_BUTTON_MIDDLE = SDL_BUTTON_MIDDLE
	};

	enum
	{
		KEY_PAGEUP = SDLK_PAGEUP,
		KEY_PAGEDOWN = SDLK_PAGEDOWN,
		KEY_F1 = SDLK_F1,
		KEY_F2 = SDLK_F2,
		KEY_F3 = SDLK_F3,
		KEY_F4 = SDLK_F4,
		KEY_F5 = SDLK_F5,
		KEY_F6 = SDLK_F6,
		KEY_F7 = SDLK_F7,
		KEY_F8 = SDLK_F8,
		KEY_F9 = SDLK_F9,
		KEY_F10 = SDLK_F10,
		KEY_F11 = SDLK_F11,
		KEY_F12 = SDLK_F12,

#if WIN32
		KEY_MOD_SHIFT = KMOD_LSHIFT
#else
		KEY_MOD_SHIFT = (KMOD_LSHIFT | KMOD_RSHIFT)
#endif
	};

	/**
	 * update viewport with current width and height values and update apsect ratio
	 */
	void updateViewport()
	{
		glViewport(0, 0, window_width, window_height);
		aspect_ratio = (float)window_width / (float)window_height;
	}


	/**
	 * set fullscreen mode or deactivate
	 */
	void setWindowFullscreenState(
			bool i_fullscreen		///< set to true, if fullscreen should be activated - false otherwise
	)
	{
		if (i_fullscreen)
		{
			if (SDL_SetWindowFullscreen(window, SDL_TRUE) == -1)
			{
				std::cerr << "Failed to initialize fullscreen mode: " << SDL_GetError();
				std::cerr << std::endl;
				return;
			}
		}
		else
		{
			if (SDL_SetWindowFullscreen(window, SDL_FALSE) == -1)
			{
				std::cerr << "Failed to initialize window mode: ";
				std::cerr << SDL_GetError();
				std::cerr << std::endl;
				return;
			}
		}

		fullscreen_active = i_fullscreen;
	}


	void checkSDLError(int line = -1)
	{
		const char *error = SDL_GetError();
		if (*error != '\0')
		{
			printf("SDL Error: %s\n", error);
//			if (line != -1)
//				printf(" + line: %i\n", line);
			SDL_ClearError();
			exit(1);
		}
	}


	/**
	 * setup an opengl render context
	 */
	bool setupOpenGLRenderContext(
			const std::string &i_initial_window_title,
			int i_initial_window_width = 800,
			int i_initial_window_height = 600,
			int i_request_opengl_major_version = 3,		///< major version of opengl context to request
			int i_request_opengl_minor_version = 0		///< minor version of opengl context to request
	)
	{
		if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_NOPARACHUTE) < 0) /* Initialize SDL's Video subsystem */
		{
			std::cerr << "Unable to initialize SDL: "; std::cerr << SDL_GetError();
			std::cerr << std::endl;
			exit(-1);
			return false;
		}

#if !WIN32
		// activate ctrl-c abort handling
		signal(SIGINT, SIG_DFL);
		signal(SIGQUIT, SIG_DFL);
#endif

		checkSDLError(__LINE__);

		/*
		 * set the OpenGL version number to create the context for (p_request_opengl_major_version, p_request_opengl_minor_version)
		 */
		if (SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, i_request_opengl_major_version) < 0)
		{
			std::cerr << "unable to set major OpenGL version to " << i_request_opengl_major_version << ": " << SDL_GetError() << std::endl;
			return false;
		}

		if (SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, i_request_opengl_minor_version) < 0)
		{
			std::cerr << "unable to set minor OpenGL version to " << i_request_opengl_minor_version << ": " << SDL_GetError() << std::endl;
			return false;
		}

		// Specifically request core mode for MacOSX systems
		if (SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE) < 0)
		{
			std::cerr << "unable to set CORE profile context" << i_request_opengl_minor_version << std::endl;
			return false;
		}
		/*
		 * Turn on double buffering with a 24bit Z buffer.
		 * You may need to change this to 16 or 32 for your system.
		 */
		SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
		SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

		// disable vsync
//		SDL_GL_SetAttribute(SDL_GL_SWAP_CONTROL, 1);

		checkSDLError(__LINE__);

//		std::cout << "Get current video driver: " << (char*)SDL_GetCurrentVideoDriver() << std::endl;

		if (i_initial_window_width < 0)
		{
			float ddpi;
			SDL_GetDisplayDPI(0, &ddpi, nullptr, nullptr);

			i_initial_window_width = 2.0*800.0*ddpi/150.0;
			i_initial_window_height = 2.0*600.0*ddpi/150.0;

			std::cout << ddpi << std::endl;
			std::cout << i_initial_window_width << ", " << i_initial_window_height << std::endl;
		}

		window = SDL_CreateWindow(
						i_initial_window_title.c_str(),
						SDL_WINDOWPOS_CENTERED,
						SDL_WINDOWPOS_CENTERED,
						i_initial_window_width,
						i_initial_window_height,
//						(fullscreen_active ?
//								(SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_FULLSCREEN | SDL_WINDOW_BORDERLESS) :
								(SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE)
//						)
			);

		checkSDLError(__LINE__);

		if (window == 0) /* Die if creation failed */
		{
			std::cerr << "Unable to create window: " << SDL_GetError();
			std::cerr << "Requested OpenGL core profile with version " << i_request_opengl_major_version << "." << i_request_opengl_minor_version << std::endl;
			std::cerr << std::endl;
			return false;
		}


		/*
		 * Create our OpenGL context and attach it to our window
		 */
		glContext = SDL_GL_CreateContext(window);
		if (!glContext)
		{
			std::cerr << "Unable to create GL context, SDL error message: " << SDL_GetError() << std::endl;
			std::cerr << "Requested OpenGL core profile with version " << i_request_opengl_major_version << "." << i_request_opengl_minor_version << std::endl;
			std::cerr << std::endl;
			SDL_DestroyWindow(window);
			return false;
		}

		checkSDLError(__LINE__);

		/*
		 * Set update intervals to 0 for immediate updates without caring about vsyncs
		 */
		SDL_GL_SetSwapInterval(0);


		if (fullscreen_active)
		{
			setWindowFullscreenState(true);
		}

		checkSDLError(__LINE__);

#if SWEET_DEBUG
		std::cout << "GL Version: " << glGetString(GL_VERSION) << std::endl;
		std::cout << "GL Vendor: " << glGetString(GL_VENDOR) << std::endl;
		std::cout << "GL Renderer: " << glGetString(GL_RENDERER) << std::endl;
		std::cout << "GL Shading language version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
#endif

		window_width = i_initial_window_width;
		window_height = i_initial_window_height;
		updateViewport();
		return true;
	}



public:
	/**
	 * initialize GUI
	 */
	RenderWindow(
			const std::string &i_initial_window_title,					///< title of window to display in the window title bar
			int p_window_width = 800,					///< initial width of window in window mode
			int p_window_height = 600,					///< initial height of window in window mode
			bool p_fullscreen = false,					///< set to true to enable fullscreen mode on window creation
			int p_request_opengl_major_version = 3,		///< major version of opengl context to request
			int p_request_opengl_minor_version = 1		///< minor version of opengl context to request
	)
	{
		windowEventCallbacksImplementation = nullptr;
		save_bitmap_thread = nullptr;
		fullscreen_active = p_fullscreen;

		initialized = setupOpenGLRenderContext(
				i_initial_window_title,
				p_window_width,
				p_window_height,
				p_request_opengl_major_version,
				p_request_opengl_minor_version
			);
	}


	/**
	 * set the window event callback
	 */
	void setWindowEventCallback(
			class WindowEventCallbacks *i_windowEventCallbacksImplementation
	)
	{
		windowEventCallbacksImplementation = i_windowEventCallbacksImplementation;
	}

	virtual ~RenderWindow()
	{
		if (initialized)
		{
			// delete context
			SDL_GL_DeleteContext(glContext);

			// destroy window
			SDL_DestroyWindow(window);

			// quit SDL
			SDL_Quit();
		}
	}



	/**
	 * set the window title
	 */
	void setWindowTitle(
			const std::string &i_title	// string of title
	)
	{
		SDL_SetWindowTitle(window, i_title.c_str());
	}



	/**
	 * this function has to be called by the rendering loop on every iteration to catch the events and
	 * forward them to the class implementing the callback functions CGuiInterface
	 */
	void eventLoop() 
	{
		SDL_Event event;
		// process pending events
		while(SDL_PollEvent(&event))
		{
			switch(event.type)
			{
				case SDL_WINDOWEVENT:

					switch (event.window.event)
					{
						case SDL_WINDOWEVENT_CLOSE:
							windowEventCallbacksImplementation->callback_quit();
							break;

						case SDL_WINDOWEVENT_RESIZED:
						case SDL_WINDOWEVENT_SIZE_CHANGED:
							window_width = event.window.data1;
							window_height = event.window.data2;

							if (window_width <= 0)	window_width = 1;
							if (window_height <= 0)	window_height = 1;

							updateViewport();
							windowEventCallbacksImplementation->callback_viewport_changed(window_width, window_height);
							break;


						case SDL_WINDOWEVENT_SHOWN:
							SDL_GetWindowSize(window, &window_width, &window_height);
							updateViewport();
							windowEventCallbacksImplementation->callback_viewport_changed(window_width, window_height);
							break;
					}
					break;

				case SDL_KEYDOWN:
					windowEventCallbacksImplementation->callback_key_down(event.key.keysym.sym, event.key.keysym.mod, event.key.keysym.scancode, 0);
					break;

				case SDL_KEYUP:
					windowEventCallbacksImplementation->callback_key_up(event.key.keysym.sym, event.key.keysym.mod, event.key.keysym.scancode, 0);
					break;

				case SDL_QUIT:
					windowEventCallbacksImplementation->callback_quit();
					break;

				case SDL_MOUSEMOTION:
					windowEventCallbacksImplementation->callback_mouse_motion(event.motion.x, event.motion.y);
					break;

				case SDL_MOUSEBUTTONDOWN:
					switch (event.button.button)
					{
					case 1:
					case 2:
					case 3:
						windowEventCallbacksImplementation->callback_mouse_button_down(event.button.button);
						break;

					case 4:
						windowEventCallbacksImplementation->callback_mouse_wheel(0, -1);
						break;

					case 5:
						windowEventCallbacksImplementation->callback_mouse_wheel(0, 1);
						break;
					}
					break;

				case SDL_MOUSEBUTTONUP:
					windowEventCallbacksImplementation->callback_mouse_button_up(event.button.button);
					break;

				case SDL_MOUSEWHEEL:
					windowEventCallbacksImplementation->callback_mouse_wheel(event.wheel.x, event.wheel.y);
					break;
			}
		}
	}

	/**
	 * return the ticks as an integer number (DEPRECATED)
	 */
	signed int getTicksInt()
	{
		return SDL_GetTicks();
	}

	/**
	 * return the ticks in seconds
	 */
	double getTicks()
	{
		return (double)SDL_GetTicks()*(1.0/1000.0);
	}

	/**
	 * save a screenshot to a bitmap file
	 */
	static
	bool saveScreenshot(
			const std::string &i_filename	///< filepath to store the screenshot to
	)
	{
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);

		int width = viewport[2];
		int height = viewport[3];

		Bitmap24 bitmap(width, height);

		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glPixelStorei(GL_PACK_ROW_LENGTH, 0);
		glPixelStorei(GL_PACK_IMAGE_HEIGHT, 0);
		glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, bitmap.data.data());

		if (!bitmap.save(i_filename))
		{
			std::cerr << "Failed to save bitmap" << std::endl;
			return false;
		}

		return true;
	}


private:
	/**
	 * thread which saves the screenshot data
	 */
	static int saveScreenshotThread(void *user_data)
	{
		RenderWindow &g = *(RenderWindow*)user_data;

		if (!g.save_bitmap.save(g.save_bitmap_filename))
		{
			std::cerr << "ERROR during saving screenshot" << std::endl;
			return false;
		}
		return true;
	}

public:
	/**
	 * save a screenshot to a bitmap file starting a thread after loading the screenshot data from framebuffer
	 */
	bool saveScreenshotWithThread(
			const std::string &filename	///< filepath to store the screenshot to
	)
	{
		if (save_bitmap_thread != nullptr)
		{
			int status;
			SDL_WaitThread(save_bitmap_thread, &status);
		}

		save_bitmap.resize(window_width, window_height);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glPixelStorei(GL_PACK_ROW_LENGTH, 0);
		glPixelStorei(GL_PACK_IMAGE_HEIGHT, 0);
		glReadPixels(0, 0, window_width, window_height, GL_BGR, GL_UNSIGNED_BYTE, save_bitmap.data.data());

		save_bitmap_filename = filename;
		save_bitmap_thread = SDL_CreateThread(&saveScreenshotThread, "saveBitmapThread", this);

		if (save_bitmap_thread == nullptr)
		{
			std::cout << "failed to create thread to store bitmap" << std::endl;
			return false;
		}

		return true;
	}

	/**
	 * this function has to be called on the end of the rendering loop to swap the front with back buffer
	 */
	void swapBuffer()
	{
		SDL_GL_SwapWindow(window);
	}
};

#endif
