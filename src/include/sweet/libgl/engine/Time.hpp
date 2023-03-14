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

#ifndef C_TIME_HPP__
#define C_TIME_HPP__

/**
 * simple method to get the time in seconds, fps and elapsed frame seconds
 * from the system
 */
class Time
{
private:
	double next_fps_time_update;
	double last_fps_time_update;
	double fps_frames;

public:
	// absolute seconds
	double elapsed_seconds;

	// seconds from last to current frame
	double frame_elapsed_seconds;

	// Simply the FPS
	double fps;

	bool fpsUpdatedInLastFrame;

	Time()
	{
		fpsUpdatedInLastFrame = false;
		elapsed_seconds = (double)SDL_GetTicks()*0.001;
		frame_elapsed_seconds = 0;
		fps = 0;
		fps_frames = 0;
		next_fps_time_update = 0;
		last_fps_time_update = 0;
	}

	~Time()
	{
	}

	void update()
	{
		double new_seconds = (double)SDL_GetTicks()*(double)0.001;

		frame_elapsed_seconds = new_seconds - elapsed_seconds;

		elapsed_seconds = new_seconds;

		fps_frames++;

		if (next_fps_time_update < elapsed_seconds)
		{
			// compute fps
			fps = fps_frames/(elapsed_seconds-last_fps_time_update);

			// set next fps time update for every second
			next_fps_time_update += 1.0;

			// if the additional second of the last step is already behind (in this case the program runs really slow), use current time value+1
			if (next_fps_time_update < elapsed_seconds)
				next_fps_time_update = elapsed_seconds+1.0;

			last_fps_time_update = elapsed_seconds;

			fps_frames = 0;

			fpsUpdatedInLastFrame = true;
		}
		else
		{
			fpsUpdatedInLastFrame = false;
		}
	}
};


#endif //__I_TIME_HPP__
