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

#ifndef __I_INPUT_STATE_MOUSE_HPP__
#define __I_INPUT_STATE_MOUSE_HPP__


/**
 * input state stores the current state about the mouse:
 *
 * its position, relative mouse movements during the last frames
 * and pressed buttons
 */
class InputStateMouse
{
public:
	bool valid;

	// mouse position in window in [-1;1] coordinates
	float mouse_x, mouse_y;

	// relative mouse movements
	float relative_mouse_x, relative_mouse_y;

	/**
	 * state of mouse buttons
	 */
	bool mouse_buttons[5];

	InputStateMouse()
	{
		mouse_x = mouse_y = 0;
		relative_mouse_x = relative_mouse_y = 0;

		for (int i = 0; i < 5; i++)
			mouse_buttons[i] = false;

		valid = false;
	}

	virtual ~InputStateMouse()
	{

	}

	void update(float x, float y)
	{
		relative_mouse_x = mouse_x - x;
		relative_mouse_y = mouse_y - y;

		mouse_x = x;
		mouse_y = y;
	}
	void clearRelativeMovement()
	{
		relative_mouse_x = 0;
		relative_mouse_y = 0;
	}
};

#endif
