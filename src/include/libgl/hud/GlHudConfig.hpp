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
 * GlHudConfig.hpp
 *
 *  Created on: Mar 22, 2010
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef CGLHUDCONFIG_HPP_
#define CGLHUDCONFIG_HPP_

#include <list>

#include <cmath>
#include <libgl/hud/GlWindow.hpp>
#include <libgl/hud/GlHudConfig.hpp>
#include <libgl/hud/GlFreeType.hpp>
#include <libgl/hud/GlRenderOStream.hpp>
#include <libgl/tools/RenderWindow.hpp>



/**
 * HudConfig cares about config variables which can be modified via a simple gui
 *
 * the config variables are not stored in this class. instead, a pointer to the storage of each config
 * variable is given as a parameter on setup.
 */
class GlHudConfig
{
public:
	/**
	 * enumeration for different types of configuration values
	 */
	enum ConfigType
	{
		CONFIG_INT = 0,
		CONFIG_UINT,
		CONFIG_LONGLONG,
		CONFIG_ULONGLONG,
		CONFIG_FLOAT,
		CONFIG_FLOAT_HI,
		CONFIG_DOUBLE,
		CONFIG_DOUBLE_HI,
		CONFIG_BOOLEAN,
		CONFIG_LINEBREAK,	///< insert a linebreak (empty line)
		CONFIG_TEXT			///< insert some text
	};

	/**
	 * configuration for integer option values
	 */
	class ConfigInt
	{
	public:
		int *value;
		int min;
		int max;
		int delta;

		void up(int c = 1)
		{
			*value = std::min<int>(max, *value + delta*c);
		}

		void down(int c = 1)
		{
			*value = std::max<int>(min, *value - delta*c);
		}
	};


	/**
	 * configuration for unsigned integer option values
	 */
	class ConfigUInt
	{
	public:
		unsigned int *value;
		unsigned int min;
		unsigned int max;
		unsigned int delta;

		void up(unsigned int c = 1)
		{
			*value = std::min<unsigned int>(max, *value + delta*c);
		}

		void down(unsigned int c = 1)
		{
			if (*value < delta*c)
			{
				*value = 0;
				return;
			}

			*value = std::max<unsigned int>(min, *value - delta*c);
		}
	};


	/**
	 * configuration for unsigned integer option values
	 */
	class ConfigULongLong
	{
	public:
		unsigned long long *value;
		unsigned long long min;
		unsigned long long max;
		unsigned long long delta;

		void up(unsigned long long c = 1)
		{
			*value = std::min<unsigned long long>(max, *value + delta*c);
		}

		void down(unsigned long long c = 1)
		{
			if (*value < delta*c)
			{
				*value = 0;
				return;
			}

			*value = std::max<unsigned long long>(min, *value - delta*c);
		}
	};



	/**
	 * configuration for unsigned integer option values
	 */
	class ConfigLongLong
	{
	public:
		long long *value;
		long long min;
		long long max;
		long long delta;

		void up(long long c = 1)
		{
			*value = std::min<long long>(max, *value + delta*c);
		}

		void down(long long c = 1)
		{
			if (*value < delta*c)
			{
				*value = 0;
				return;
			}

			*value = std::max<long long>(min, *value - delta*c);
		}
	};


	/**
	 * configuration for float option values
	 */
	class ConfigFloat
	{
	public:
		float *value;
		float min;
		float max;
		float delta;

		void up(int c = 1)
		{
			*value = std::min<float>(max, *value + delta*(float)c);
		}

		void down(int c = 1)
		{
			*value = std::max<float>(min, *value - delta*(float)c);
		}
	};



	/**
	 * configuration for float option values
	 */
	class ConfigDouble
	{
	public:
		double *value;
		double min;
		double max;
		double delta;

		void up(int c = 1)
		{
			*value = std::min<double>(max, *value + delta*(double)c);
		}

		void down(int c = 1)
		{
			*value = std::max<double>(min, *value - delta*(double)c);
		}
	};

	/**
	 * configuration for float option values with half increment (HI)
	 *
	 * the larger the value, the larger is the change
	 *
	 * for simplicity, the change is half of the current variables value
	 */
	class ConfigFloatHI
	{
	public:
		float *value;
		float min;
		float max;
		float delta;

		void up(int c = 1)
		{
			for (int i = 0; i < c; i++)
			{
				*value = std::min<float>(max, *value*(1.0+0.5*delta));
			}
		}

		void down(int c = 1)
		{
			for (int i = 0; i < c; i++)
			{
				*value = std::min<float>(max, *value*(1.0-0.5*delta));
			}
		}
	};

	/**
	 * configuration for float option values with half increment (HI)
	 *
	 * the larger the value, the larger is the change
	 *
	 * for simplicity, the change is half of the current variables value
	 */
	class ConfigDoubleHI
	{
	public:
		double *value;
		double min;
		double max;
		double delta;

		void up(int c = 1)
		{
			for (int i = 0; i < c; i++)
			{
				*value = std::min<double>(max, *value*(1.0+0.5*delta));
			}
		}

		void down(int c = 1)
		{
			for (int i = 0; i < c; i++)
			{
				*value = std::min<double>(max, *value*(1.0-0.5*delta));
			}
		}
	};

	/**
	 * configuration for boolean option values
	 */
	class ConfigBoolean
	{
	public:
		bool *value;

		void change()
		{
			*value = !*value;
		}
	};

	/**
	 * option handler
	 */
	class COption
	{
	public:
		std::string description;
		ConfigType type;
		int user_value;

		// private handlers for rendering / gui interactions
		int id;			///< unique id
		bool in_area;	///< true, if mouse is within the option area
		bool button_left_down;	///< true, if left mouse button is pressed (but not yet released)
		bool button_right_down;	///< true, if right mouse button is pressed (but not yet released)

		int area_left;
		int area_right;
		int area_top;
		int area_bottom;

		int description_render_left;
		int description_render_top;

		int value_render_left;
		int value_render_top;

		// callback handler which is called, when the value was changed
		void (*callback_changed)(void *user_ptr);
		void *callback_user_ptr;

		COption()	:
			area_left(0),
			area_right(0),
			area_top(0),
			area_bottom(0)
		{

		}

		union
		{
			ConfigInt configInt;
			ConfigUInt configUInt;
			ConfigLongLong configLongLong;
			ConfigULongLong configULongLong;
			ConfigFloat configFloat;
			ConfigFloatHI configFloatHI;
			ConfigDouble configDouble;
			ConfigDoubleHI configDoubleHI;
			ConfigBoolean configBoolean;
		};

		COption &setupBoolean(	const char *p_description,
								bool *p_value,
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_BOOLEAN;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configBoolean.value = p_value;

			callback_changed = nullptr;
			return *this;
		}


		/**
		 * single precision
		 */
		COption &setupFloat(	const char *p_description,
								float *p_value,
								float p_delta = 0.1,
								float p_min = -std::numeric_limits<float>::infinity(),
								float p_max = std::numeric_limits<float>::infinity(),
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_FLOAT;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configFloat.value = p_value;
			configFloat.delta = p_delta;
			configFloat.min = p_min;
			configFloat.max = p_max;

			callback_changed = nullptr;
			return *this;
		}



		/**
		 * double precision
		 */
		COption &setupFloat(	const char *p_description,
								double *p_value,
								double p_delta = 0.1,
								double p_min = -std::numeric_limits<double>::infinity(),
								double p_max = std::numeric_limits<double>::infinity(),
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_DOUBLE;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configDouble.value = p_value;
			configDouble.delta = p_delta;
			configDouble.min = p_min;
			configDouble.max = p_max;

			callback_changed = nullptr;
			return *this;
		}


		/**
		 * single precision
		 */
		COption &setupFloatHI(	const char *p_description,
								float *p_value,
								float p_delta = 0.1,
								float p_min = -std::numeric_limits<float>::infinity(),
								float p_max = std::numeric_limits<float>::infinity(),
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_FLOAT_HI;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configFloat.value = p_value;
			configFloat.delta = p_delta;
			configFloat.min = p_min;
			configFloat.max = p_max;

			callback_changed = nullptr;
			return *this;
		}


		/**
		 * double precision
		 */
		COption &setupFloatHI(	const char *p_description,
								double *p_value,
								double p_delta = 0.1,
								double p_min = -std::numeric_limits<float>::infinity(),
								double p_max = std::numeric_limits<float>::infinity(),
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_DOUBLE_HI;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configDouble.value = p_value;
			configDouble.delta = p_delta;
			configDouble.min = p_min;
			configDouble.max = p_max;

			callback_changed = nullptr;
			return *this;
		}


		COption &setupInt(		const char *p_description,
								int *p_value,
								int p_delta = 1,
								int p_min = -std::numeric_limits<int>::max(),
								int p_max = std::numeric_limits<int>::max(),
								int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_INT;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configInt.value = p_value;
			configInt.delta = p_delta;
			configInt.min = p_min;
			configInt.max = p_max;

			callback_changed = nullptr;
			return *this;
		}

		COption &setupUInt(	const char *p_description,
							unsigned int *p_value,
							unsigned int p_delta = 1,
							unsigned int p_min = -std::numeric_limits<int>::max(),
							unsigned int p_max = std::numeric_limits<int>::max(),
							unsigned int p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_UINT;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configUInt.value = p_value;
			configUInt.delta = p_delta;
			configUInt.min = p_min;
			configUInt.max = p_max;

			callback_changed = nullptr;
			return *this;
		}

		COption &setupULongLong(
							const char *p_description,
							unsigned long long *p_value,
							unsigned long long p_delta = 1,
							unsigned long long p_min = -std::numeric_limits<unsigned long long>::max(),
							unsigned long long p_max = std::numeric_limits<unsigned long long>::max(),
							unsigned long long p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_ULONGLONG;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configULongLong.value = p_value;
			configULongLong.delta = p_delta;
			configULongLong.min = p_min;
			configULongLong.max = p_max;

			callback_changed = nullptr;
			return *this;
		}

		COption &setupLongLong(
							const char *p_description,
							long long *p_value,
							long long p_delta = 1,
							long long p_min = -std::numeric_limits<long long>::max(),
							long long p_max = std::numeric_limits<long long>::max(),
							long long p_user_value = 0)
		{
			description = p_description;
			type = CONFIG_LONGLONG;
			user_value = p_user_value;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			configLongLong.value = p_value;
			configLongLong.delta = p_delta;
			configLongLong.min = p_min;
			configLongLong.max = p_max;

			callback_changed = nullptr;
			return *this;
		}

		COption &setupLinebreak()
		{
			description = "";
			type = CONFIG_LINEBREAK;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			callback_changed = nullptr;
			return *this;
		}


		COption &setupText(const char *p_description)
		{
			description = p_description;
			type = CONFIG_TEXT;
			in_area = false;
			button_left_down = false;
			button_right_down = false;

			callback_changed = nullptr;
			return *this;
		}
	};
private:
	std::list<COption> option_list;	///< list with all options

	bool visible;		///< true, if hud is rendered

	int option_id_counter;	///< incrementing counter to set unique ids for the options

	/**
	 * pointer to currently activated option for switching.
	 * when this pointer is valid and another option button is clicked with the left mouse button,
	 * both values are changes. this is useful for comparisons.
	 */
	COption *switch_option;

public:
	COption *active_option;	///< pointer to currently activated or highlighted option

	int area_left;		/// overall configuration area for mouse interactions
	int area_right;
	int area_top;
	int area_bottom;

	int area_width;
	int area_height;

	GlHudConfig()
	{
		visible = true;
		option_id_counter = 0;
		active_option = nullptr;
		switch_option = nullptr;
	}

	/**
	 * hide hud configuration
	 */
	void hide()
	{
		visible = false;
	}

	/**
	 * shot hud configuration
	 */
	void show()
	{
		visible = true;
	}

	/**
	 * setup a callback handler which is called, when the value changed
	 */
	void set_callback(	void *value_ptr,							///< pointer to value
						void (*callback_handler)(void *user_ptr),	///< callback handler
						void *user_ptr
					)
	{
		std::list<COption>::iterator i;
		bool found = false;
		for (i = option_list.begin(); i != option_list.end(); i++)
		{
			COption &o = *i;
			switch(o.type)
			{
				case CONFIG_BOOLEAN:
					if (o.configBoolean.value == value_ptr)
						found = true;
					break;

				case CONFIG_INT:
					if (o.configInt.value == value_ptr)
						found = true;
					break;

				case CONFIG_UINT:
					if (o.configUInt.value == value_ptr)
						found = true;
					break;

				case CONFIG_LONGLONG:
					if (o.configLongLong.value == value_ptr)
						found = true;
					break;

				case CONFIG_ULONGLONG:
					if (o.configULongLong.value == value_ptr)
						found = true;
					break;

				case CONFIG_FLOAT:
					if (o.configFloat.value == value_ptr)
						found = true;
					break;

				case CONFIG_FLOAT_HI:
					if (o.configFloat.value == value_ptr)
						found = true;
					break;

				case CONFIG_DOUBLE:
					if (o.configDouble.value == value_ptr)
						found = true;
					break;

				case CONFIG_DOUBLE_HI:
					if (o.configDouble.value == value_ptr)
						found = true;
					break;

				default:
					break;
			}

			if (found)
				break;
		}

		if (!found)
			return;

		(*i).callback_changed = callback_handler;
		(*i).callback_user_ptr = user_ptr;
	}

	/**
	 * insert new option into configuration list
	 */
	void insert(COption &p_option)
	{
		p_option.id = option_id_counter;
		option_list.push_back(p_option);
		option_id_counter++;
	}

	GlFreeType *cGlFreeType;
	GlRenderOStream *cGlRenderOStream;

	/**
	 * setup area_* variables to speedup rendering and mouse interactions
	 */
	void setup(
		GlFreeType &io_cGlFreeType,
		GlRenderOStream &io_cGlRenderOStream,
		int i_pos_x,	// x position to start rendering
		int i_pos_y		// y position to start rendering - this is the left upper corner of the first character!!!
	)
	{
		cGlFreeType = &io_cGlFreeType;
		cGlRenderOStream = &io_cGlRenderOStream;

		// get maximum text length
		int max_text_length = 0;
		for (std::list<COption>::iterator i = option_list.begin(); i != option_list.end(); i++)
			max_text_length = std::max<int>(cGlFreeType->getTextLength((*i).description.c_str()), max_text_length);

		int description_value_spacing = cGlFreeType->font_size/2;

		area_left = 0xfffffff;
		area_right = -0xfffffff;
		area_bottom = 0xfffffff;
		area_top = -0xfffffff;

		for (std::list<COption>::iterator i = option_list.begin(); i != option_list.end(); i++)
		{
			COption &o = *i;

			int text_length = cGlFreeType->getTextLength(o.description.c_str());
			o.area_left = i_pos_x + max_text_length - text_length;
			o.area_top = i_pos_y;

			// relative (window coords) render position of description text
			o.description_render_left = i_pos_x + max_text_length - text_length;
			o.description_render_top = i_pos_y - cGlFreeType->font_size;

			// relative position to render value
			o.value_render_left = i_pos_x + max_text_length + description_value_spacing;
			o.value_render_top = i_pos_y - cGlFreeType->font_size;

			int advance_x = 0;
			// render values / control
			switch(o.type)
			{
				case CONFIG_BOOLEAN:
				case CONFIG_INT:
				case CONFIG_UINT:
				case CONFIG_LONGLONG:
				case CONFIG_ULONGLONG:
				case CONFIG_FLOAT:
				case CONFIG_FLOAT_HI:
				case CONFIG_DOUBLE:
				case CONFIG_DOUBLE_HI:
					advance_x = cGlFreeType->font_size*4;
					break;

				default:
					break;
			}

			o.area_right = i_pos_x + max_text_length + advance_x;
			o.area_bottom = i_pos_y - cGlFreeType->font_size;

			area_left = std::min<int>(area_left, o.area_left);
			area_right = std::max<int>(area_right, o.area_right);
			area_bottom = std::min<int>(area_bottom, o.area_bottom);
			area_top = std::max<int>(area_top, o.area_top);

			// next line
			i_pos_y -= cGlFreeType->font_size;
		}

		area_width = area_right - area_left;
		area_height = area_top - area_bottom;
	}

	void renderConfigContent()
	{
		if (!visible)
			return;

		GlRenderOStream &cGlRenderOStream_ = *cGlRenderOStream;

		cGlFreeType->setColor(GLSL::vec3(1,1,1));
		CGlErrorCheck();

		for (std::list<COption>::iterator i = option_list.begin(); i != option_list.end(); i++)
		{
			COption &o = *i;

			if (o.button_left_down)
				cGlFreeType->setColor(GLSL::vec3(0.7,1,0.7));
			else if (o.in_area)
				cGlFreeType->setColor(GLSL::vec3(1,0.5,0.5));


			std::array<float,2> tmp;
			tmp = {(float)o.description_render_left, (float)o.description_render_top};
			cGlFreeType->setPosition(tmp);

			cGlRenderOStream_ << o.description << std::flush;

			tmp = {(float)o.value_render_left, (float)o.value_render_top};
			cGlFreeType->setPosition(tmp);

			switch(o.type)
			{
				case CONFIG_BOOLEAN:
					cGlRenderOStream_ << (*o.configBoolean.value ? "X" : "O") << std::flush;
					if (&o == switch_option)
						cGlRenderOStream_ << " S" << std::flush;
					break;

				case CONFIG_INT:
					cGlRenderOStream_ << *o.configInt.value << std::flush;
					break;

				case CONFIG_UINT:
					cGlRenderOStream_ << *o.configUInt.value << std::flush;
					break;

				case CONFIG_LONGLONG:
					cGlRenderOStream_ << *o.configLongLong.value << std::flush;
					break;

				case CONFIG_ULONGLONG:
					cGlRenderOStream_ << *o.configULongLong.value << std::flush;
					break;

				case CONFIG_FLOAT:
					cGlRenderOStream_ << *o.configFloat.value << std::flush;
					break;

				case CONFIG_FLOAT_HI:
					cGlRenderOStream_ << *o.configFloatHI.value << std::flush;
					break;

				case CONFIG_DOUBLE:
					cGlRenderOStream_ << *o.configDouble.value << std::flush;
					break;

				case CONFIG_DOUBLE_HI:
					cGlRenderOStream_ << *o.configDoubleHI.value << std::flush;
					break;

				default:
					break;
			}
			CGlErrorCheck();

			if (o.button_left_down || o.in_area)
				cGlFreeType->setColor(GLSL::vec3(1,1,1));

			CGlErrorCheck();
		}
	}


	int old_mouse_x;
	int old_mouse_y;

	/**
	 * CALLBACK FUNCTION: mouse button is pressed
	 */
	void mouse_button_down(int button)
	{
		if (active_option != nullptr)
		{
			if (button == RenderWindow::MOUSE_BUTTON_LEFT)
				active_option->button_left_down = true;

			if (button == RenderWindow::MOUSE_BUTTON_RIGHT)
				active_option->button_right_down = true;
		}
	}

	/**
	 * CALLBACK FUNCTION: mouse button is released
	 */
	void mouse_wheel(int x, int y)
	{

		if (active_option == nullptr)
			return;

		switch(active_option->type)
		{
			case CONFIG_BOOLEAN:
				if (y & 1)
				{
					active_option->configBoolean.change();
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				}
				break;

			case CONFIG_INT:
				if (y > 0)
					active_option->configInt.up();
				else if (y < 0)
					active_option->configInt.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;

			case CONFIG_UINT:
				if (y > 0)
					active_option->configUInt.up();
				else if (y < 0)
					active_option->configUInt.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;

			case CONFIG_LONGLONG:
				if (y > 0)
					active_option->configLongLong.up();
				else if (y < 0)
					active_option->configLongLong.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;

			case CONFIG_ULONGLONG:
				if (y > 0)
					active_option->configULongLong.up();
				else if (y < 0)
					active_option->configULongLong.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;

			case CONFIG_FLOAT:
				if (y > 0)
					active_option->configFloat.up();
				else if (y < 0)
					active_option->configFloat.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;


			case CONFIG_FLOAT_HI:
				if (y > 0)
					active_option->configFloatHI.up();
				else if (y < 0)
					active_option->configFloatHI.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;



			case CONFIG_DOUBLE:
				if (y > 0)
					active_option->configDouble.up();
				else if (y < 0)
					active_option->configDouble.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;


			case CONFIG_DOUBLE_HI:
				if (y > 0)
					active_option->configDoubleHI.up();
				else if (y < 0)
					active_option->configDoubleHI.down();

				if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
				break;

			default:
				break;
		}
	}


	/**
	 * CALLBACK FUNCTION: mouse button is released
	 */
	void mouse_button_up(int button)
	{
		if (active_option == nullptr)
			return;

		mouse_motion(old_mouse_x, old_mouse_y);

		if (active_option->button_left_down == true)
		{
			active_option->button_left_down = false;

			switch(active_option->type)
			{
				case CONFIG_BOOLEAN:
					active_option->configBoolean.change();

					if (switch_option != nullptr && switch_option != active_option)
						switch_option->configBoolean.change();
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_INT:
					break;

				case CONFIG_UINT:
					break;

				case CONFIG_LONGLONG:
					break;

				case CONFIG_ULONGLONG:
					break;

				case CONFIG_FLOAT:
					break;

				default:
					break;
			}
			return;
		}

		if (active_option->button_right_down == true)
		{
			active_option->button_right_down = false;

			if (switch_option == active_option)
			{
				// deactivate this option as the switch option button
				switch_option = nullptr;
				return;
			}
			else
			{
				// otherwise set the current option to the switch option
				if (active_option->type != CONFIG_BOOLEAN)
					return;

				switch_option = active_option;
				return;
			}
		}
	}

	/**
	 * CALLBACK FUNCTION: mouse is moved
	 *
	 * \return true, if the mouse motion is handled by the configuration hud
	 */
	bool mouse_motion(int x, int y)
	{
		if (active_option != nullptr)
		{
			int dx = (old_mouse_x - x) + (old_mouse_y - y);
//			int dy = old_mouse_y - y;

			if (active_option->button_left_down)
			{
				switch(active_option->type)
				{
				case CONFIG_FLOAT:
					if (dx > 0)			active_option->configFloat.down(dx);
					else if (dx < 0)	active_option->configFloat.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_FLOAT_HI:
					if (dx > 0)			active_option->configFloatHI.down(dx);
					else if (dx < 0)	active_option->configFloatHI.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;


				case CONFIG_DOUBLE:
					if (dx > 0)			active_option->configDouble.down(dx);
					else if (dx < 0)	active_option->configDouble.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_DOUBLE_HI:
					if (dx > 0)			active_option->configDoubleHI.down(dx);
					else if (dx < 0)	active_option->configDoubleHI.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;


				case CONFIG_INT:
					if (dx > 0)			active_option->configInt.down(dx);
					else if (dx < 0)	active_option->configInt.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_UINT:
					if (dx > 0)			active_option->configUInt.down(dx);
					else if (dx < 0)	active_option->configUInt.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_LONGLONG:
					if (dx > 0)			active_option->configLongLong.down(dx);
					else if (dx < 0)	active_option->configLongLong.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				case CONFIG_ULONGLONG:
					if (dx > 0)			active_option->configULongLong.down(dx);
					else if (dx < 0)	active_option->configULongLong.up(-dx);
					if (active_option->callback_changed)	active_option->callback_changed(active_option->callback_user_ptr);
					break;

				default:
					break;
				}
				old_mouse_x = x;
				old_mouse_y = y;
				return true;
			}
		}

		old_mouse_x = x;
		old_mouse_y = y;

		/**
		 * we only like to change the state of an option, when we are inside the configuration area
		 */
		if (area_left <= x && area_right >= x && area_top >= y && area_bottom <= y)
		{
			int c = 0;
			for (std::list<COption>::iterator i = option_list.begin(); i != option_list.end(); i++)
			{
				COption &o = *i;

				// check, if we are inside the current option area
				if (o.area_left <= x && o.area_right > x && o.area_top > y && o.area_bottom <= y)
				{
					// we are still in in the same valid area => do nothing
					if (active_option == &o)
						return true;

					if (active_option != nullptr)
					{
						active_option->in_area = false;
						active_option->button_left_down = false;
					}

					if (o.type == CONFIG_TEXT || o.type == CONFIG_LINEBREAK)
					{
						active_option = nullptr;
						return true;
					}

					// activate area
					o.in_area = true;
					active_option = &o;
					return true;
				}
				c++;
			}
		}

		if (active_option == nullptr)
			return true;

		active_option->in_area = false;
		active_option->button_left_down = false;
		active_option = nullptr;

		return false;
	}
};

#endif /* CGLHUDCONFIG_HPP_ */
