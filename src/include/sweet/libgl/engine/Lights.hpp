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
 * CGlLights.hpp
 *
 *  Created on: Mar 27, 2010
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <array>

#ifndef CLIGHTS_HPP_
#define CLIGHTS_HPP_

class Lights
{
public:
	class Light
	{
	public:
		bool enabled;

		std::array<float,3> light_view_pos3;		///< position of light in view space

		std::array<float,3> world_pos3;		///< position of light in world space
		std::array<float,3> ambient_color3;	///< ambient color
		std::array<float,3> diffuse_color3;	///< diffuse color
		std::array<float,3> specular_color3;	///< specular color
	};

	std::array<Light,1> lights;

	Lights()
	{
		lights[0].world_pos3 = {{2, 3, 8}};

		lights[0].ambient_color3 = {{0.8, 0.8, 0.8}};
		lights[0].diffuse_color3 = {{0.8, 0.8, 0.8}};
		lights[0].specular_color3 = {{1, 1, 1}};
	}

};

#endif /* CGLLIGHT_HPP_ */
