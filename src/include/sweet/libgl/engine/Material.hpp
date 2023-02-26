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
 * CGlMaterial.hpp
 *
 *  Created on: Mar 27, 2010
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef CGLMATERIAL_HPP_
#define CGLMATERIAL_HPP_

#include <sweet/libgl/core/GlTexture.hpp>
#include <sweet/libmath/CGlSlMath.hpp>

class GlMaterial
{
public:
	std::string name;			///< name for material (just for convenience)

	GLSL::vec3 ambient_color3;	///< ambient color
	GLSL::vec3 diffuse_color3;	///< diffuse color
	float specular_exponent;	///< specular exponent
	GLSL::vec3 specular_color3;	///< specular color

	GlTexture *texture0;
	GlTexture *normal0;

	GlMaterial()	:
		ambient_color3(0.2,0.2,0.2),
		diffuse_color3(0.6,0.6,0.6),
		specular_exponent(20),
		specular_color3(1.0,1.0,1.0),
		texture0(nullptr),
		normal0(nullptr)
	{
	}

	void loadTexture0(const std::string &file)
	{
		if (texture0)
			delete texture0;

		texture0 = new GlTexture;
		texture0->loadFromFile(file);
	}


	void loadNormal0(const std::string &file)
	{
		if (normal0)
			delete normal0;

		normal0 = new GlTexture;
		normal0->loadFromFile(file);
	}

	void cleanup()
	{
		if (texture0)
		{
			delete texture0;
			texture0 = nullptr;
		}

		if (normal0)
		{
			delete normal0;
			normal0 = nullptr;
		}
	}

	~GlMaterial()
	{
		cleanup();
	}
};

#endif
