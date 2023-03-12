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
 * GlCommonShaderPrograms.hpp
 *
 *  Created on: Mar 22, 2010
 * Author: Martin SCHREIBER <schreiberx@gmail.com> (schreiberx@gmail.com)
 */

#ifndef CCOMMON_SHADERPROGRAMS_HPP_
#define CCOMMON_SHADERPROGRAMS_HPP_

#include <sweet/libgl/shaders/shader_blinn/CShaderBlinn.hpp>
#include <sweet/libgl/shaders/shader_blinn_shadow_map/CShaderBlinnShadowMap.hpp>
#include <sweet/libgl/shaders/shader_cube_map_mirror/CShaderCubeMapMirror.hpp>
#include <sweet/libgl/shaders/shader_texturize/CShaderTexturize.hpp>
#include <sweet/libgl/shaders/shader_texturize_rainbowmap/CShaderTexturizeRainbowmap.hpp>
//#include <sweet/libgl/shaders/shader_color/CShaderColor.hpp>
#include <sweet/libgl/shaders/shader_height_color_blinn/CShaderHeightColorBlinn.hpp>



/**
 * \brief class to handle different shaders stored in shaders
 */
class GlCommonShaderPrograms
{
//	bool verbose;

public:
	GlShaderCubeMapMirror		shaderCubeMapMirror;		///< shader to render cubemap
	GlShaderTexturize			shaderTexturize;			///< shader for simple texturization
	GlShaderTexturizeRainbowmap	shaderTexturizeRainbowmap;	///< shader for simple texturization using a rainbow map
	GlShaderBlinn				shaderBlinn;				///< blinn shader
	GlShaderBlinnShadowMap		shaderBlinnShadowMap;		///< blinn shader with shadow map
//	GlShaderColor				shaderColor;				///< color shader

	GlCommonShaderPrograms(
			bool i_verbose = true
	)//	:
//			verbose(i_verbose)
	{
	}
};

#endif
