#ifndef CGL_SHADER_BLINN_SHADOW_AND_CAUSTIC_MAP_HPP
#define CGL_SHADER_BLINN_SHADOW_AND_CAUSTIC_MAP_HPP

#include <sweet/libgl/core/GlTexture.hpp>
#include <sweet/libgl/core/CGlError.hpp>


/**
 * blinn shader supporting shadow maps
 */
#include <sweet/libgl/core/GlProgram.hpp>
#include "shaders/shader_blinn/CShaderBlinnSkeleton.hpp"


class CShaderBlinnShadowAndCausticMap	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	///< uniform to enable and disable texturing
	GlUniform shadow_map_matrix_uniform;	///< shadow map matrix


	CShaderBlinnShadowAndCausticMap()
	{
		initVertFragShadersFromDirectory("shader_blinn_shadow_and_caustic_map");
		CError_PtrAppendReturn(this);

		// link programs
		link();
		if (error())
		{
			std::cerr << "info Log: during linking: " << getInfoLog() << std::endl;
			return;
		}

		setupUniform(texture0_enabled, "texture0_enabled");
		setupUniform(shadow_map_matrix_uniform, "shadow_map_matrix");

		use();
		setUniform1i("texture0", 0);
		setUniform1i("texture_shadow_map", 1);
		setUniform1i("texture_caustic_map", 2);
		setUniform1i("texture_caustic_depth_map", 3);
		disable();

		initBlinnSkeleton(*this);
	}

	~CShaderBlinnShadowAndCausticMap()
	{
	}

	/**
	 * setup the light for rendering
	 */
	void setupUniforms(	GlMaterial	&material,
						CGlLights &lights,
						const GLSL::vec3 &light_view_pos3,
						const GLSL::mat4 &shadow_map_matrix
	)
	{
		CShaderBlinnSkeleton::setupUniformsLights(lights, light_view_pos3);
		setupUniformsMaterial(material);
		setupUniformsShadowMapping(shadow_map_matrix);
	}

	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsShadowMapping(	const GLSL::mat4 &shadow_map_matrix	)
	{

		shadow_map_matrix_uniform.set(shadow_map_matrix);
	}

	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsMaterial(	GlMaterial	&material	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};

#endif
