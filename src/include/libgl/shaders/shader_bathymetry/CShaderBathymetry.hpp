#ifndef CGL_SHADER_BATHYMETRY_HPP
#define CGL_SHADER_BATHYMETRY_HPP

#include "libgl/shaders/CDefaultShaderDir.hpp"
#include "libgl/core/GlTexture.hpp"
#include "libgl/core/CGlError.hpp"


/**
 * general blinn shader to use for rendering vertices
 */
#include "libgl/core/GlProgram.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp"


class CShaderBathymetry	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	///< uniform to enable and disable texturing

	CShaderBathymetry()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_bathymetry");
		attachFragShader(SHADER_GLSL_DEFAULT_DIR"shader_blinn/fragment_shader_skeleton.glsl");

		// link programs
		link();
		if (error())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		initBlinnSkeleton(*this);
	}

	~CShaderBathymetry()
	{
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(
			CGlMaterial	&material,
			CGlLights &lights,
			const GLSL::vec3 &light_view_pos3
	)
	{
		CShaderBlinnSkeleton::setupUniforms(material, lights, light_view_pos3);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniformsMaterial(
			CGlMaterial	&material
	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};


#endif
