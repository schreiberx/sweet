#ifndef CGL_SHADER_BLINN_HPP
#define CGL_SHADER_BLINN_HPP

#include "libgl/shaders/CDefaultShaderDir.hpp"
#include "libgl/core/CGlTexture.hpp"
#include "libgl/core/CGlError.hpp"


/**
 * general blinn shader to use for rendering vertices
 */
#include "libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp"
#include "libgl/core/CGlProgram.hpp"


class GlShaderBlinn	:
	public CGlProgram,
	public CShaderBlinnSkeleton
{
public:
	CGlUniform texture0_enabled;	///< uniform to enable and disable texturing
	CGlUniform vertex_color;		///< uniform to basic vertex color of fragment

	GlShaderBlinn()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_blinn");
		attachFragShader(SHADER_GLSL_DEFAULT_DIR"shader_blinn/fragment_shader_skeleton.glsl");

		// link programs
		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(texture0_enabled, "texture0_enabled");
		setupUniform(vertex_color, "vertex_color");

		initBlinnSkeleton(*this);
	}

	~GlShaderBlinn()
	{
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(	CGlMaterial	&material,
						Lights &lights,
						const GLSL::vec3 &light_view_pos3
	)
	{
		CShaderBlinnSkeleton::setupUniforms(material, lights, light_view_pos3);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(	CGlMaterial	&material,
						Lights &lights
	)
	{
		CShaderBlinnSkeleton::setupUniforms(material, lights);

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
