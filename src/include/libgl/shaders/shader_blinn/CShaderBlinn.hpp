#ifndef CGL_SHADER_BLINN_HPP
#define CGL_SHADER_BLINN_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlProgram.hpp>
#include "libgl/shaders/CDefaultShaderDir.hpp"
//#include "libgl/core/GlTexture.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp"


class GlShaderBlinn	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	///< uniform to enable and disable texturing
	GlUniform vertex_color;		///< uniform to basic vertex color of fragment

	GlShaderBlinn()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_blinn");
		attachFragShader((ShaderDir::getDirectory()+"shader_blinn/fragment_shader_skeleton.glsl").c_str());

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
	void setupUniforms(	GlMaterial	&material,
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
	void setupUniforms(	GlMaterial	&material,
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
			GlMaterial	&material
	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};


#endif
