#ifndef CGL_SHADER_HEIGHT_COLOR_BLINN_HPP
#define CGL_SHADER_HEIGHT_COLOR_BLINN_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlProgram.hpp>
#include <libgl/core/GlTexture.hpp>
#include "libgl/shaders/CDefaultShaderDir.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp"


class GLShaderHeightColorBlinn	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	///< uniform to enable and disable texturing
	GlUniform vertex_color;		///< uniform to basic vertex color of fragment
	GlUniform height_color_scale;		///< uniform to setup scaling factor of height
	GlUniform height_color_offset;		///< uniform to setup offset of height applied before scaling

	GLShaderHeightColorBlinn()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_height_color_blinn");
		attachFragShader(SHADER_GLSL_DEFAULT_DIR"shader_blinn/fragment_shader_skeleton.glsl");

		// link programs
		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(height_color_scale, "height_color_scale");
		setupUniform(height_color_offset, "height_color_offset");

		initBlinnSkeleton(*this);

		setupColorScaleAndOffset();
	}


	/**
	 * setup the offset value and scaling factor for the height as input
	 * value for the color map
	 */
	void setupColorScaleAndOffset(
			float offset = 0.0,
			float scale = 1.0
		)
	{

		use();
		height_color_scale.set1f(scale);
		height_color_offset.set1f(offset);
		disable();
	}


	~GLShaderHeightColorBlinn()
	{
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(
			GlMaterial	&material,
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
	void setupUniformsMaterial(
			GlMaterial	&material
	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};


#endif
