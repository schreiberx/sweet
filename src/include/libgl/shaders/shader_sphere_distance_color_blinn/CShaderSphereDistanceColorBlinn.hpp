#ifndef CGL_SHADER_SPHERE_DISTANCE_COLOR_BLINN_HPP
#define CGL_SHADER_SPHERE_DISTANCE_COLOR_BLINN_HPP

#include "libgl/shaders/CDefaultShaderDir.hpp"
#include "libgl/core/GlTexture.hpp"
#include "libgl/core/CGlError.hpp"


/**
 * general blinn shader to use for rendering vertices
 */
#include "libgl/core/GlProgram.hpp"
#include "libgl/shaders/shader_blinn/CShaderBlinnSkeleton.hpp"


class CShaderSphereDistanceColorBlinn	:
	public GlProgram,
	public CShaderBlinnSkeleton
{
public:
	GlUniform texture0_enabled;	///< uniform to enable and disable texturing
	GlUniform vertex_color;		///< uniform to basic vertex color of fragment
	GlUniform height_color_scale;		///< uniform to setup scaling factor of height
	GlUniform height_color_offset;		///< uniform to setup offset of height applied before scaling

	CShaderSphereDistanceColorBlinn()
	{
		std::string infoLog;

		initVertFragShadersFromDirectory("shader_sphere_distance_color_blinn");
		attachFragShader(ShaderDir::getDirectory()+"shader_blinn/fragment_shader_skeleton.glsl");

		// link programs
		link();
		if (error())
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
			float i_offset = 1.0,
			float i_scale = 1.0
		)
	{
		use();
		height_color_offset.set1f(i_offset);
		height_color_scale.set1f(i_scale);
		disable();
	}


	virtual ~CShaderSphereDistanceColorBlinn()
	{
	}


	/**
	 * setup the uniforms for rendering
	 */
	void setupUniforms(
			GlMaterial	&material,
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
			GlMaterial	&material
	)
	{
		CShaderBlinnSkeleton::setupUniformsMaterial(material);

		texture0_enabled.set1b(material.texture0 != nullptr);
	}
};


#endif
