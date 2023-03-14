#ifndef CGL_SHADER_CUBE_MAP_HPP
#define CGL_SHADER_CUBE_MAP_HPP

#include <sweet/libgl/core/GlError.hpp>
#include <sweet/libgl/core/GlTexture.hpp>

/**
 * Shader to implement reflections using cube maps
 */
class GlShaderCubeMapMirror	: public GlProgram
{
public:
	GlUniform pvm_matrix;				///< uniform to pvm matrix
	GlUniform view_model_normal_matrix3;	///< uniform to view-model normal matrix
	GlUniform view_model_matrix;			///< uniform to view-model matrix

//	GlProgram::CUniform view_matrix;
	GlUniform transposed_view_matrix3;	///< uniform to transposed of view matrix

	GlUniform light0_enabled;			///< uniform to enable or disable uniforms
	GlUniform light0_view_pos3;			///< uniform to position of light

	GlUniform light0_ambient_intensity;	///< uniform to ambient intensity
	GlUniform light0_ambient_color3;		///< uniform to ambient color

	GlUniform light0_diffuse_intensity;	///< uniform to diffuse intensity
	GlUniform light0_diffuse_color3;		///< uniform to diffuse color

	GlUniform light0_specular_intensity;	///< uniform to specular intensity
	GlUniform light0_specular_exponent;	///< uniform to specular exponent
	GlUniform light0_specular_color3;	///< uniform to specular color

	GlShaderCubeMapMirror()
	{
		if (!initVertFragShadersFromDirectory("shader_cube_map_mirror"))
			return;

		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(pvm_matrix, "pvm_matrix");
		setupUniform(view_model_normal_matrix3, "view_model_normal_matrix3");
		setupUniform(view_model_matrix, "view_model_matrix");

		setupUniform(transposed_view_matrix3, "transposed_view_matrix3");

		setupUniform(light0_enabled, "light0_enabled");
		setupUniform(light0_view_pos3, "light0_view_pos3");
		setupUniform(light0_ambient_intensity, "light0_ambient_intensity");
		setupUniform(light0_ambient_color3, "light0_ambient_color3");
		setupUniform(light0_diffuse_intensity, "light0_diffuse_intensity");
		setupUniform(light0_diffuse_color3, "light0_diffuse_color3");
		setupUniform(light0_specular_intensity, "light0_specular_intensity");
		setupUniform(light0_specular_exponent, "light0_specular_exponent");
		setupUniform(light0_specular_color3, "light0_specular_color3");
	}

	~GlShaderCubeMapMirror()
	{
	}
};


#endif
