#version 150

/**
 * blinn shader
 */


/*
 * input from vertex shader
 */
in vec3 frag_normal3;
in vec3 frag_position3;

// LIGHT0
uniform bool light0_enabled = true;
uniform vec3 light0_view_pos3 = vec3(10,10,10);	// position in view space

uniform float light0_ambient_intensity = 0.2f;
uniform vec3 light0_ambient_color3 = vec3(1.f,1.f,1.f);

uniform float light0_diffuse_intensity = 0.5f;
uniform vec3 light0_diffuse_color3 = vec3(1.f,1.f,1.f);

uniform float light0_specular_intensity = 0.8f;
uniform float light0_specular_exponent = 20.0f;
uniform vec3 light0_specular_color3 = vec3(1.f,1.f,1.f);

vec4 blinnShading(vec4 material_color)
{
	if (!light0_enabled)
		return material_color;

	vec3 frag_normal3_normalized = normalize(frag_normal3);

	// compute vector from surface point to light
	vec3 light_vec = normalize(light0_view_pos3 - frag_position3);

	vec3 viewer_vec = -normalize(frag_position3);
	vec3 halfway_vec = normalize(light_vec+viewer_vec);

	return material_color * vec4(
					(
					light0_ambient_intensity*(light0_ambient_color3)	+
					light0_diffuse_intensity*max(0.0f, dot(light_vec, frag_normal3_normalized))*light0_diffuse_color3
					), 1.0f)
					+
					vec4(	(light0_specular_intensity*
							max(0.0f, pow(dot(halfway_vec,frag_normal3_normalized), light0_specular_exponent))
							*light0_specular_color3), 1.0f);
}
