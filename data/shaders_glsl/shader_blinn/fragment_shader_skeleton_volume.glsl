#version 150

// LIGHT0
uniform bool light0_enabled = true;
uniform vec3 light0_view_pos3 = vec3(10,10,10);	// position in view space

uniform vec3 light0_ambient_color3 = vec3(1.f,1.f,1.f);
uniform vec3 light0_diffuse_color3 = vec3(1.f,1.f,1.f);
uniform vec3 light0_specular_color3 = vec3(1.f,1.f,1.f);

uniform vec3 material_ambient_color3 = vec3(1.f,1.f,1.f);
uniform vec3 material_diffuse_color3 = vec3(1.f,1.f,1.f);
uniform vec3 material_specular_color3 = vec3(1.f,1.f,1.f);
uniform float material_specular_exponent = 20.0f;

vec4 blinnShading(	vec4 material_color,	///< color of material
					vec3 view_normal3,		///< normal of fragment in viewspace
					vec3 view_position3		///< position of fragment in viewspace
)
{
	if (!light0_enabled)
		return material_color;

	view_normal3 = normalize(view_normal3);

	// compute vector from surface point to light
	vec3 light_vec = normalize(light0_view_pos3 - view_position3);

	vec3 viewer_vec = -normalize(view_position3);
	vec3 halfway_vec = normalize(light_vec+viewer_vec);

	return material_color * vec4(
					(
						light0_ambient_color3*material_ambient_color3	+
						max(0.0f, dot(light_vec, view_normal3))*light0_diffuse_color3*material_diffuse_color3
					), 1.0f)
					+
					vec4(	(pow(max(0.0f, dot(halfway_vec, view_normal3)), material_specular_exponent)
							*light0_specular_color3*material_specular_color3), 1.0f);
}
