#version 150
/*
 * input from vertex shader
 */
in vec3 frag_normal3;
in vec3 frag_position3;

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


vec4 blinnShading(vec4 material_color)
{
	if (!light0_enabled)
		return vec4(0,0,0,0);

	vec3 frag_normal3_normalized = normalize(frag_normal3);

	// compute vector from surface point to light
	vec3 light_vec = normalize(light0_view_pos3 - frag_position3);

	vec3 viewer_vec = -normalize(frag_position3);
	vec3 halfway_vec = normalize(light_vec+viewer_vec);

#if 0
	// distance attenuation
	float dist = length(light0_view_pos3 - frag_position3);

	float attenuation = min(1.0/(dist*dist*0.005) + 1.0/(dist*0.05), 1.0);

//return vec4(material_ambient_color3);
	return 		(	material_color*vec4(
						(
						light0_ambient_color3*material_ambient_color3	+
						max(0.0f, dot(light_vec, frag_normal3_normalized))*light0_diffuse_color3*material_diffuse_color3
						), 1.0f)
					+
					vec4(	(pow(max(0.0f, dot(halfway_vec, frag_normal3_normalized)), material_specular_exponent)
							*light0_specular_color3*material_specular_color3), 1.0f)
				)*attenuation;
#else
	return 		material_color*vec4(
					(
					light0_ambient_color3*material_ambient_color3	+
					max(0.0f, dot(light_vec, frag_normal3_normalized))*light0_diffuse_color3*material_diffuse_color3
					), 1.0f)
				+
				vec4(	(pow(max(0.0f, dot(halfway_vec, frag_normal3_normalized)), material_specular_exponent)
						*light0_specular_color3*material_specular_color3), 1.0f);
#endif
}
