#version 150


in vec3 Normal3;
in vec2 TexCoord2;
in vec3 frag_position3;

uniform samplerCube cube_map;

uniform vec4 light_pos = vec4(1,1,1,0);
uniform vec3 light_color3 = vec3(1,1,1);

out vec4 frag_data;
//uniform mat3 normal_matrix3;
uniform mat3 transposed_view_matrix3;
//uniform mat4 view_matrix;


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

void main(void)
{
	vec3 frag_normal3_normalized = normalize(Normal3);
	vec3 reflected_vector3 = reflect(normalize(frag_position3), frag_normal3_normalized);

	// we use the transpose cz. we invert the normal which is invert(transp(invert(M))) = transp(M)
	vec4 frag_color = texture(cube_map, transposed_view_matrix3*reflected_vector3);

	if (light0_enabled)
	{
		// compute vector from surface point to light
		vec3 light_vec = normalize(light0_view_pos3 - frag_position3);

		vec3 viewer_vec = -normalize(frag_position3);
		vec3 halfway_vec = normalize(light_vec+viewer_vec);

		frag_color = frag_color * vec4(
						(
						light0_ambient_intensity*(light0_ambient_color3)	+
						light0_diffuse_intensity*max(0.0f, dot(light_vec, frag_normal3_normalized))*light0_diffuse_color3
						), 1.0f)
						+
						vec4(	(light0_specular_intensity*
								max(0.0f, pow(dot(halfway_vec,frag_normal3_normalized), light0_specular_exponent))
								*light0_specular_color3), 1.0f);
	}
	frag_data = frag_color;
}
