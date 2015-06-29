#version 150

// TEXTURE0
uniform bool texture0_enabled = false;
uniform sampler2D texture0;
uniform sampler2DShadow texture_shadow_map;
uniform mat4 shadow_map_matrix;

in vec3 frag_position3;

in vec2 frag_texture_coord2;

out vec4 frag_data;

vec4 blinnShading(vec4 material_color);
uniform vec3 material_ambient_color3 = vec3(1.f,1.f,1.f);

void main(void)
{
	vec4 frag_color;
	if (texture0_enabled)	frag_color = texture(texture0, frag_texture_coord2);
	else					frag_color = vec4(1);

	frag_data = blinnShading(frag_color);

	vec4 p = shadow_map_matrix*vec4(frag_position3, 1);
	p.xyz /= p.w;

	if (p.w > 0)
		if (texture(texture_shadow_map, p.xyz) > gl_FragCoord.z)
			frag_data.xyz = material_ambient_color3*frag_color.xyz;
//	frag_data.xyz = vec3(gl_FragCoord.z);
}
