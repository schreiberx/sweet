#version 150

// TEXTURE0
uniform bool texture0_enabled = false;
uniform sampler2D texture0;

in vec2 frag_texture_coord2;

out vec4 frag_data;

vec4 blinnShading(vec4 material_color);

void main(void)
{
	vec4 frag_color;
	if (texture0_enabled)	frag_color = texture(texture0, frag_texture_coord2);
	else					frag_color = vec4(1);

	frag_data = blinnShading(frag_color);
}
