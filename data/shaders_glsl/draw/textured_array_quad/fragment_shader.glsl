#version 150

in vec2 rast_texture_coord;

out vec4 frag_data;

uniform sampler2DArray sampler;
uniform int texture_array_layer;

void main(void)
{
	frag_data = texture(sampler, vec3(rast_texture_coord, texture_array_layer));
}
