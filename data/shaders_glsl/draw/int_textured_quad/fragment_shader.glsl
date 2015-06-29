#version 150

in vec2 rast_texture_coord;

out vec4 frag_data;

uniform usampler2D sampler;

void main(void)
{
	frag_data = vec4(texture(sampler, rast_texture_coord));
}
