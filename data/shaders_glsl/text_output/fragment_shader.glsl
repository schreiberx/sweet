#version 150

in vec2 rast_texture_coord;
out vec4 frag_data;

uniform sampler2D sampler;

void main(void)
{
	vec4 col = texture(sampler, rast_texture_coord);
	frag_data = col*col.a;
}
