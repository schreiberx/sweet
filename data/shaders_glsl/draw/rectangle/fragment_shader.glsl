#version 150

in vec2 rast_texture_coord;

out vec4 frag_data;

uniform sampler2DRect slice_rectangle;

void main(void)
{
	frag_data = texture(slice_rectangle, rast_texture_coord);
}
