#version 150

in vec3 rast_texture_coord;
out vec4 frag_data;

uniform samplerCube sampler_cube;

void main(void)
{
	frag_data = texture(sampler_cube, rast_texture_coord);
//	frag_data = texture(sampler_cube, vec3(1,0,0));
}
