#version 150

in vec3 rast_texture_coord;
out vec4 frag_data;
uniform samplerCube sampler;

void main(void)
{
	frag_data = texture(sampler, rast_texture_coord);
}
