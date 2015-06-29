#version 150

in vec3 vertex_position;

out vec3 rast_texture_coord;

uniform mat4 pvm_matrix;

void main(void)
{
	gl_Position = pvm_matrix*vec4(vertex_position, 1);
	rast_texture_coord = vertex_position;
}
