#version 150

in vec4 vertex_position;
in vec2 vertex_texture_coord;
out vec2 rast_texture_coord;

//uniform mat4 pvm_matrix;

void main(void)
{
//	gl_Position = pvm_matrix*vertex_position;
	gl_Position = vertex_position;
	rast_texture_coord = vertex_texture_coord;
}
