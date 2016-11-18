#version 150

in vec4 vertex_position4;
in vec3 vertex_normal3;
in vec2 vertex_texture_coord2;

out vec3 frag_vertex_normal3;
out vec2 frag_texture_coord2;

uniform mat4 pvm_matrix;
//uniform mat4 normal_matrix;

void main(void)
{
	gl_Position = pvm_matrix*vertex_position4;

	frag_vertex_normal3 = vertex_normal3;
	frag_texture_coord2 = vertex_texture_coord2;
}
