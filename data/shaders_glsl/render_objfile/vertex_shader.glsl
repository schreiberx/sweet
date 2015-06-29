#version 150

in vec4 vertex_position;
in vec3 vertex_normal;

uniform mat4 pvm_matrix;
uniform mat3 normal_matrix;

out vec3 frag_normal;

void main(void)
{
	gl_Position = pvm_matrix*vertex_position;
	frag_normal = vec3(normal_matrix*vertex_normal);
}
