#version 150

in vec4 vertex_position;
in vec3 vertex_normal3;

uniform mat4 pvm_matrix;
uniform mat3 view_model_normal_matrix3;	// transposed, inverse of pvm matrix
uniform mat4 view_model_matrix;

out vec3 Normal3;
out vec3 frag_position3;

void main(void)
{
	gl_Position = pvm_matrix*vertex_position;

	Normal3 = view_model_normal_matrix3*vertex_normal3;

	frag_position3 = vec3(view_model_matrix*vertex_position);
}
