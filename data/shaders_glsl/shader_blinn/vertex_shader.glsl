#version 150

in vec4 vertex_position;
in vec3 vertex_normal3;
in vec2 vertex_texture_coord2;

uniform mat4 pvm_matrix;
uniform mat3 view_model_normal_matrix3;	// transposed, inverse of pvm matrix
uniform mat4 view_model_matrix;

out vec3 frag_normal3;
out vec2 frag_texture_coord2;
out vec3 frag_position3;

out vec4 screen_position;

void main(void)
{
	gl_Position = pvm_matrix*vertex_position;

	frag_normal3 = view_model_normal_matrix3*vertex_normal3;
	frag_texture_coord2 = vertex_texture_coord2;
	frag_position3 = vec3(view_model_matrix*vertex_position);
}
