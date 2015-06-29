#version 150

in vec4 vertex_position;
in vec3 vertex_texture_coord;
smooth out vec4 rast_texture_coord;

uniform mat4 pvm_matrix;

void main(void)
{
	gl_Position = pvm_matrix*vertex_position;

	// avoid clipping!
//	float positive = float(gl_Position.z < 0.0);
//	gl_Position.z = gl_Position.w*positive + (1.0f-positive)*gl_Position.z;

	rast_texture_coord = vec4(vertex_texture_coord, 1.0);
}
