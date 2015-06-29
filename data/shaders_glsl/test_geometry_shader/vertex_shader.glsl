#version 150

//in vec4 frag_pos;
uniform mat4 pvm_matrix;

void main(void)
{
	gl_Position = pvm_matrix*vec4(-1,1,0,1);
	return;
	if (gl_VertexID == 0)
		gl_Position = pvm_matrix*vec4(-1,1,0,1);
	else
		gl_Position = pvm_matrix*vec4(1,1,0,1);
//	gl_Position = pvm_matrix*frag_pos;
}

