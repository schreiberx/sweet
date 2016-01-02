// use opengl 3.2 core version!
#version 150

in vec4 vertex_attrib;

void main(void)
{
	gl_Position = vertex_attrib;
}
