#version 150

out vec4 frag_data;

uniform vec4 color = vec4(1,0,0,0);

void main(void)
{
	frag_data = color;
}
