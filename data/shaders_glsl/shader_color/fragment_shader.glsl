#version 150

// COLOR of vertex if texturing is disabled
uniform vec4 frag_color = vec4(1.f,1.f,1.f,1.f);

out vec4 frag_data;

void main(void)
{
	frag_data = frag_color;
}
