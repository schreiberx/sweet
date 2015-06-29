#version 150

in vec3 frag_normal;
out vec4 frag_data;

uniform vec3 color = vec3(1,1,0);

void main(void)
{
#if 0
	frag_data = vec4(frag_normal, 0.0f);
#else
	float light = dot(normalize(frag_normal), normalize(vec3(1,1,1)));
	frag_data = vec4(color*light, 0.0f);
#endif
}
