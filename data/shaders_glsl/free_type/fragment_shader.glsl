// use opengl 3.2 core shaders
#version 150

in vec2 texel_pos;

out vec4 frag_data;

uniform sampler2D sampler;
uniform vec3 color;

void main(void)
{
	float a = texelFetch(sampler, ivec2(texel_pos), 0).r;

	if (a <= 0.0)
		discard;

	frag_data = vec4(color,a);
}
