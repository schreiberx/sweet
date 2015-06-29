#version 330

out vec4 frag_data;
in vec2 gl_PointCoord;
in float surface_height;

void main(void)
{
	float a = length(gl_PointCoord.xy-vec2(0.5))*2.0;
	if (a >= 1.0)	discard;
	a = 1.0-a;

/*
 * multiplier can be used to use texture format with reduced accuracy with range [0;256] 
 * to increase final accuracy since all a values above 1 will be padded to 1
 */
#define MULTIPLIER	0.25
	frag_data = vec4(surface_height*a*MULTIPLIER, a*MULTIPLIER, 0, 0);
}
