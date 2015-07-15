#version 150

in vec2 rast_texture_coord;

out vec3 frag_data;

uniform sampler2D sampler;

/**
 * computation of rainbow color map from
 * http://http.developer.nvidia.com/GPUGems/gpugems_ch08.html
 */
vec3 blend2ndOrder(vec3 f)
{
	vec3 y = 1.0-f*f;
	return max(y, vec3(0));
}


void main(void)
{
	float red_value = texture(sampler, rast_texture_coord).r;
//	red_value = rast_texture_coord.x;

	float limited_height = max(0.0, min(1.0, red_value));
#if 1
	frag_data.xyz = blend2ndOrder(
					vec3(3.3)*(limited_height-vec3(0.75, 0.5, 0.25))
				);
#else
	frag_data.xyz = blend2ndOrder(
					vec3(3.0)*(limited_height-vec3(1.0, 0.5, 0.0))
				);
#endif
	return;
}
